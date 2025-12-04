# ---------------------------------------------------------------------------
# #         ****        0. Bioconductor or CRAN installs      ****
# ---------------------------------------------------------------------------

# Make sure BiocManager is installed

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

#  Packages you need (CRAN + Bioconductor)

packages <- c(
  # Bioconductor
  "DiffBind", "ComplexHeatmap", "InteractiveComplexHeatmap",
  "ChIPseeker", "AnnotationDbi", "org.Hs.eg.db",
  "GenomeInfoDb", "GenomicRanges", "GenomicAlignments",
  "rtracklayer", "biomaRt", "BSgenome.Hsapiens.UCSC.hg38",
  "Rsubread", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  
  # CRAN
  "DescTools", "tidyverse", "pbapply", "pbmcapply",
  "profileplyr", "vroom"
)

#  Install any missing packages

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Decide which repo to use
    if (pkg %in% rownames(installed.packages())) next   # already there (just in case)
    if (pkg %in% BiocManager::available()) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}


# ---------------------------------------------------------------------------
#                  ****       1. LIBRARY IMPORTS        ****
# ---------------------------------------------------------------------------
# Load Standard Libraries

libraries_to_load <- c(
  "DiffBind", "DescTools", "tidyverse", "ComplexHeatmap",
  "InteractiveComplexHeatmap", "pbapply", "pbmcapply",
  "profileplyr", "ChIPseeker", "AnnotationDbi", "org.Hs.eg.db",
  "GenomeInfoDb", "GenomicRanges", "GenomicAlignments",
  "rtracklayer", "biomaRt", "BSgenome.Hsapiens.UCSC.hg38",
  "Rsubread", "TxDb.Hsapiens.UCSC.hg38.knownGene", "vroom"
)
lapply(libraries_to_load, library, character.only = TRUE)

# --------------------------------------
# 2. Read CSVs with vroom (faster than read.csv)
# --------------------------------------
dexCsv_macs2  <- vroom('csv/macs2/dex_csv_macs2.csv')
dexCsv_seacr_v1 <- vroom('csv/seacr_v1/dex_csv_seacr1.csv')
dexCsv_seacr_v2 <- vroom('csv/seacr_v2/dex_csv_seacr2.csv')
# ---------------------------------------------------------------------------
#                    ****      3.  HOUSEKEEPING        ****
# ---------------------------------------------------------------------------
# Grab names of samples from CSV files
seacr1_names <- dexCsv_seacr_v1$SampleID
seacr2_names <- dexCsv_seacr_v2$SampleID
macs2_names <- dexCsv_macs2$SampleID

#Load genomic databases/references
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

pseudocount <- 1

# ---------------------------------------------------------------------------
#                     ****      4.  FUNCTIONS           ****
# ---------------------------------------------------------------------------
union_matrices <- function(...) {
  mats   <- list(...)                 # all supplied matrices
  if (length(mats) == 0) stop("No matrices supplied")
  # 1. Find the set of all row names that appear in any matrix
  all_rn <- unique(unlist(lapply(mats, rownames)))  # vector of unique row names
  # 2. Create an empty matrix that can hold all data
  n_rows <- length(all_rn)
  n_cols <- ncol(mats[[1]])          # all matrices are assumed to have the same number of columns
  
  combined <- matrix(NA_real_, nrow = n_rows, ncol = n_cols)
  rownames(combined) <- all_rn
  colnames(combined) <- colnames(mats[[1]])
  # 3. Fill the combined matrix with the data from each input matrix
  for (m in mats) {
    # Find where each row from the current matrix fits into the combined one
    idx <- match(rownames(m), all_rn)   # returns indices (or NA if not found)
    
    # Put the values in – NA rows are left untouched
    combined[idx, ] <- m
  }
  
  return(combined)
}

avg_replicates <- function(df) {
      df$A1 = (df$A1_R1 + df$A1_R2) / 2
      df$A2 = (df$A2_R1 + df$A2_R2) / 2
      df$A3 = (df$A3_R1 + df$A3_R2) / 2
      df$A4 = (df$A4_R1 + df$A4_R2) / 2
      df$A5 = (df$A5_R1 + df$A5_R2) / 2
      df$A6 = (df$A6_R1 + df$A6_R2) / 2
      df$A7 = (df$A7_R1 + df$A7_R2) / 2
      df$A8 = (df$A8_R1 + df$A8_R2) / 2
      df2 <- df[, c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8")]
      return(df2)
}

# Function to perform the column subtraction and renaming for a single matrix
process_matrix <- function(matrix_name) {
  # Perform column subtraction
  col1 <- matrix_name[, "A1"] - matrix_name[, "A5"]
  col2 <- matrix_name[, "A2"] - matrix_name[, "A6"]
  col3 <- matrix_name[, "A3"] - matrix_name[, "A7"]
  col4 <- matrix_name[, "A4"] - matrix_name[, "A8"]
  
  # Create a new matrix (or data frame) from the subtracted columns
  new_matrix <- cbind(col1, col2, col3, col4) 
  # Rename the columns
  colnames(new_matrix) <- c("55-3 FC", "55-7 FC", "57-2 FC", "Par FC")
  
  return(new_matrix)
}



PrepareAnnotationGrangesObject <- function(TxDb, map_to_TSS_only = TRUE, map_to_genes = TRUE) {
  
  if (map_to_genes) {
    GR <- genes(TxDb, single.strand.genes.only = TRUE)
    GR <- keepSeqlevels(GR, standardChromosomes(GR)[1:23], pruning.mode = "coarse") %>%
      sortSeqlevels() %>% 
      sort(ignore.strand = TRUE) %>%
      as.data.frame()
    
    GR$ENSEMBL <- suppressMessages(mapIds(org.Hs.eg.db, keys = GR$gene_id, keytype="ENTREZID", column = "ENSEMBL"))
    GR <- subset(GR, !is.na(ENSEMBL))
    GR$SYMBOL <- suppressMessages(mapIds(org.Hs.eg.db, keys = GR$gene_id, keytype="ENTREZID", column = "SYMBOL"))
    GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = TRUE)
    
  } else {
    GR <- transcripts(TxDb)
    GR <- keepSeqlevels(GR, standardChromosomes(GR)[1:23], pruning.mode = "coarse") %>%
      sortSeqlevels() %>% 
      sort(ignore.strand = TRUE) %>%
      as.data.frame()
    
    GR$tx_name <- gsub("\\..*","", GR$tx_name)
    GR$ENSEMBL <- suppressMessages(mapIds(org.Hs.eg.db, keys = GR$tx_name, keytype="ENSEMBLTRANS", column = "ENSEMBL"))
    GR <- subset(GR, !is.na(ENSEMBL))
    GR$SYMBOL <- suppressMessages(mapIds(org.Hs.eg.db, keys = GR$tx_name, keytype="ENSEMBLTRANS", column = "SYMBOL"))
    GR <- makeGRangesFromDataFrame(GR, keep.extra.columns = TRUE)
    
  }
  
  if (map_to_TSS_only) {
    GR <- resize(GR, 1)
  }
  
  return(GR)
  
}


makeMatrixFromDFList <- function(df_list, col_names){
  processed_dfs <-lapply(seq_along(df_list), function(i){
    df <- df_list[[i]][, c("id", "Score")]   # keep only id and Score
    rownames(df) <- df$id
    df$id <- NULL
    setNames(df, paste0("Score", i))
  })
  
  single_matrix <- as.matrix(do.call(cbind, processed_dfs))
  colnames(single_matrix) <- col_names
  return(single_matrix)
}



GetFlankingGenes <- function(query, ref) {
  
  #Get upstream nearest gene
  upstreamGenes <- GenomicRanges::follow(query, unstrand(ref), ignore.strand = TRUE)
  upstreamGenes2 <- upstreamGenes[!is.na(upstreamGenes)]
  q2 <- query[which(!is.na(upstreamGenes))]
  upstreamDist <- GenomicRanges::distance(q2, GR[upstreamGenes2])
  
  #Get downstream nearest gene
  downstreamGenes <- GenomicRanges::precede(query, unstrand(ref), ignore.strand = TRUE)
  downstreamGenes2 <- downstreamGenes[!is.na(downstreamGenes)]
  q2 <- query[which(!is.na(downstreamGenes))]
  downstreamDist <- GenomicRanges::distance(q2, GR[downstreamGenes2])
  
  nearestGenes <- GenomicRanges::nearest(query, unstrand(ref), ignore.strand = TRUE)
  nearestGenes2 <- nearestGenes[!is.na(nearestGenes)]
  q2 <- query[which(!is.na(nearestGenes))]
  nearestDist <- GenomicRanges::distance(q2, GR[nearestGenes2])
  
  query <- as.data.frame(query)
  ref <- as.data.frame(ref)
  
  #Populate data frame
  query$UpstreamGene <- ref$SYMBOL[upstreamGenes]
  query$UpstreamDistance <- NA
  query$UpstreamDistance[which(!is.na(query$UpstreamGene))] <- upstreamDist*-1
  
  query$DownstreamGene <- ref$SYMBOL[downstreamGenes]
  query$DownstreamDistance <- NA
  query$DownstreamDistance[which(!is.na(query$DownstreamGene))] <- downstreamDist
  
  query$NearestGene <- ref$SYMBOL[nearestGenes]
  query$NearestDistance <- NA
  query$NearestDistance[which(!is.na(query$NearestGene))] <- nearestDist
  
  return(query)
  
}

makeDFListfromDBAforMatrix <- function(DBAobj){
  peak_list <- DBAobj$peaks
  edited_list <- lapply(peak_list, function(df){
    # Concatenate row-wise
    df$id <- paste(df$seqnames, df$start, df$end, sep=":")
    return(df)
  })
}

sumAllMatrixRows <- function(mat, n){
  sum_rows <- apply(mat, MARGIN = 1, sum) #could just use rowSums func
  ordered <- sum_rows[order(-sum_rows)]
  filtered_mat <- mat[names(ordered[1:n]),]
  return(filtered_mat)
}

varAllMatrixRows <- function(mat, n){
  var_rows <- apply(mat, MARGIN = 1, var) #could use rowVar #names(ordered(1:n))
  ordered <- var_rows[order(-var_rows)]
  filtered_mat <- mat[names(ordered[1:n]),]
  return(filtered_mat)
}

varIKZF1vsWTRownames <- function(mat, n){
  # Perform column addition
  col1 <- mat[, "55.3 FC"] + mat[, "55.7 FC"]
  col2 <- mat[, "Par FC"] + mat[, "57.2 FC"]
  
  # Create a new matrix (or data frame) from the added columns
  new_matrix <- cbind(col1, col2) 
  
  var_rows <- apply(new_matrix, MARGIN = 1, var)
  ordered <- var_rows[order(-var_rows)]
  filtered_mat <- new_matrix[names(ordered[1:n]),]
  return(rownames(filtered_mat))
}

sumIKZF1vsWTRownames <- function(mat, n){
  # Perform column addition
  col1 <- mat[, "55.3 FC"] + mat[, "55.7 FC"]
  col2 <- mat[, "Par FC"] + mat[, "57.2 FC"]
  
  # Create a new matrix (or data frame) from the added columns
  new_matrix <- cbind(col1, col2) 
  
  sum_rows <- apply(new_matrix, MARGIN = 1, sum)
  ordered <- sum_rows[order(-sum_rows)]
  filtered_mat <- new_matrix[names(ordered[1:n]),]
  return(rownames(filtered_mat))
}

# ---------------------------------------------------------------------------
#     ****    5. CREATING DBA OBJECTS WITH DIFFERENT PEAKSETS       **** 
# ---------------------------------------------------------------------------

data_macs2 <- dba(sampleSheet = dexCsv_macs2) %>%
  dba.count(bParallel = TRUE, summits = 200) 
  dba.normalize(data_macs2)

data_seacr1 <- dba(sampleSheet = dexCsv_seacr_v1) %>%
  dba.count(bParallel = TRUE, summits = 200) %>%
  dba.normalize() %>%
  dba.contrast(contrast=c("Treatment", "DMSO", "50nM_Dex")) %>% 
  dba.analyze(method = DBA_DESEQ2) 

data_seacr2 <- dba(sampleSheet = dexCsv_seacr_v2) %>%
  dba.count(bParallel = TRUE, summits = 200) %>%
  dba.normalize() %>%
  dba.contrast(contrast=c("Treatment", "DMSO", "50nM_Dex")) %>% 
  dba.analyze(method = DBA_DESEQ2)


# ---------------------------------------------------------------------------
#       ****      6.  MAKING MATRICES FROM DBA OBJECTS        ****
# ---------------------------------------------------------------------------

# 1. Extract the list of data frames from peaks column inside of DBA object
seacr1_peaklist <- makeDFListfromDBAforMatrix(data_seacr1)
seacr2_peaklist <- makeDFListfromDBAforMatrix(data_seacr2)

# makeDFListfromDBA function doesn't work on macs2 peakcaller =>

macs2_peaklist <- lapply(data_macs2$peaks, function(df){
  df$id <- paste(df$Chr, df$Start, df$End, sep=":")
  df$sample <- 
  return(df)
})

# 2. Set the rownames for the data frames as the "id" that we generated in the last step
# * This step also merges the list of data frames into one single matrix
# ** This step requires all row names be shared across all data frames!!

seacr1_matrix <- makeMatrixFromDFList(seacr1_peaklist, seacr1_names)
seacr2_matrix <- makeMatrixFromDFList(seacr2_peaklist, seacr2_names)
macs2_matrix <- makeMatrixFromDFList(macs2_peaklist, macs2_names)


# ---------------------------------------------------------------------------
#            ****   7.  COMBINING PEAK CALLERS INTO ONE DF        ****
# ---------------------------------------------------------------------------

# Filter matrices for each peak caller to get 75-100 rows of usable data

filtered_macs2_sum <- sumAllMatrixRows(macs2_matrix, 150)
filtered_seacr1_sum <- sumAllMatrixRows(seacr1_matrix, 150)
filtered_seacr2_sum <- sumAllMatrixRows(seacr2_matrix, 150)

filtered_macs2_var <- varAllMatrixRows(macs2_matrix, 150)
filtered_seacr1_var <- varAllMatrixRows(seacr1_matrix, 150)
filtered_seacr2_var <- varAllMatrixRows(seacr2_matrix, 150)

# Average and Combine Replicates

filtered_seacr1_sum <- avg_replicates(as.data.frame(filtered_seacr1_sum))
filtered_seacr1_var <- avg_replicates(as.data.frame(filtered_seacr1_var))

filtered_seacr2_sum <- avg_replicates(as.data.frame(filtered_seacr2_sum))
filtered_seacr2_var <- avg_replicates(as.data.frame(filtered_seacr2_var))

filtered_macs2_sum <- avg_replicates(as.data.frame(filtered_macs2_sum))
filtered_macs2_var <- avg_replicates(as.data.frame(filtered_macs2_var))

seacr1_mat_union <- union_matrices(as.matrix(filtered_seacr1_sum), as.matrix(filtered_seacr1_var))
seacr2_mat_union <- union_matrices(as.matrix(filtered_seacr2_sum), as.matrix(filtered_seacr2_var))
macs2_mat_union <- union_matrices(as.matrix(filtered_macs2_sum), as.matrix(filtered_macs2_var))

# Calculate Fold Changes

log_seacr1 <- log2(seacr1_mat_union+pseudocount)
log_seacr2 <- log2(seacr2_mat_union+pseudocount)
log_macs2 <- log2(macs2_mat_union+pseudocount)

seacr1_FC <- process_matrix(log_seacr1)
seacr2_FC <- process_matrix(log_seacr2)
macs2_FC <- process_matrix(log_macs2)

# Combine Fold Change Matrices into one Matrix

union_FC <- union_matrices(macs2_FC, seacr1_FC, seacr2_FC)


# Turn the Fold Change Matrix into a Data Frame for Gene Annotating Step

union_FC_df <- as.data.frame(union_FC, row.names=rownames(union_FC))
union_FC_df <- union_FC_df %>%
  rownames_to_column('row_names') %>% 
  # Next step is to get rid of any rows that don't fit the "seqname_start_end" 
  #   format.
  filter(str_detect(row_names, ":") & stringr::str_count(row_names, ":") == 2) %>%
  separate_wider_delim(
    cols = row_names,
    delim = ":",
    names = c('seqnames', 'start', 'end'))


# ----------------------------------------------------------------------------
# ****  8.  CREATING BED FILES FROM PEAK CALLERS TO VIEW ON IGV * OPTIONAL   ****
# ----------------------------------------------------------------------------

#Output as GRANGES Objects

union_gr <- makeGRangesFromDataFrame(
                union_FC_df,
                keep.extra.columns = TRUE,
                seqnames.field = 'seqnames',
                start.field = 'start',
                end.field = 'end'
)


# Export as bed files
union_bed <- 'bed/union_peaks_granges.bed'
export(union_gr, union_bed)

# ---------------------------------------------------------------------------
#                ****       9. ANNOTATING GENES          *****
# ---------------------------------------------------------------------------

GR <- PrepareAnnotationGrangesObject(txdb)

#Annotate regions
peaks <- makeGRangesFromDataFrame(union_FC_df, keep.extra.columns = TRUE) %>%
  ChIPseeker::annotatePeak(TxDb = txdb, 
                           level = "gene", 
                           annoDb="org.Hs.eg.db", 
                           verbose = F, 
                           addFlankGeneInfo = F) %>%
  as.data.frame() %>%
  #subset(grepl("Promoter", annotation) | grepl("Distal Intergenic", annotation)) %>%
  #subset(grepl("Promoter", annotation)) %>%
  #subset(!is.na(ENSEMBL)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  sortSeqlevels() %>%
  sort(ignore.strand = TRUE) %>%
  as.data.frame() %>%
  mutate(ID = paste(seqnames, start,end, sep = ":")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  GetFlankingGenes(GR)




# ---------------------------------------------------------------------------
#     ****       10. ANNOTATING THE HEATMAP WITH GENE_LIST          *****
# ---------------------------------------------------------------------------

peaks_mat <- as.matrix(peaks[, c("X55.3.FC", 'X55.7.FC', 'X57.2.FC', 'Par.FC')])
rownames(peaks_mat)<-peaks$ID
colnames(peaks_mat)<- c('55.3 FC', '55.7 FC', '57.2 FC', "Par FC")

peaks_htmap <- Heatmap(
  peaks_mat,
  show_row_dend = FALSE,
  use_raster = FALSE,
  show_row_names = FALSE,
  cluster_rows = TRUE,
  column_order = c("55.3 FC", '55.7 FC', 'Par FC', '57.2 FC'),
  show_column_dend = FALSE,
  name = 'Score',
  column_title = "Fold Change Heatmap"
)
draw(peaks_htmap)

row_order <- row_order(peaks_htmap)
clustered_matrix <- peaks_mat[row_order,]

var_50_rows <- varIKZF1vsWTRownames(clustered_matrix, 50)
sum_50_rows <- sumIKZF1vsWTRownames(clustered_matrix, 50)


mask <- peaks$ID %in% sum_50_rows
gene_peaks <- filter(peaks, mask)
htmap_mark <- peaks$ID
mark_index <- match(htmap_mark, rownames(clustered_matrix))
htmap_label <- peaks$SYMBOL

heatmap_anno <- rowAnnotation(foo = anno_mark(at = mark_index, labels = htmap_label))

clust_htmap <- Heatmap(
  clustered_matrix,
  show_row_dend = FALSE,
  use_raster = FALSE,
  show_row_names = FALSE,
  cluster_rows = F,
  column_order = c("55.3 FC", '55.7 FC', 'Par FC', '57.2 FC'),
  show_column_dend = FALSE,
  right_annotation = heatmap_anno,
  name = 'Score',
  column_title = "Clustered Heatmap"
)
ht_shiny(clust_htmap)

# Reorder columns instead of clustering columns to separate IKZF1 KO and WT
# Find literature that 



# ---------------------------------------------------------------------------
#            ****       SAVING/LOADING SESSION          ****
# ---------------------------------------------------------------------------


union_FC_df_v2 <- readRDS('rds/union_FC_df_v2.rds')
gene_names <- readRDS('rds/gene_names.rds')
union_FC_df <- readRDS('rds/union_FC_df.rds')
union_FC <- readRDS('rds/union_FC.rds')
macs2_FC <- readRDS('rds/macs2_FC.rds')
seacr2_FC <- readRDS('rds/seacr2_FC.rds')
seacr1_FC <- readRDS('rds/seacr1_FC.rds')

macs2_matrix <- readRDS('rds/macs2_matrix.rds')
seacr_v1_matrix <- readRDS('rds/seacr1_matrix.rds')
seacr_v2_matrix <- readRDS('rds/seacr2_matrix.rds')

macs2_df <- readRDS('rds/macs2_df.rds')
seacr1_df <- readRDS('rds/seacr1_df.rds')
seacr2_df <- readRDS('rds/seacr2_df.rds')

data_seacr1 <- readRDS('rds/data_seacr1.rds')
data_seacr2 <- readRDS('rds/data_seacr2.rds')
data_macs2 <- readRDS('rds/data_macs2.rds')

saveRDS(union_FC_df, 'rds/union_FC_df.rds')
saveRDS(union_FC, 'rds/union_FC.rds')
saveRDS(macs2_FC, 'rds/macs2_FC.rds')
saveRDS(seacr2_FC, 'rds/seacr2_FC.rds')
saveRDS(seacr1_FC, 'rds/seacr1_FC.rds')

saveRDS(merged_macs2_df, "merged_macs2_df.rds")
saveRDS(merged_sv1_df, 'merged_sv1_df.rds')
saveRDS(merged_sv2_df, 'merged_sv2_df.rds')

saveRDS(macs2_df, "macs2_df.rds")
saveRDS(seacr1_df, "seacr1_df.rds")
saveRDS(seacr2_df, "seacr2_df.rds")

saveRDS(data_seacr1, "data_seacr1.rds")
saveRDS(data_seacr2, "data_seacr2.rds")
saveRDS(data_macs2, "data_macs2.rds")
saveRDS(macs2_matrix, "macs2_matrix.rds")
saveRDS(seacr_v1_matrix, "seacr1_matrix.rds")
saveRDS(seacr_v2_matrix, "seacr2_matrix.rds")
