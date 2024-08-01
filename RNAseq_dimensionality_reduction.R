# VARIABLES
case <- ""
sample <- ""
salmon_path <- ""
date <- Sys.Date()
# Unfortunately you will need to download the Treehouse tpms
#treehouse_tpm_file_location <- paste0(getwd(), "/transcript_tpm_data/Treehouse_tpms/tpm_pedaya_treehouse.tsv")
treehouse_tpm_file_location <- paste0(getwd(),"/transcript_tpm_data/Treehouse_tpms/TPMFILENAMEHERE.tsv")


# IMPORT LIBRARIES
library(tidyverse)
library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)
library(vroom)
library(umap)
library(rsvg)
library(ggrepel)
library(RColorBrewer)
library(tibble)
library(UpSetR)
library(matrixStats)
library(scales)
library(tximport)


# CREATE RESULTS DIRECTORIES
results.dir <- paste0(getwd(), "/results/")
case.dir <- paste0(results.dir,case,"/")
working.dir <- paste0(case.dir,sample,"/")
matrix.dir <- paste0(working.dir,"TPM_analysis/")
outlier.dir <- paste0(working.dir,"Protein_Coding_Outliers/")
dimred.dir <- paste0(working.dir,"Dimensionality_Reduction/")
top.dir <- getwd()
meta.dir <- paste0(top.dir,"/metadata/")
top.tpm.dir <- paste0(top.dir,"/transcript_tpm_data/")
tpm.dir <- paste0(top.dir,"/transcript_tpm_data/Internal_tpms/")
tpm.th.dir <- paste0(top.dir,"/transcript_tpm_data/Treehouse_tpms/")
sample.tpm.dir <- paste0(tpm.dir, sample,".cov")
# Make this directory if it does not already exist
ifelse(!dir.exists(file.path(results.dir)), dir.create(file.path(results.dir)), FALSE)
ifelse(!dir.exists(file.path(case.dir)), dir.create(file.path(case.dir)), FALSE)
ifelse(!dir.exists(file.path(working.dir)), dir.create(file.path(working.dir)), FALSE)
ifelse(!dir.exists(file.path(matrix.dir)), dir.create(file.path(matrix.dir)), FALSE)
ifelse(!dir.exists(file.path(outlier.dir)), dir.create(file.path(outlier.dir)), FALSE)
ifelse(!dir.exists(file.path(dimred.dir)), dir.create(file.path(dimred.dir)), FALSE)
ifelse(!dir.exists(file.path(top.tpm.dir)), dir.create(file.path(top.tpm.dir)), FALSE)
ifelse(!dir.exists(file.path(tpm.dir)), dir.create(file.path(tpm.dir)), FALSE)
ifelse(!dir.exists(file.path(tpm.th.dir)), dir.create(file.path(tpm.th.dir)), FALSE)
ifelse(!dir.exists(file.path(top.tpm.dir)), dir.create(file.path(top.tpm.dir)), FALSE)


# IMPORT INTERNAL COHORT METADATA
meta <- read_xlsx(paste0(meta.dir, "metadata_internal_cohort.xlsx"))
meta$source <- "Internal"
dir.exists(file.path(sample.tpm.dir))
# Copy file over if it does not exist
ifelse(!dir.exists(file.path(sample.tpm.dir)), R.utils::copyDirectory(salmon_path, sample.tpm.dir, recursive=TRUE), FALSE)


# IMPORT INTERNAL COHORT TPMS
salmon_paths <-file.path(tpm.dir, meta$salmon_gencode, "quant.sf")
names(salmon_paths) <- meta$salmon_gencode
salmon_paths <- as.matrix(salmon_paths)
salmon_paths_df <- data.frame(salmon_paths)
salmon_paths_df$present <- ifelse((file.exists(file.path(salmon_paths))), print("TRUE"), print(paste0("TPM FILE MISSING")))
salmon_paths_with_tpm <- subset(salmon_paths_df,present==TRUE)
salmon_paths_with_tpm$salmon_gencode <- row.names(salmon_paths_with_tpm)
salmon_paths_with_tpm_df <- salmon_paths_with_tpm[, c("salmon_gencode", "salmon_paths")]
salmon_paths_with_tpm_final <- with(salmon_paths_with_tpm_df, setNames(salmon_paths, salmon_gencode))
salmon_paths <- salmon_paths_with_tpm_final
tx_annotation <- vroom::vroom(paste0(getwd(),"/accessory_files/salmon_tx_annotation.tsv"))
tx2gene <- tx_annotation %>%
  dplyr::select(Name, hgnc_symbol) %>%
  distinct()
txi_salmon_paths <- as.vector(salmon_paths)
names(txi_salmon_paths) <- names(salmon_paths)
txi <- tximport::tximport(
  txi_salmon_paths,
  type = "salmon",
  tx2gene = tx2gene,
  dropInfReps = TRUE,
  importer = function(x)
    vroom::vroom(
      x,
      col_types = cols(
        Name = col_character(),
        Length = col_double(),
        EffectiveLength = col_double(),
        TPM = col_double(),
        NumReads = col_double()
      )
    )
)
write_rds(txi, paste0(matrix.dir, sample, "_tximport_data.rds"))
tpm <- txi$abundance


# TRANSFORM INTERNAL COHORT TPMS
tpm_nonlog <- tpm
tpm %>% as_tibble(rownames = "Gene")
tpm <- log2(tpm + 1) %>%  as_tibble(rownames = "Gene")
names(tpm)=chartr(".", "-", as.character(names(tpm)))
names(tpm) = gsub(pattern = ".salmon.cov", replacement = "", x = names(tpm))
names(tpm) = gsub(pattern = ".cov", replacement = "", x = names(tpm))


# VERIFY INTERNAL COHORT MATCHES
meta <- dplyr::inner_join(meta,salmon_paths_with_tpm_df)
meta$metadata_id = chartr(".", "-", as.character(meta$metadata_id))
pedaya_samp_ids <- meta$metadata_id
ped <- names(tpm)
print(fromList(list(
  tpm = colnames(tpm)[-1],
  meta = meta$metadata_id
)) %>%
  upset(., text.scale = 1.6))
tpm <- tpm %>% dplyr::select(all_of(c("Gene", meta$metadata_id)))


# IMPORT TREEHOUSE METADATA
th_meta <- vroom::vroom(paste0(meta.dir,"treehouse_metadata.txt"))
# set Diagnosis (since thats what we use for IGM) as the column header for disease
cancer_map_path <- dir(meta.dir, "treehouse_v11_v9_pedaya_CNS_diagnoses.xlsx",
                         full.names = TRUE)
names(th_meta)[names(th_meta) == "disease"] <- "Diagnosis"
cancer_map <- readxl::read_xlsx(cancer_map_path)
names(cancer_map)[names(cancer_map) == "disease"] <- "Diagnosis"
# join subtyping to external metadata
th_meta <- th_meta %>% left_join(cancer_map)
#rename the th_sampleid to metadata_id so that it can fit better with the IGM metadata
names(th_meta)[names(th_meta) == "th_sampleid"] <- "metadata_id"
th_meta$source <- "Treehouse"


# IMPORT TREEHOUSE TPMS
# You will need to download the following data from the source. It is too large for github.
treehouse_tpm <- vroom::vroom(treehouse_tpm_file_location)
th_tpm <- treehouse_tpm
CNS_th_ids <- th_meta$metadata_id
th_log2_tpm <- th_tpm %>% dplyr::select(all_of(c("Gene", CNS_th_ids)))


# JOIN INTERNAL AND EXTERNAL COHORT TPMS AND METADATA
# Create final TPM with all data
joined_tpm <- tpm %>%
  inner_join(tpm) %>%
  inner_join(th_log2_tpm)
# Join metadata
th_meta <- as.data.frame(th_meta)
meta <- as.data.frame(meta)
joined_meta <- rbind.fill(th_meta, meta)
# Check sources
headername <- names(joined_tpm)
pedaya_samp_ids <- joined_meta$metadata_id
pedaya_log2_tpm <- joined_tpm %>%
  dplyr::select(all_of(c("Gene", joined_meta$metadata_id)))
sources <- dplyr::select(joined_meta, source, cancer_grp)
sources_total <- sources %>% 
  dplyr::group_by(source) %>% 
  dplyr::summarize(n=n()) %>% 
  dplyr::mutate(Frequency=n)
sources_total
# Select only protein coding genes
protein_coding <- read_csv(paste0(getwd(), "/accessory_files/Gene_Type_Salmon.csv"))
pedaya_log2_tpm_pc_join <- pedaya_log2_tpm %>%
  dplyr::inner_join(protein_coding)
pedaya_log2_tpm_pc <- pedaya_log2_tpm_pc_join[which(pedaya_log2_tpm_pc_join$Type=="protein_coding"),]
# Remove type column, save file
pedaya_log2_tpm_pc <- dplyr::select(pedaya_log2_tpm_pc, -c(Type))
write_tsv(pedaya_log2_tpm_pc, paste0(matrix.dir, date, "_log2_tpm_all_pc_genes.tsv"))


# IDENTIFY OUTLIERS
pedaya_log2_tpm_mat <- as.matrix(pedaya_log2_tpm_pc[,-1])
row.names(pedaya_log2_tpm_mat) <- pedaya_log2_tpm_pc[[1]]
outlier_detection <- function(tpm_mat, sample_id){
  # Quantile normalize data
  tpm_dimnames <- dimnames(pedaya_log2_tpm_mat)
  norm_tpm <- pedaya_log2_tpm_mat %>%
    preprocessCore::normalize.quantiles()
  dimnames(norm_tpm) <- tpm_dimnames
  # Compute fold change for sample and cbind to TPM
  gene_medians <- matrixStats::rowMedians(norm_tpm[, !colnames(norm_tpm) %in% sample_id])
  fc_sample <- norm_tpm[, sample_id] - gene_medians
  tpm_fc_mat <- cbind(sample_id=norm_tpm[, sample_id], fc_sample)
  # Compute mahalanobis distance
  mahal_dist <- mahalanobis(tpm_fc_mat, colMeans(tpm_fc_mat), cov(tpm_fc_mat))
  # Compute pvalue
  pval <- pchisq(mahal_dist, df = 2, lower.tail=FALSE)
  # Combine all in data frame
  dist_df <- as.data.frame(tpm_fc_mat) %>%
    mutate(mahal_dist = mahal_dist,
           pvalue = pval,
           fold_change = fc_sample)
}

outliers <- outlier_detection(pedaya_log2_tpm_mat,sample)
outliers <- outliers %>% dplyr::rename(tpm = sample_id)
# Add column for outlier status
outliers$outlier_status <- "NotOutlier"
outliers$outlier_status[(outliers$fc_sample < -3.5 | outliers$fc_sample > 3.5) & outliers$pvalue < 0.005] <- "Outlier"
# Create a file for all protein coding genes and outlier status
total_outliers <- outliers
total_outliers <- tibble::rownames_to_column(total_outliers, "Gene")
write_csv(total_outliers, paste0(outlier.dir,sample,"_all_pc_genes_outlier_status.csv"))


# PRINCIPAL COMPONENTS ANALYSIS
pca <- preprocessCore::normalize.quantiles(pedaya_log2_tpm_mat)
is.numeric(pca)
ntop = 500
Pvars <- rowVars(pca)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
pcaData <- prcomp(t(pca[select, ]))
pca.proportionvariances <- ((pcaData$sdev^2) / (sum(pcaData$sdev^2)))*100
pca_df_out <- as.data.frame(pcaData$x)
rownames(pca_df_out) <- colnames(pedaya_log2_tpm_mat)
pca_tibble <- pca_df_out %>%
  as_tibble(rownames = "metadata_id")
joined_pca_data <- inner_join(pca_tibble, joined_meta, by=c("metadata_id"="metadata_id"))
write_csv(joined_pca_data, paste0(dimred.dir,sample,"_pca.csv"))


# UMAP
ID <- sample
# Import PCA data
PCs <- joined_pca_data
PCs$Diagnosis <- tolower(PCs$Diagnosis)
PCs$select <- NA
PCs$select <- ifelse(grepl(case, PCs$case), 1, 0)
# Use first 34 PCs to create matrix for further analysis
PC.data = as.matrix(PCs[,2:35])
PC.data.matrix <- data.matrix(PC.data)
colnames(PC.data.matrix) <- NULL

# Begin umap work
neighbors <- 15
run.umap <- umap(PC.data.matrix, init = "random", n_neighbors = neighbors, pca_center = FALSE, local_connectivity = 5, min_dist=0.01)
umap.dataframe <- data.frame(UMAP1 = run.umap$layout[,1],
                             UMAP2 = run.umap$layout[,2],
                             Diagnosis = PCs$Diagnosis,
                             ID = PCs$metadata_id,
                             Case = PCs$case,
                             Select = PCs$select)
umap.plot <- ggplot(umap.dataframe, aes(UMAP1, UMAP2, colour = Diagnosis)) + 
  theme_bw() + geom_hline(yintercept=0, size=0.2, color = "darkgray") + geom_vline(xintercept=0, size=0.2, color = "darkgray") +
  theme(axis.title.x = element_text(color = "black", size = 12, angle = 00, face = "bold")) +
  theme(axis.title.y = element_text(color = "black", size = 12, angle = 90, face = "bold"))  + 
  ggtitle(paste0("Cancer Sample Clustering by UMAP: CNS Tumors")) + guides(col = guide_legend(nrow = 45)) +
  geom_point(alpha = 0.6, size = 1.5) +
    theme(legend.key.size = unit(0.15, "cm"))
umap.plot <- umap.plot + guides(colour = guide_legend(ncol = 1))
subset.umap <- subset(umap.dataframe, umap.dataframe$Select == 1)
label.umap <- list(
  x = subset.umap$UMAP1,
  y = subset.umap$UMAP2,
  xref = "x",
  yref = "y",
  text = subset.umap$ID,
  showarrow = TRUE,
  arrowhead = 4,
  arrowsize = 1,
  ax = 20,
  ay = -40
)
# Create interactive plot
interactive.umap <- plotly::ggplotly(umap.plot)
interactive.umap.annotated <- interactive.umap %>% plotly::layout(annotations = label.umap)
#create interactive umap html file, but name correctly (check CH number)
htmlwidgets::saveWidget((interactive.umap.annotated), paste0(dimred.dir,"interactive_umap_",sample,".html"))
# Create static plot that is labeled by the select column
umap.static.annotated <- umap.plot + geom_text_repel(
  data = subset.umap,
  aes(label = ID),
  size = 4.5, color = "black", fontface="bold",
  box.padding = unit(3.5, "lines"),
  point.padding = unit(0.2, "lines")
)
# Save as png
ggsave(paste0(dimred.dir,sample,"_UMAP.png"), height = 10, width = 20)


# FINAL MESSAGE FOR CONSOLE
print(paste0("Results for this sample are located at ", working.dir))
