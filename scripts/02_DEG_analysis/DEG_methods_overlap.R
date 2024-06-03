library(dplyr)
library(readr)

############
## DeSeq2 ##
############

# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/DESeq2/M_MUT_and_WT_M_P30_CORT")
# List all the CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")
# Read each CSV file and combine them into a single dataframe
DeSeq2_P30 <- do.call(rbind, lapply(csv_files, read.csv))
DeSeq2_P30 <- DeSeq2_P30 %>%
  distinct(X, .keep_all = TRUE)
DeSeq2_P30$SYMBOL <- DeSeq2_P30$X

# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/DESeq2/M_MUT_and_WT_M_P60_CORT")
# List all the CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")
# Read each CSV file and combine them into a single dataframe
DeSeq2_P60 <- do.call(rbind, lapply(csv_files, read.csv))
DeSeq2_P60 <- DeSeq2_P60 %>%
  distinct(X, .keep_all = TRUE)
DeSeq2_P60$SYMBOL <- DeSeq2_P60$X

# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/DESeq2/M_MUT_and_WT_M_P120_CORT")
# List all the CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")
# Read each CSV file and combine them into a single dataframe
DeSeq2_P120 <- do.call(rbind, lapply(csv_files, read.csv))
DeSeq2_P120 <- DeSeq2_P120 %>%
  distinct(X, .keep_all = TRUE)
DeSeq2_P120$SYMBOL <- DeSeq2_P120$X

############
## EdgeR ##
############
# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/EdgeR/M_MUT_and_WT_M_P30_CORT")
# List all the CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")
# Read each CSV file and combine them into a single dataframe
EdgeR_P30 <- do.call(rbind, lapply(csv_files, read.csv))
EdgeR_P30 <- EdgeR_P30 %>%
  distinct(X, .keep_all = TRUE)
EdgeR_P30$SYMBOL <- EdgeR_P30$X

# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/EdgeR/M_MUT_and_WT_M_P60_CORT")
# List all the CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")
# Read each CSV file and combine them into a single dataframe
EdgeR_P60 <- do.call(rbind, lapply(csv_files, read.csv))
EdgeR_P60 <- EdgeR_P60 %>%
  distinct(X, .keep_all = TRUE)
EdgeR_P60$SYMBOL <- EdgeR_P60$X

# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/EdgeR/M_MUT_and_WT_M_P120_CORT")
# List all the CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")
# Read each CSV file and combine them into a single dataframe
EdgeR_P120 <- do.call(rbind, lapply(csv_files, read.csv))
EdgeR_P120 <- EdgeR_P120 %>%
  distinct(X, .keep_all = TRUE)
EdgeR_P120$SYMBOL <- EdgeR_P120$X

############
## Limma ##
############
# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_P30_CORT")
csv_files <- list.files(pattern = "\\.csv$")
df_list <- list()
for (file in csv_files) {
  tryCatch({
    df <- read.csv(file, encoding = "UTF-8", stringsAsFactors = FALSE)
    df_list[[file]] <- df
  }, error = function(e) {
    cat("Error reading file:", file, "\n")
    print(e)
  })
}
Limma_P30 <- do.call(rbind, df_list)
Limma_P30 <- Limma_P30 %>%
  distinct(X, .keep_all = TRUE)
Limma_P30$SYMBOL <- Limma_P30$X
rm(df, df_list)

#P60
# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_P60_CORT")
csv_files <- list.files(pattern = "\\.csv$")
df_list <- list()
for (file in csv_files) {
  lines <- readLines(file, warn = FALSE)
  tryCatch({
    df <- read.csv(text = lines, header = TRUE)
    df_list[[file]] <- df
  }, error = function(e) {
    cat("Error reading file:", file, "\n")
    print(e)
  })
}
Limma_P60 <- do.call(rbind, df_list)
Limma_P60 <- Limma_P60 %>%
  distinct(X, .keep_all = TRUE)
Limma_P60$SYMBOL <- Limma_P60$X
rm(df, df_list)

#P120
# Set the working directory to the path where your CSV files are located
setwd("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/Limma/M_MUT_and_WT_M_P120_CORT")
csv_files <- list.files(pattern = "\\.csv$")
df_list <- list()
for (file in csv_files) {
  lines <- readLines(file, warn = FALSE)
  tryCatch({
    df <- read.csv(text = lines, header = TRUE)
    df_list[[file]] <- df
  }, error = function(e) {
    cat("Error reading file:", file, "\n")
    print(e)
  })
}
Limma_P120 <- do.call(rbind, df_list)
Limma_P120 <- Limma_P120 %>%
  distinct(X, .keep_all = TRUE)
Limma_P120$SYMBOL <- Limma_P120$X
rm(df, df_list)

#################
## LimmaVoomCC ##
#################
# Set the working directory to the path where your CSV files are located
Limmavoomcc <- read.csv("/Users/osman/Documents/GitHub/original-snRNA-seq-pipeline/DEG_data/master_df.csv")
filtered_df <- Limmavoomcc %>%
  filter(method == "limmaVoomCC" & sex == "M" & adj.P.Val <= 0.05)
limmvoomcc_P30 <- filtered_df %>%
  filter(timepoint == "P30")
limmvoomcc_P60 <- filtered_df %>%
  filter(timepoint == "P60")
limmvoomcc_P120 <- filtered_df %>%
  filter(timepoint == "P120")
rm(filtered_df, Limmavoomcc)


library(ggVennDiagram)
P30 <- list(DEseq2 = DeSeq2_P30$SYMBOL,
            EdgeR = EdgeR_P30$SYMBOL,
            Limma = Limma_P30$SYMBOL,
            LimmaVoomCC = limmvoomcc_P30$SYMBOL)
ggVennDiagram(P30) + ggtitle("Overlap of significant DEGs P30")
ggplot2::ggsave("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/method_comparison/DEG_method_overlap_P30.pdf",
                device = NULL,
                height = 10,
                width = 12)
P60 <- list(DEseq2 = DeSeq2_P60$SYMBOL,
            EdgeR = EdgeR_P60$SYMBOL,
            Limma = Limma_P60$SYMBOL,
            LimmaVoomCC = limmvoomcc_P60$SYMBOL)
ggVennDiagram(P60) + ggtitle("Overlap of significant DEGs P60")
ggplot2::ggsave("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/method_comparison/DEG_method_overlap_P60.pdf",
                device = NULL,
                height = 10,
                width = 12)
P120 <- list(DEseq2 = DeSeq2_P120$SYMBOL,
            EdgeR = EdgeR_P120$SYMBOL,
            Limma = Limma_P120$SYMBOL,
            LimmaVoomCC = limmvoomcc_P120$SYMBOL)
ggVennDiagram(P120) + ggtitle("Overlap of significant DEGs P120")
ggplot2::ggsave("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/method_comparison/DEG_method_overlap_P120.pdf",
                device = NULL,
                height = 10,
                width = 12)

test <- read.csv("/Users/osman/Downloads/total_sig_mouse_DEGs_limmaVoom.csv")
filtered_df <- test %>%
  filter(Sex == "males")

modules <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/08_hdWGCNA_analysis/modules.csv")
DEG_wgcna_overlap <- list(LimmaVoomCC_DEGs = test$SYMBOL,
                          hdWGCNA_Genes = modules$gene_name)
ggVennDiagram(DEG_wgcna_overlap) + ggtitle("Overlap of DEGs and hdWGCNA genes")
ggplot2::ggsave("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/method_comparison/DEG_hdWGCNA_overlap.pdf",
                device = NULL,
                height = 10,
                width = 12)
