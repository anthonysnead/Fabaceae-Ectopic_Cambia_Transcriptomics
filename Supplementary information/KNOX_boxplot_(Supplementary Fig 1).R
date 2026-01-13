# Script written by Israel L. Cunha-Neto
# AI-assisted development using ChatGPT (OpenAI)

install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("limma", force = TRUE)
BiocManager::install("BiocParallel")
library(edgeR)
library(tibble)
library(tidyr)
library(dplyr)
library(tidyverse)

setwd("")

#Extracting table with expression data for all genes from RData
load("Ectopic.RData") #file containing all clusters from the DGE analysis
ls()
class(Ortho.DGE.Data)
dim(Ortho.DGE.Data)
head(Ortho.DGE.Data)
counts_matrix <- Ortho.DGE.Data$counts
counts_df <- counts_matrix %>% #Create a table with clusters for each sample
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Cluster")

write.csv(counts_df, file = "Ortho_DGE_Counts.csv", row.names = FALSE)

#Read table
expr_matrix <- read.csv("Ortho_DGE_Counts.csv", row.names = 1)

#Filter the clusters to KNOX genes only
clusters_to_keep <- c("Cluster_9937", "Cluster_2082", "Cluster_4072", "Cluster_1332")
subset_matrix <- expr_matrix[rownames(expr_matrix) %in% clusters_to_keep, ]

#Reformat for plotting with KNOX genes only 
df_long <- subset_matrix %>%
  rownames_to_column(var = "Cluster") %>%
  pivot_longer(-Cluster, names_to = "Sample", values_to = "Expression") 

#Adding section number
section_info <- read.csv("section_info.csv")
df_long.3 <- left_join(df_long.2, section_info, by = "Sample")
write.csv(df_long.2, file = "KNOX.Expression.Sections.csv")

#Plotting boxplots per cambium types
ggplot(df_long.3, aes(x = Section, y = Expression, fill = SampleType)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = SampleType),
              width = 0.1, height = 0,
              size = 1.5, alpha = 0.7) +
  facet_wrap(~Cluster, scales = "free_y") +
  theme_minimal() +
  labs(x = "Composite Section Number",
       y = "Normalized Expression",
       title = "Expression Pattern Across Cryosectioned Material")

#End
