# ========================================================
# Transcriptomic Analysis of Basal Breast Cancer
# GEO Accession: GSE45827
# 
# Objective:
# Identification of genes differentially expressed between
# basal breast tumors and normal breast tissue.
# ========================================================


# =========================================================
# Required Packages:
# =========================================================

library(tidyverse)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)


# =========================================================
# Project Directories:
# =========================================================

dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

# =========================================================
# GEO Data Acquisition:
# =========================================================

gset <- getGEO(
  "GSE45827",
  GSEMatrix = TRUE
)

# =========================================================
# Expression Matrix:
# =========================================================

expr_matrix <- exprs(gset[[1]])

# ========================================================
# Phenotype Metadata:
# ========================================================

pheno <- pData(gset[[1]])

# =========================================================
# Dataset Overview:
# =========================================================

dim(expr_matrix)

head(expr_matrix[, 1:5])

table(pheno$`tumor subtype:ch1`)

# =========================================================
# Tumor Subtype Annotation:
# =========================================================

tumor_subtype <- pheno$`tumor subtype:ch1`

tumor_subtype[tumor_subtype == "N/A"] <- "Normal"
tumor_subtype[tumor_subtype == "Luminal A"] <- "LuminalA"
tumor_subtype[tumor_subtype == "Luminal B"] <- "LuminalB"

tumor_subtype <- factor(tumor_subtype)

# ========================================================
# Subtype distribution:
# ========================================================

table(tumor_subtype)

levels(tumor_subtype)

# =========================================================
# Sample Filtering:
# =========================================================

valid_sample <- !is.na(tumor_subtype)

expr_matrix_filtered <- expr_matrix[, valid_sample]

tumor_subtype_filtered <- factor(tumor_subtype[valid_sample])

# ========================================================
# Phenotype Encoding:
# ========================================================

levels(tumor_subtype_filtered) <- c(
  "Basal",
  "Her2",
  "LuminalA",
  "LuminalB",
  "Normal"
)

# ========================================================
# Experimental Design Matrix:
# ========================================================

design <- model.matrix(~0 + tumor_subtype_filtered)

colnames(design) <- levels(tumor_subtype_filtered)

design

# =======================================================
# Linear Model Fitting:
# =======================================================

fit_model <- lmFit(
  expr_matrix_filtered,
  design
)

# ======================================================
# Contrast Specification:
# Basal Breast Cancer vs Normal Tissue:
# ======================================================

contrast_matrix <- makeContrasts(
  BasalvsNormal = Basal - Normal,
  levels = design
)

contrast_matrix

# ======================================================
# Empirical Bayes Moderation:
# ======================================================

fit_contrast <- contrasts.fit(
  fit_model,
  contrast_matrix
)

fit_contrast <- eBayes(fit_contrast)

# ======================================================
# Differential Expression Results:
# ======================================================

deg_results <- topTable(
  fit_contrast,
  adjust = "fdr",
  number = Inf
)

head(deg_results)

# =======================================================
# Significant Differentially Expressed Genes:

# Thresholds:
# Adjusted P-value < 0.05
# Absolute Log2 Fold Change > 1
# =======================================================

significant_genes <- deg_results %>%
  filter(
    adj.P.Val < 0.05 &
      abs(logFC) > 1
  )

  dim(significant_genes)
  
  head(significant_genes)

# =====================================================
# Export Differential Expression Results:
# =====================================================

write.csv(
  deg_results,
  "outputs/full_differential_expression_results.csv"
)
  
write.csv(
  significant_genes,
  "outputs/significant_genes.csv"
)

# ======================================================
# Volcano Plot Visualization:
# ======================================================

deg_results$Significance <- as.factor(
  deg_results$adj.P.Val < 0.05 & 
    abs(deg_results$logFC) > 1
)

volcano_plot <- ggplot(
  deg_results,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Significance
  )
) + 
  geom_point(
    alpha = 0.7,
    size = 1.5
  ) + 
  theme_minimal() +
  labs(
    title = "Volcano Plot: Basal vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

  volcano_plot

# ========================================================
# Volcano Plot Export:
# ========================================================
  
ggsave(
  "figures/volcano_plot.png",
  volcano_plot,
  width = 8,
  height = 6
)

ggsave(
  "figures/volcano_plot.pdf",
  volcano_plot,
  width = 8,
  height = 6
)

# ========================================================
# Top Differentially Expressed Genes
# ========================================================

top_genes <- rownames(
  significant_genes
) [1:50]

heatmap_matrix <- expr_matrix_filtered[
  top_genes,
]

annotation_data <- data.frame(
  Subtype = tumor_subtype_filtered
  )

rownames(annotation_data) <- colnames(
  heatmap_matrix
)

# ========================================================
# Heatmap Visualization:
# ========================================================

heatmap_plot <- pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_data,
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top Differentially Expressed Genes"
)

# ========================================================
# Heatmap Export:
# ========================================================

png(
  "figures/heatmap_top_genes.png",
  width = 1200,
  height = 1000
)

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_data,
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top Differentially Expressed Genes"
)

dev.off()

pdf(
  "figures/heatmap_top_genes.pdf",
  width = 10,
  height = 8
)

pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_data,
  show_rownames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Top Differentially Expressed Genes"
)

dev.off()

# =========================================================
# Session Information:
# =========================================================

sessionInfo()

# Note:
# Analysis performed on microarray-based transcriptomic data.
# Experimental validation was not inlcuded in this workflow.
