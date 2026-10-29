# =====================================================================
# Project: Single-cell RNA-seq Analysis Pipeline
# Description: Standard Seurat workflow for single-cell transcriptomic data
# =====================================================================

rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)

setwd("analysis/1st/")

# Define the sample names for analysis (modify here or loop through all samples)
# Sample list:  "c14_1", "c14_2", "c14_3", "e14_1", "e14_2", "e14_3"
sample_name <- "c14_1"  

# Read 10X cell-gene expression matrix
data_dir <- paste0("analysis/QC_Report_2026_09_21/1.MobiVision/ZB_Abx_gill_1/filtered_cell_gene_matrix")
pbmc <- Read10X(data.dir = data_dir)
pbmc <- CreateSeuratObject(counts = pbmc, project = sample_name)

# Quality control and filtering (QC)
# Note: the mitochondrial gene matching pattern may vary by species (mouse/zebrafish usually use "mt-" or "MT-")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-") 

# Filtering criteria: adjust according to the specific QC report
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 

# Normalization and selection of highly variable genes
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Data scaling and PCA dimensionality reduction
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Clustering and nonlinear dimensionality reduction (UMAP / t-SNE)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

# Find marker genes for each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Save key results
dir.create("output", showWarnings = FALSE)
write.csv(top10, file = paste0("output/", sample_name, "_top10_markers.csv"), row.names = FALSE)
saveRDS(pbmc, file = paste0("output/", sample_name, "_seurat_object.rds"))


# =====================================================================
# Project: Single-Cell RNA-Seq Analysis Pipeline (C vs E Groups)
# =====================================================================

# Initialize environment and load dependencies
rm(list = ls())
setwd("single cell/new_anlysis/ALL")

library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)

# Read multi-sample data and assign prefixes
sample_mapping <- c(
    "single cell/rowdata/c14_1/filtered_feature_bc_matrix" = "c1",
    "single cell/rowdata/C14_2/filtered_feature_bc_matrix" = "c2",
    "single cell/rowdata/C14_3/filtered_feature_bc_matrix" = "c3",
    "single cell/rowdata/E14_1/filtered_feature_bc_matrix" = "e1",
    "single cell/rowdata/E14_2/filtered_feature_bc_matrix" = "e2",
    "single cell/rowdata/E14_3/filtered_feature_bc_matrix" = "e3"
)

sample_ids <- c("c1", "c2", "c3", "e1", "e2", "e3")

seurat_list <- lapply(names(sample_mapping), function(dir) {
    data <- Read10X(data.dir = dir)
    CreateSeuratObject(counts = data, project = sample_mapping[dir], min.cells = 3, min.features = 200)
})

# Add barcode prefixes
for (i in seq_along(seurat_list)) {
    colnames(seurat_list[[i]]) <- paste0(sample_ids[i], "_", colnames(seurat_list[[i]]))
}

#  Merge data and perform basic preprocessing
combined <- merge(seurat_list[[1]], y = seurat_list[2:6])

# Remove mitochondrial genes
all_genes <- rownames(combined)
mt_genes <- grep("^mt-", all_genes, value = TRUE, ignore.case = TRUE)
combined <- subset(combined, features = setdiff(all_genes, mt_genes))

# Normalization, highly variable genes, and PCA
combined <- combined %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 50, verbose = FALSE)

# Add group information and perform batch correction with Harmony
combined$group <- ifelse(combined$orig.ident %in% c("c1", "c2", "c3"), "C_group", "E_group")

combined <- combined %>%
    RunHarmony(group.by.vars = "orig.ident", assay.use = "RNA", plot_convergence = FALSE)

# Dimensionality reduction and clustering (UMAP & t-SNE) and cluster identification
combined <- combined %>%
    RunUMAP(reduction = "harmony", dims = 1:10) %>%
    RunTSNE(reduction = "harmony", dims = 1:10) %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = 0.4)

# Visualize and save UMAP results
p_group <- DimPlot(combined, group.by = "group", cols = c("#99FFCC", "#FF99CC")) +
    ggtitle("Cells by Group (C vs E)")
ggsave("UMAP_by_group.png", plot = p_group, width = 8, height = 6, dpi = 300)

p_cluster <- DimPlot(combined, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle("Lymphocyte Clusters")
ggsave("UMAP_clusters.png", plot = p_cluster, width = 8, height = 6, dpi = 300)

# Find marker genes and export
combined <- JoinLayers(combined)
all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

top10_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top10_markers, file = "top10_markers_per_cluster.csv", row.names = FALSE, quote = FALSE)
write.csv(top20_markers, file = "top20_markers_per_cluster.csv", row.names = FALSE, quote = FALSE)

# Plot and save marker heatmap
top_genes <- unique(top10_markers$gene)
combined <- ScaleData(combined, features = top_genes)

p_heatmap <- DoHeatmap(combined, features = top_genes, group.by = "seurat_clusters") + 
    theme(axis.text.y = element_text(size = 6))
ggsave("top10_markers_heatmap.pdf", plot = p_heatmap, width = 16, height = 8)

# Save the final object
saveRDS(combined, "combined_final.rds")

# =====================================================================
# Project: Cell Type Annotation, Proportion Analysis & Differential Expression
# =====================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(randomcoloR)
combined <- readRDS("combined_final.rds")
# Fine-grained cell annotation
new.cluster.ids <- c("CD8 T cell", "CD3z T cell", "cd79a B cell", "ccr7 T cell", "Th17 cell", 
                     "CD4 T cell", "ano6 T cell", "M1 cell", "lta T cell", "MCM T cell",
                     "SOX4 cells", "SOX13 cells", "other1", "other2")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined$cell_type <- Idents(combined)

# DotPlot of core marker genes (bubble plot)
genes_to_plot <- c("ccl38.6", "cd8a", "cd247l", "cd79a", "ccr7", "il17a/f3", 
                   "cd4-1", "ano6", "lta", "mcm4", "sox4a", "sox13", "epcam")

p_dot <- DotPlot(combined, features = genes_to_plot, group.by = "cell_type", cols = c("lightgrey", "red"), dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 14)) +
    labs(x = "Cell Types", y = "Genes")

ggsave("DotPlot_markers.png", plot = p_dot, width = 12, height = 6, dpi = 300)

# Broad cell type categorization
combined$broad_celltype <- case_when(
    combined$cell_type %in% c("CD8 T cell", "CD3z T cell", "ccr7 T cell", "Th17 cell", 
                              "CD4 T cell", "ano6 T cell", "lta T cell", "MCM T cell") ~ "T Cell",
    combined$cell_type %in% c("cd79a B cell") ~ "B Cell",
    combined$cell_type %in% c("M1 cell") ~ "Macrophage",
    combined$cell_type %in% c("SOX13 cells", "SOX4 cells") ~ "SOX cells", 
    TRUE ~ "Others"
)

# Cell proportion statistics and stacked bar plot
prop_data <- combined@meta.data %>% 
    group_by(group, broad_celltype) %>% 
    summarise(n = n(), .groups = "drop") %>% 
    mutate(prop = n / sum(n))

p_prop <- ggplot(prop_data, aes(x = group, y = prop, fill = broad_celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8, color = "black", linewidth = 0.1) +
    coord_flip() +
    scale_fill_manual(values = c("T Cell" = "#1F77B4", "B Cell" = "#FF7F0E", 
                                 "Macrophage" = "#2CA02C", "SOX cells" = "#D62728", "Others" = "grey")) +
    labs(x = "Sample Group", y = "Proportion", fill = "Cell Type") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.line = element_line(color = "black", linewidth = 1),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          panel.grid = element_blank())

ggsave(filename = "Cell_Proportion_Barplot.png", plot = p_prop, width = 8, height = 4, dpi = 300)

# Differential expression analysis between groups (C_group vs E_group)
Idents(combined) <- "group"
de_genes <- FindMarkers(combined, ident.1 = "C_group", ident.2 = "E_group", 
                        test.use = "wilcox", logfc.threshold = 0.25, min.pct = 0.1)

de_genes_sorted <- de_genes[order(de_genes$p_val_adj), ]
write.csv(de_genes_sorted, "C_vs_E_differential_genes.csv")
