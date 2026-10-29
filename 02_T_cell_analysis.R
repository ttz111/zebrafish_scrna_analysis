# =====================================================================
# Project: T-Cell Subcluster, Trajectory & DEG/KEGG Enrichment Analysis
# =====================================================================

# Environment initialization and data loading
rm(list = ls()); gc()
setwd("single cell/new_anlysis/20251201/tcell")

library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Dr.eg.db)
library(monocle)
library(ggrepel)
library(pheatmap)

T_cell <- readRDS("tcell_clean.rds")
T_cell <- JoinLayers(T_cell)

# Cluster re-integration and cell type annotation
cluster_mapping <- c("0"="0", "1"="1", "2"="2", "3"="3", "4"="4",
                     "5"="5", "6"="6", "7"="7", "8"="8", "9"="1", "10"="1", "11"="9")
T_cell$seurat_clusters <- factor(as.numeric(cluster_mapping[as.character(T_cell$seurat_clusters)]), levels = 0:9)
Idents(T_cell) <- factor(Idents(T_cell), levels = as.character(0:9))

new.cluster.ids <- c("NK cell", "Naïve CD8⁺ T cell", "ccr7+ T cell", "rorc+ T cell",
                     "mcm2+ T cell", "CD4+ T cellls", "Effector T cell",
                     "Treg cell", "lta+ T cell", "sema7a+ T cell")
names(new.cluster.ids) <- levels(T_cell)
T_cell <- RenameIdents(T_cell, new.cluster.ids)
T_cell$cell_type <- Idents(T_cell)

my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#A65628",
               "#FF7F00", "#FFFF33", "#8DD3C7", "#F781BF", "#3300FF")

# Core expression and composition visualizations
# UMAP dimensionality reduction plot
p_umap <- DimPlot(T_cell, label = FALSE, cols = my_colors) + ggtitle("T Cell Subclusters") + theme_bw()
ggsave("T_umap.png", plot = p_umap, width = 8, height = 6, dpi = 300)

# Marker gene DotPlot
genes_to_plot <- c("nitr12", "cd8a", "ccr7", "rorc", "mcm2", "cd4-1", "gzma", "foxp3a", "lta", "sema7a")
p_dot <- DotPlot(T_cell, features = genes_to_plot, group.by = "cell_type", cols = c("lightgrey", "red"), dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("T_cell_dotplot.png", plot = p_dot, width = 10, height = 5, dpi = 300)

# Cell proportion analysis
cell_prop_df <- T_cell@meta.data %>%
    group_by(group, cell_type) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(group) %>%
    mutate(Proportion = Count / sum(Count) * 100)

p_prop <- ggplot(cell_prop_df, aes(x = group, y = Proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "fill", width = 0.4) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = my_colors) +
    theme_minimal() + labs(x = "Group", y = "Proportion")
ggsave("T_cell_proportion.png", plot = p_prop, width = 6, height = 6, dpi = 300)

# Pseudotime trajectory analysis (Monocle2)
counts_matrix <- GetAssayData(T_cell, assay = "RNA", layer = "counts")
pd <- new('AnnotatedDataFrame', data = T_cell@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(counts_matrix), row.names = rownames(counts_matrix)))

mycds <- newCellDataSet(as(counts_matrix, "sparseMatrix"), phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds)
mycds <- detectGenes(mycds, min_expr = 0.1)

diff_test_res <- differentialGeneTest(mycds, fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
mycds <- setOrderingFilter(mycds, ordering_genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)

p_traj <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") + 
    facet_wrap(~group) + theme_bw() + scale_color_manual(values = my_colors)
ggsave("T_trajectory_by_group.png", plot = p_traj, width = 11, height = 5, dpi = 300)

# Differential gene expression analysis between groups and volcano plot (E_group vs C_group)
Idents(T_cell) <- "group"
deg_all <- FindMarkers(T_cell, ident.1 = "E_group", ident.2 = "C_group", logfc.threshold = 0, test.use = "wilcox")
deg_all$gene <- rownames(deg_all)
deg_all$label <- ifelse(deg_all$p_val_adj < 0.05 & deg_all$avg_log2FC > 0.25, "Up",
                        ifelse(deg_all$p_val_adj < 0.05 & deg_all$avg_log2FC < -0.25, "Down", "Stable"))

select_genes <- c("ly6m2", "rgs2", "mhc1zca", "tapbp.2", "ftr32.1", "ncf1", "zbtb16a", "sox4b")
top_deg <- deg_all[deg_all$gene %in% select_genes, ]

p_volcano <- ggplot(deg_all, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = label), alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Stable" = "grey80", "Down" = "#377EB8")) +
    geom_text_repel(data = top_deg, aes(label = gene)) +
    geom_vline(xintercept = c(-0.25, 0.25), lty = "dashed") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    theme_bw() + labs(title = "Volcano Plot: E_group vs C_group", x = "log2 Fold Change", y = "-log10 Adj P-value")
ggsave("Volcano_T_E_vs_C.png", plot = p_volcano, width = 5, height = 4, dpi = 300)

# KEGG functional enrichment analysis of Effector T cells
effector_cells <- subset(T_cell, subset = cell_type == "Effector T cell")
Idents(effector_cells) <- "group"
eff_markers <- FindMarkers(effector_cells, ident.1 = "E_group", ident.2 = "C_group", min.pct = 0.25, logfc.threshold = 0.25)

top_genes <- rownames(eff_markers %>% arrange(p_val_adj) %>% head(200))
gene_ids <- bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)

kegg_res <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "dre", pvalueCutoff = 0.05)
p_kegg <- dotplot(kegg_res, showCategory = 10) + ggtitle("Effector T cell: KEGG Enrichment (E vs C)")
ggsave("KEGG_Effector_Tcell.png", plot = p_kegg, width = 6.5, height = 6, dpi = 300)

# Save key data object
saveRDS(T_cell, "T_cell_final.rds")