
# =====================================================================
# Project: CD8+ T Cell Subcluster, Trajectory & KEGG Enrichment Analysis
# =====================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle)
library(clusterProfiler)
library(org.Dr.eg.db)

# Extract the CD8+ cell subset and filter out CD4-positive contaminating cells
seurat_obj <- readRDS("single cell/new_anlysis/20251201/tcell/tcell_clean.rds")
cd8_subset <- subset(seurat_obj, subset = `cd8a` > 0)
cd8_cleaned <- subset(cd8_subset, subset = `cd4-1` <= 0.1)

# Re-run the Seurat dimensionality reduction and clustering workflow
cd8_cleaned <- FindVariableFeatures(cd8_cleaned, selection.method = "vst", nfeatures = 2000)
cd8_cleaned <- ScaleData(cd8_cleaned)
cd8_cleaned <- RunPCA(cd8_cleaned, npcs = 30)
cd8_cleaned <- FindNeighbors(cd8_cleaned, dims = 1:15)
cd8_cleaned <- FindClusters(cd8_cleaned, resolution = 0.1)
cd8_cleaned <- RunUMAP(cd8_cleaned, dims = 1:15)
cd8_cleaned <- JoinLayers(cd8_cleaned)

# Extract core clusters and perform cell type annotation
cd8_final <- subset(cd8_cleaned, idents = c("0", "1", "2", "3"))
cd8_final <- FindVariableFeatures(cd8_final, selection.method = "vst", nfeatures = 2000)
cd8_final <- ScaleData(cd8_final)
cd8_final <- RunPCA(cd8_final, npcs = 30)
cd8_final <- FindNeighbors(cd8_final, dims = 1:15)
cd8_final <- FindClusters(cd8_final, resolution = 0.1)
cd8_final <- RunUMAP(cd8_final, dims = 1:15)
cd8_final <- JoinLayers(cd8_final)

new_cluster_labels <- c(
    "0" = "Resident_CD8_T",
    "1" = "Regulatory_CD8_T",
    "2" = "Naive CD8_T",
    "3" = "Effector_CD8_T"
)
cd8_final <- RenameIdents(cd8_final, new_cluster_labels)
cd8_final$cell_type <- Idents(cd8_final)

# Save the annotated object and the marker gene table
saveRDS(cd8_final, file = "cd8_final.rds")
all_markers <- FindAllMarkers(cd8_final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")
top20_markers <- all_markers %>% group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)
write.csv(top20_markers, "cd8_top20_markers.csv", row.names = FALSE)

# ---------------------------------------------------------------------
# Core figure 1: UMAP distribution plot
# ---------------------------------------------------------------------
p_umap <- DimPlot(cd8_final, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 1) +
    theme(axis.title = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12), 
          axis.line = element_line(linewidth = 1.2, color = "black"), 
          axis.ticks = element_line(size = 1.2),
          legend.text = element_text(size = 12))
ggsave("CD8_T_Cell_UMAP.png", plot = p_umap, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------
# Core figure 2: Feature marker bubble plot (DotPlot)
# ---------------------------------------------------------------------
features_to_plot <- c(
    "ccr7", "tcf7", "ccr9a", "cxcr4a", 
    "il17a/f3", "foxp3a", "rorc", "il23r", 
    "eomesb", "gzm3", "ifng1", "gzma", 
    "ccl38.6", "ccl36.1", "ccl38a.5", "plxnc1"
)

p_dot <- DotPlot(cd8_final, features = features_to_plot, group.by = "cell_type", dot.scale = 8, cols = c("lightgrey", "red")) +
    RotatedAxis() +
    labs(title = "Marker Expression across CD8+ T Cell Subsets", x = "Feature Genes", y = "Cell Type") +
    theme(
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, face = "italic"),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(fill = NA, linewidth = 1.2)
    )
ggsave("Markercd8_DotPlot.png", plot = p_dot, width = 8, height = 5, dpi = 300)

# ---------------------------------------------------------------------
# Core figure 3: Cell proportion donut plot (Control vs Infected)
# ---------------------------------------------------------------------
plot_data <- cd8_final@meta.data %>%
    group_by(group, cell_type) %>%
    tally() %>%
    group_by(group) %>%
    mutate(pct = n / sum(n), 
           ymax = cumsum(pct), 
           ymin = c(0, head(ymax, n = -1)), 
           label_pos = (ymax + ymin) / 2, 
           label_text = paste0(round(pct * 100, 1), "%"))

my_colors <- c("#F8766D", "#0cb702", "#00bfc4", "#ff61cc")

p_prop <- ggplot(plot_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = cell_type)) +
    geom_rect(color = "white") +
    geom_text(x = 4.5, aes(y = label_pos, label = label_text), size = 3.5) +
    coord_polar(theta = "y") +
    xlim(c(2.5, 4.5)) +
    facet_wrap(~group) +
    theme_void() +
    theme(legend.position = "right", strip.text = element_text(size = 14, face = "bold")) + 
    scale_fill_manual(values = my_colors) +
    labs(title = "CD8+ Cell Type Proportion: C vs E Group", fill = "Cell Type")
ggsave("CD8_Cell_Proportion_Donut.png", plot = p_prop, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------
# Core figure 4: Monocle2 pseudotime developmental trajectory
# ---------------------------------------------------------------------
data <- GetAssayData(cd8_final, assay = "RNA", slot = "counts")
pd <- new('AnnotatedDataFrame', data = cd8_final@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = row.names(data), row.names = row.names(data)))
my_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)

diff_test_res <- differentialGeneTest(my_cds, fullModelFormulaStr = "~cell_type", cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
my_cds <- setOrderingFilter(my_cds, ordering_genes)
my_cds <- reduceDimension(my_cds, max_components = 2, method = 'DDRTree')
my_cds <- orderCells(my_cds)

p_traj <- plot_cell_trajectory(my_cds, color_by = "cell_type") + scale_color_manual(values = my_colors)
ggsave("CD8_Trajectory.png", plot = p_traj, width = 7, height = 6, dpi = 300)

# ---------------------------------------------------------------------
# Core figure 5: Mirror bar plot of KEGG enrichment for between-group DEGs
# ---------------------------------------------------------------------
Idents(cd8_final) <- cd8_final$group
cd8_degs <- FindMarkers(object = cd8_final, ident.1 = "E_group", ident.2 = "C_group", test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.25)

up_genes <- rownames(subset(cd8_degs, avg_log2FC > 0.25 & p_val_adj < 0.05))
down_genes <- rownames(subset(cd8_degs, avg_log2FC < -0.25 & p_val_adj < 0.05))

up_entrez <- bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dr.eg.db, drop=TRUE)$ENTREZID
down_entrez <- bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Dr.eg.db, drop=TRUE)$ENTREZID

kegg_up <- enrichKEGG(gene = up_entrez, organism = "dre", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_entrez, organism = "dre", pvalueCutoff = 0.05)

df_up <- as.data.frame(kegg_up) %>% top_n(8, wt = -pvalue) %>% mutate(Group = "Infected (E_group)", Direction = 1, LogP = -log10(pvalue))
df_down <- as.data.frame(kegg_down) %>% top_n(2, wt = -pvalue) %>% mutate(Group = "Control (C_group)", Direction = -1, LogP = -log10(pvalue))

duplicate_pathways <- intersect(df_up$Description, df_down$Description)
df_plot <- rbind(
    df_up %>% filter(!Description %in% duplicate_pathways),
    df_down %>% filter(!Description %in% duplicate_pathways)
) %>% 
    mutate(PlotValue = LogP * Direction) %>% 
    arrange(PlotValue)
df_plot$Description <- factor(df_plot$Description, levels = unique(df_plot$Description))

sci_colors <- c("Control (C_group)" = "#1F78B4", "Infected (E_group)" = "#E31A1C")
p_cd8_kegg <- ggplot(df_plot, aes(x = Description, y = PlotValue, fill = Group)) +
    geom_bar(stat = "identity", width = 0.65, color = "black", size = 0.2) +
    coord_flip() +
    scale_y_continuous(labels = function(x) abs(x), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = sci_colors) +
    theme_classic() +
    labs(x = "Enriched KEGG Pathways", y = "-log10(p-value)", fill = "Conditions", title = "Functional Transition of CD8+ T Cells Post-Infection") +
    theme(
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", color = "black"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 13, face = "bold", color = "#333333"),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        panel.grid = element_blank()
    )

ggsave("cd8_kegg.png", plot = p_cd8_kegg, width = 8.5, height = 6, dpi = 300)