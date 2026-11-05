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

# Load data and rename cell subpopulations
setwd("single cell/new_anlysis/20251201/bcell")
B_cells_subset_data <- readRDS("single cell/mono/B_cells/B_cells_clustered.rds")

new.cluster.ids <- c("scos3b+ B cell", "proliferating B cell", "active B cell",
                     "ccl38+ B cell", "hsp70+ B cell", "IgDhi B cell", 
                     "foxp3a+ B cell", "IgMhi B cell", "other")
names(new.cluster.ids) <- levels(B_cells_subset_data)
B_cells_subset_data <- RenameIdents(B_cells_subset_data, new.cluster.ids)
B_cells_subset_data$cell_identity <- Idents(B_cells_subset_data)

# ---------------------------------------------------------------------
#  UMAP and FeaturePlot
# ---------------------------------------------------------------------
p_clean <- DimPlot(B_cells_subset_data, cols = my_colors, group.by = "seurat_clusters", label = TRUE, label.size = 0) + 
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(linewidth = 1),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 12, face = "bold")) + 
    ggtitle("B cell")
ggsave("B_umap.png", plot = p_clean, width = 6, height = 5, dpi = 300)

p6 <- FeaturePlot(B_cells_subset_data, features = c("cd79a", "ighd", "ighz", "ighv1-4.1", "socs3b", "dkc1", "mpeg1.1", "ifng1r", "foxp3a"))
ggsave("features_B.png", plot = p6, width = 9, height = 8, dpi = 300)

#  Between-group proportion analysis of B cell subpopulations (Composition)
# ---------------------------------------------------------------------
cell_prop_df <- as.data.frame(table(B_cells_subset_data$group, B_cells_subset_data$cell_identity)) %>%
    rename(Group = Var1, Cell_type = Var2, Count = Freq) %>%
    group_by(Group) %>%
    mutate(Proportion = Count / sum(Count) * 100) %>%
    ungroup()

p3 <- ggplot(cell_prop_df, aes(x = Group, y = Proportion, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "fill", width = 0.4) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = my_colors) +
    labs(title = "Cell Type Composition by Group", x = "Group", y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")
ggsave("B_cell_prop.png", plot = p3, width = 6, height = 6, dpi = 300)

# Extract the CD4+ cell subset and filter out contaminating cells
seurat_obj <- readRDS("single cell/new_anlysis/20251201/tcell/tcell_clean.rds")
cd4_subset <- subset(seurat_obj, subset = `cd4-1` > 0)
cd4_cleaned <- subset(cd4_subset, subset = cd8a <= 0.5)

# Re-run the standard Seurat dimensionality reduction and clustering workflow
cd4_cleaned <- FindVariableFeatures(cd4_cleaned, selection.method = "vst", nfeatures = 2000)
cd4_cleaned <- ScaleData(cd4_cleaned)
cd4_cleaned <- RunPCA(cd4_cleaned, npcs = 30)
cd4_cleaned <- FindNeighbors(cd4_cleaned, dims = 1:15)
cd4_cleaned <- FindClusters(cd4_cleaned, resolution = 0.2)
cd4_cleaned <- RunUMAP(cd4_cleaned, dims = 1:15)
cd4_cleaned <- JoinLayers(cd4_cleaned)

# Extract core clusters and rename cell types
cd4_final <- subset(cd4_cleaned, idents = c("0", "1", "2", "3", "4"))
cd4_final <- FindVariableFeatures(cd4_final, selection.method = "vst", nfeatures = 2000)
cd4_final <- ScaleData(cd4_final)
cd4_final <- RunPCA(cd4_final, npcs = 30)
cd4_final <- FindNeighbors(cd4_final, dims = 1:15)
cd4_final <- FindClusters(cd4_final, resolution = 0.2)
cd4_final <- RunUMAP(cd4_final, dims = 1:15)
cd4_final <- JoinLayers(cd4_final)

new_cluster_labels <- c(
    "0" = "Activated_T", "1" = "Th1", "2" = "T_Effector",
    "3" = "Th17_Th22", "4" = "Cycling_T", "5" = "Innate_like_T"
)
cd4_final <- RenameIdents(cd4_final, new_cluster_labels)
cd4_final$cell_type <- Idents(cd4_final)

# Core figure 1: UMAP distribution plot
p_umap <- DimPlot(cd4_final, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 1) +
    theme(axis.title = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12), 
          axis.line = element_line(linewidth = 1.2, color = "black"),
          legend.text = element_text(size = 12))
ggsave(filename = "CD4_T_Cell_UMAP.png", plot = p_umap, width = 8, height = 6, dpi = 300)

# Core figure 3: Cell proportion donut/stacked bar plot (C vs E Group)
plot_data <- cd4_final@meta.data %>%
    group_by(group, cell_type) %>%
    tally() %>%
    group_by(group) %>%
    mutate(pct = n / sum(n), ymax = cumsum(pct), ymin = c(0, head(ymax, n = -1)), label_pos = (ymax + ymin) / 2, label_text = paste0(round(pct * 100, 1), "%"))

my_colors <- c("#F8766D", "#cd9600", "#0cb702", "#00bfc4", "#00a9ff", "#ff61cc")
p_prop <- ggplot(plot_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = cell_type)) +
    geom_rect(color = "white") +
    geom_text(x = 4.5, aes(y = label_pos, label = label_text), size = 3.5) +
    coord_polar(theta = "y") +
    xlim(c(2.5, 4.5)) +
    facet_wrap(~group) +
    theme_void() +
    theme(legend.position = "right", strip.text = element_text(size = 14, face = "bold")) +
    scale_fill_manual(values = my_colors) +
    labs(title = "Cell Type Proportion: C vs E Group", fill = "Cell Type")
ggsave("CD4_Cell_Proportion.png", plot = p_prop, width = 8, height = 6, dpi = 300)


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
#UMAP distribution plot
# ---------------------------------------------------------------------
p_umap <- DimPlot(cd8_final, reduction = "umap", group.by = "cell_type", label = FALSE, pt.size = 1) +
    theme(axis.title = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12), 
          axis.line = element_line(linewidth = 1.2, color = "black"), 
          axis.ticks = element_line(size = 1.2),
          legend.text = element_text(size = 12))
ggsave("CD8_T_Cell_UMAP.png", plot = p_umap, width = 8, height = 6, dpi = 300)

# Cell proportion donut plot (Control vs Infected)
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

# Load T cell and B cell data and remap the layout
B_cell <- readRDS("single cell/new_anlysis/20251201/bcell/new/B_final_cell.rds")
T_cell <- readRDS("single cell/new_anlysis/20251201/tcell/tcell_clean.rds")

# T cell subcluster renaming and classification mapping
t_digit_mapping <- c(
    "Th17" = "T1", "T-cd4" = "T1", "T-foxp3a" = "T1",
    "T-cd8" = "T2", "T-gzma" = "T2", "T-hsp70" = "T2", "T-itgae.2" = "T2",
    "NK cells" = "T3", "T-ccr7" = "T4", "T-mcm2" = "T5", "T-lta" = "T6", "T-sema7a" = "T7"
)
T_cell@meta.data[["cellchat_layout"]] <- t_digit_mapping[as.character(Idents(T_cell))]
B_cell$cellchat_layout <- paste0("B", B_cell$seurat_clusters)

# Merge immune cell objects and create the CellChat object
immune_merged <- merge(T_cell, y = B_cell, add.cell.ids = c("T", "B"))
immune_merged <- JoinLayers(immune_merged)
Idents(immune_merged) <- "cellchat_layout"
immune_merged$ident <- Idents(immune_merged)

cellchat <- createCellChat(object = immune_merged, group.by = "ident")
cellchat@DB <- CellChatDB.zebrafish

# Run the core CellChat computation workflow
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Communication bubble plot from T cells (Bubble Plot)
# ---------------------------------------------------------------------
p_T <- netVisual_bubble(cellchat, 
                        sources.use = c("T1", "T2"), 
                        targets.use = c("B0", "B1", "B2", "B3", "B4", "B6"), 
                        remove.isolate = FALSE)

ggsave("Zebrafish_Immune_Bubble_Plot.png", plot = p_T, width = 7, height = 4, dpi = 300)