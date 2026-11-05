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
    scale_fill_manual(values = c("T cells" = "#377EB8", 
                                 "B cells" = "#66C2A5",  
                                 "DCs" = "#F781BF",     
                                 "others" = "#FFFF33")) +
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