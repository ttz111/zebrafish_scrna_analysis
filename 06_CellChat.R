# =====================================================================
# Project: CellChat Immune Cell-Cell Communication Analysis
# =====================================================================

library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(tidyverse)

setwd("single cell/new_anlysis/20260603/cellchat")

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

# ---------------------------------------------------------------------
# Core figure 1: Global cell communication circle plot (Circle Plot - Count & Strength)
# ---------------------------------------------------------------------
png("Zebrafish_Immune_Circle_Count.png", width = 1500, height = 1500, res = 300)
netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)), 
                 weight.scale = TRUE, label.edge = FALSE, vertex.label.cex = 0.00001, 
                 title.name = "Number of interactions")
dev.off()

png("Zebrafish_Immune_Circle_Strength.png", width = 1500, height = 1500, res = 300)
netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)), 
                 weight.scale = TRUE, label.edge = FALSE, vertex.label.cex = 0.00001, 
                 title.name = "Interaction strength")
dev.off()

# ---------------------------------------------------------------------
# Core figure 2: Communication bubble plot from T cells (Bubble Plot)
# ---------------------------------------------------------------------
p_T <- netVisual_bubble(cellchat, 
                        sources.use = c("T1", "T2"), 
                        targets.use = c("B0", "B1", "B2", "B3", "B4", "B6"), 
                        remove.isolate = FALSE)

ggsave("Zebrafish_Immune_Bubble_Plot.png", plot = p_T, width = 7, height = 4, dpi = 300)

# ---------------------------------------------------------------------
# Core figure 3: Customized advanced bubble plot based on actual p-values
# ---------------------------------------------------------------------
bubble_data <- p_T$data
net_df <- subsetCommunication(cellchat) 

merged_bubble_data <- merge(
    bubble_data, 
    net_df[, c("source", "target", "ligand", "receptor", "pval")], 
    by = c("source", "target", "ligand", "receptor"),
    all.x = TRUE,
    suffixes = c(".plot", ".original")
)

merged_bubble_data$pair <- paste0(merged_bubble_data$source, " -> ", merged_bubble_data$target)
merged_bubble_data$lr_pair <- merged_bubble_data$interaction_name_2
merged_bubble_data$p_group <- factor(ifelse(merged_bubble_data$pval.original < 0.05, "< 0.05", "> 0.05"), levels = c("> 0.05", "< 0.05"))

p_bubble <- ggplot(merged_bubble_data, aes(x = pair, y = lr_pair)) +
    geom_point(aes(size = p_group, color = prob)) +
    scale_color_gradientn(
        colors = c("#3A54A4", "#499FA4", "#F7ECC7", "#EE7648", "#A41646"), 
        name = "Commum.Prob"
    ) +
    scale_size_manual(
        values = c("> 0.05" = 1.5, "< 0.05" = 5), 
        name = "p-value",
        drop = FALSE
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title = element_blank(),
        panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.title = element_text(size = 11, face = "plain"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(),
        legend.key = element_blank()
    )

ggsave("CellChat_Custom_Bubble.png", plot = p_bubble, width = 9, height = 3.5, dpi = 300)