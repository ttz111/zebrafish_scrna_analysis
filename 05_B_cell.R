# =====================================================================
# Project: B Cell Subcluster, Trajectory & KEGG Enrichment Analysis
# =====================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle)
library(clusterProfiler)
library(org.Dr.eg.db)

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
# Core figure 1: UMAP and FeaturePlot
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

# ---------------------------------------------------------------------
# Core figure 2: Marker gene DotPlot
# ---------------------------------------------------------------------
gene_list <- c("socs3b", "tagln2", "ncf2", "npm1a", "fbl", "gnl3", "egr1", "fosb", "nr4a1",
               "ifng1r", "faslg", "il2rb", "ddit4", "hsp70l", "rps26", "cd79a", "ighd", "mhc1zba",
               "foxp3a", "ikzf4", "cd27", "fkbp11", "ighv1-4.1", "tent5c", "spi1a", "cfb", "nrg1")

pb <- DotPlot(B_cells_subset_data, features = gene_list, group.by = "seurat_clusters", cols = c("lightgrey", "red"), dot.scale = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16)) +
    labs(x = "Genes", y = "Cell Types")
ggsave("b_cell_dotplot.png", plot = pb, width = 10, height = 5, dpi = 300)

# ---------------------------------------------------------------------
# Core figure 3: Monocle2 pseudotime developmental trajectory analysis
# ---------------------------------------------------------------------
data <- LayerData(B_cells_subset_data, assay = "RNA", layer = "counts")
pd <- new('AnnotatedDataFrame', data = B_cells_subset_data@meta.data)
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = rownames(data), row.names = rownames(data)))

mycds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds)
mycds <- detectGenes(mycds, min_expr = 0.1)

mycds_expressed_genes <- row.names(subset(fData(mycds), num_cells_expressed >= 10))
mycds <- mycds[mycds_expressed_genes, ]

diff_test_res <- differentialGeneTest(mycds[mycds_expressed_genes, ], fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
mycds <- setOrderingFilter(mycds, ordering_genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)

p1 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") + scale_color_manual(values = my_colors)
p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("trajectory_cell_type.png", plot = p1, width = 4, height = 5, dpi = 300)
ggsave("trajectory_pseudotime.png", plot = p2, width = 4, height = 4, dpi = 300)

# Pseudotime expression plot of a key gene
cds_subset <- mycds["ighv1-4.1", ]
p_gene_pseudo <- plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters", ncol = 1) +
    scale_color_manual(values = my_colors) + ggtitle("ighv1-4.1 Expression along Pseudotime")
ggsave("ighm_pseudotime.png", plot = p_gene_pseudo, width = 5, height = 4, dpi = 300)

# ---------------------------------------------------------------------
# Core figure 4: Between-group proportion analysis of B cell subpopulations (Composition)
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

# ---------------------------------------------------------------------
# Core figure 5: KEGG bubble enrichment plot (Top 5 pathways per cluster)
# ---------------------------------------------------------------------
B_cells_subset_data <- JoinLayers(B_cells_subset_data)
DEGs <- FindAllMarkers(B_cells_subset_data, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = "wilcox")

top_genes <- DEGs %>% group_by(cluster) %>% filter(p_val_adj < 0.05) %>% top_n(n = 200, wt = avg_log2FC)
gene_entrez <- bitr(top_genes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)

kegg_results <- lapply(split(top_genes, top_genes$cluster), function(df) {
    entrez_ids <- gene_entrez$ENTREZID[gene_entrez$SYMBOL %in% df$gene]
    enrichKEGG(gene = entrez_ids, organism = "dre", pvalueCutoff = 0.05, qvalueCutoff = 0.1, keyType = "kegg")
})

kegg_df <- purrr::map2_dfr(kegg_results[!sapply(kegg_results, is.null)], names(kegg_results)[!sapply(kegg_results, is.null)], function(enrich_obj, clust_name) {
    as.data.frame(enrich_obj) %>% mutate(cluster = clust_name) %>% filter(p.adjust < 0.1)
})

kegg_top <- kegg_df %>% group_by(cluster) %>% arrange(p.adjust) %>% slice_head(n = 5) %>% ungroup()

p4 <- ggplot(kegg_top, aes(x = cluster, y = reorder(Description, -p.adjust), size = Count, color = -log10(p.adjust))) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() +
    labs(title = "Top 5 KEGG Pathways per Cluster", x = "Cell Cluster", y = "KEGG Pathway", color = "-log10(p.adj)", size = "Gene Count") +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
ggsave("KEGG_top5.png", plot = p4, width = 9, height = 9, dpi = 300)
