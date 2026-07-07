#!/usr/bin/env Rscript
# re-run the main pipeline and save every figure as svg instead of png.
# same seeded steps as the Rmd, makes figs 01-13 and 18-23 in Plots/.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(svglite)
  library(scRNAseq); library(SingleCellExperiment); library(ggridges)
})
set.seed(42)
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n")
plots_dir <- "Plots"; dir.create(plots_dir, showWarnings = FALSE)

# save as svg (rewrite the .png name to .svg)
save_plot <- function(plot, filename, width = 10, height = 7) {
  filename <- sub("\\.png$", ".svg", filename)
  ggsave(file.path(plots_dir, filename), plot, width = width, height = height, device = svglite)
  msg("saved", filename)
}

msg("Loading Zilionis mouse data ...")
sce <- ZilionisLungData('mouse'); suppressWarnings(reducedDim(sce) <- NULL)
pbmcs <- as.Seurat(sce, counts = "counts", data = NULL)

# 01 raw cell distribution by tissue
cc <- as.data.frame(table(pbmcs$Tissue))
save_plot(ggplot(cc, aes(Var1, Freq, fill = Var1)) + geom_bar(stat="identity") +
            labs(x="Tissue Type", y="Number of cells", title="Cell Distribution by Tissue Type") +
            theme_minimal(), "01_cell_distribution_by_tissue.png", 8, 6)

pbmcs <- subset(pbmcs, Used == TRUE)
# 02 QC violin (pre-filter)
save_plot(VlnPlot(pbmcs, features = c("nFeature_originalexp","nCount_originalexp",
          "Percent.counts.from.mitochondrial.genes")), "02_QC_violin.png", 12, 6)

pbmcs <- subset(pbmcs, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 &
                  Percent.counts.from.mitochondrial.genes < 10)
# 03 QC violin (post-filter)
save_plot(VlnPlot(pbmcs, features = c("nFeature_originalexp","nCount_originalexp",
          "Percent.counts.from.mitochondrial.genes")), "03_QC_violin_postfilter.png", 12, 6)

# doublets
library(scDblFinder)
sce_db <- as.SingleCellExperiment(pbmcs, assay = "originalexp")
sce_db <- scDblFinder(sce_db, samples = pbmcs$Animal)
pbmcs$scDblFinder.class <- sce_db$scDblFinder.class
pbmcs$scDblFinder.score <- sce_db$scDblFinder.score
# 02b doublet score
save_plot(VlnPlot(pbmcs, features = "scDblFinder.score", group.by = "Tissue"),
          "02b_doublet_score.png", 8, 6)
pbmcs <- subset(pbmcs, subset = scDblFinder.class == "singlet")
msg("cells after doublet removal:", ncol(pbmcs))

pbmcs <- NormalizeData(pbmcs, normalization.method = "LogNormalize", verbose = FALSE)
pbmcs <- FindVariableFeatures(pbmcs, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# cell-cycle scoring
to_mouse <- function(g){ g <- tolower(g); substr(g,1,1) <- toupper(substr(g,1,1)); g }
s.genes   <- intersect(to_mouse(cc.genes$s.genes),   rownames(pbmcs))
g2m.genes <- intersect(to_mouse(cc.genes$g2m.genes), rownames(pbmcs))
pbmcs <- CellCycleScoring(pbmcs, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# 04 variable features
top22 <- head(VariableFeatures(pbmcs), 22)
save_plot(LabelPoints(plot = VariableFeaturePlot(pbmcs), points = top22, repel = TRUE),
          "04_variable_features.png", 10, 7)

pbmcs <- ScaleData(pbmcs, vars.to.regress = c("S.Score","G2M.Score"), verbose = FALSE)
pbmcs <- RunPCA(pbmcs, features = VariableFeatures(pbmcs), verbose = FALSE)
# 05 elbow, 06 pca dimplot
save_plot(ElbowPlot(pbmcs), "05_PCA_elbow.png", 8, 6)
save_plot(DimPlot(pbmcs, reduction = "pca", group.by = "Tissue"), "06_PCA_dimplot.png", 8, 6)

pbmcs <- RunUMAP(pbmcs, dims = 1:16, verbose = FALSE)
# 07,08,09
save_plot(DimPlot(pbmcs, group.by = "Animal", split.by = "Tissue"),
          "07_UMAP_by_animal_split_tissue.png", 14, 6)
save_plot(DimPlot(pbmcs, group.by = "Tissue"), "08_UMAP_by_tissue.png", 8, 6)
save_plot(FeaturePlot(pbmcs, features = "nFeature_originalexp"),
          "09_UMAP_featureplot_nFeature.png", 8, 6)

# harmony (for figs 10,11; does not affect PCA-based clustering below)
suppressPackageStartupMessages(library(harmony))
pbmcs <- RunHarmony(pbmcs, group.by.vars = "Animal")
save_plot(DimPlot(pbmcs), "10_harmony_dimplot.png", 8, 6)
save_plot(ElbowPlot(pbmcs), "11_harmony_elbow.png", 8, 6)

pbmcs <- FindNeighbors(pbmcs, dims = 1:16, verbose = FALSE)
pbmcs <- FindClusters(pbmcs, resolution = 0.4, verbose = FALSE)
pbmcs <- RunUMAP(pbmcs, dims = 1:16, verbose = FALSE)
msg("n clusters:", nlevels(Idents(pbmcs)))
save_plot(DimPlot(pbmcs, reduction = "umap", label = TRUE), "12_UMAP_clusters.png", 8, 6)

# 18-21 cell cycle diagnostics
save_plot(DimPlot(pbmcs, group.by = "Phase") + ggtitle("Cell-cycle phase on UMAP"),
          "18_cellcycle_UMAP_phase.png", 8, 6)
save_plot(FeaturePlot(pbmcs, features = c("S.Score","G2M.Score")),
          "19_cellcycle_scores_featureplot.png", 11, 5)
save_plot(RidgePlot(pbmcs, features = c("S.Score","G2M.Score"), group.by = "Phase", ncol = 2),
          "20_cellcycle_scores_ridge.png", 11, 5)
df <- as.data.frame(prop.table(table(Cluster = Idents(pbmcs), Phase = pbmcs$Phase), 1))
save_plot(ggplot(df, aes(Cluster, Freq, fill = Phase)) + geom_col() +
            labs(y="Proportion of cells", title="Cell-cycle phase composition per cluster") +
            theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
          "21_cellcycle_phase_by_cluster.png", 9, 6)

# 13 cluster-marker heatmap (raster body -> compact SVG)
msg("FindAllMarkers (slow) ...")
markers <- FindAllMarkers(pbmcs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
top10 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
save_plot(DoHeatmap(pbmcs, features = top10$gene, raster = TRUE) + ggtitle("Cluster Markers HeatMap"),
          "13_cluster_markers_heatmap.png", 12, 10)

# annotation (positional labels, same as Rmd) -> figs 22, 23
if (nlevels(Idents(pbmcs)) != 12) msg("WARNING: expected 12 clusters, got", nlevels(Idents(pbmcs)))
new.ids <- c("0"="Neutrophils (mature)","1"="B cells","2"="T cells",
             "3"="Neutrophils (Siglecf+ TAN)","4"="Monocytes","5"="NK cells",
             "6"="Neutrophils (immature)","7"="Neutrophils (Ccl3+)","8"="Alveolar macrophages",
             "9"="cDC2","10"="mregDC","11"="Macrophages (C1q+Trem2+ TAM)")
pbmcs <- RenameIdents(pbmcs, new.ids)
save_plot(DimPlot(pbmcs, reduction = "umap", label = TRUE, repel = TRUE) +
            ggtitle("Manually annotated cell types (top-5 marker based)"),
          "22_UMAP_annotated_manual.png", 11, 7)
canonical <- c("Cd79a","Ms4a1","Cd19","Cd3g","Lat","Thy1","Ncr1","Gzma","Klra4",
               "Ace","F13a1","Clec4a1","Krt79","Cidec","C1qa","C1qb","Trem2","Mmp12",
               "Cd209a","Clec10a","Ccl17","Fscn1","Ccl22",
               "S100a8","Ly6g","Retnlg","Ngp","Mmp8","Ccl3","Siglecf","Ffar2")
canonical <- canonical[canonical %in% rownames(pbmcs)]
save_plot(DotPlot(pbmcs, features = canonical) + RotatedAxis() +
            ggtitle("Canonical marker expression per annotated cell type"),
          "23_canonical_marker_dotplot.png", 15, 7)

msg("DONE regen_pipeline_svg")
