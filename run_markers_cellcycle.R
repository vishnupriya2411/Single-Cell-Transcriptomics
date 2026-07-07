#!/usr/bin/env Rscript
# rebuild the clustered object, export markers and make the cell-cycle figures.
# same as the Rmd (PCA, dims 1:16, res 0.4) just seeded and scaling only the HVGs
# to save memory - clusters come out the same.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
  library(scRNAseq); library(SingleCellExperiment); library(ggridges)
})
set.seed(42)
t0 <- Sys.time(); msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(),"%H:%M:%S")), ..., "\n")

plots_dir <- "Plots"; dir.create(plots_dir, showWarnings = FALSE)
save_plot <- function(p, f, w=10, h=7, dpi=200)
  ggsave(file.path(plots_dir, f), p, width=w, height=h, dpi=dpi)

msg("Loading Zilionis mouse data ...")
sce <- ZilionisLungData('mouse')
suppressWarnings(reducedDim(sce) <- NULL)
pbmcs <- as.Seurat(sce, counts = "counts", data = NULL)
pbmcs <- subset(pbmcs, Used == TRUE)
msg("cells after Used==TRUE:", ncol(pbmcs))

# QC filter (same thresholds as Rmd)
pbmcs <- subset(pbmcs, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 &
                  Percent.counts.from.mitochondrial.genes < 10)
msg("cells after QC filter:", ncol(pbmcs))

# Doublet removal (best-effort)
tryCatch({
  library(scDblFinder)
  sce_db <- as.SingleCellExperiment(pbmcs, assay = "originalexp")
  sce_db <- scDblFinder(sce_db, samples = pbmcs$Animal)
  pbmcs$scDblFinder.class <- sce_db$scDblFinder.class
  pbmcs <- subset(pbmcs, subset = scDblFinder.class == "singlet")
  msg("cells after doublet removal:", ncol(pbmcs))
}, error = function(e) msg("scDblFinder skipped:", conditionMessage(e)))

# Normalize + HVG
pbmcs <- NormalizeData(pbmcs, normalization.method = "LogNormalize", verbose = FALSE)
pbmcs <- FindVariableFeatures(pbmcs, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# cell cycle scoring (mouse-cased gene names)
to_mouse <- function(g){ g <- tolower(g); substr(g,1,1) <- toupper(substr(g,1,1)); g }
s.genes   <- intersect(to_mouse(cc.genes$s.genes),   rownames(pbmcs))
g2m.genes <- intersect(to_mouse(cc.genes$g2m.genes), rownames(pbmcs))
pbmcs <- CellCycleScoring(pbmcs, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
msg("Phase table:"); print(table(pbmcs$Phase))

# Scale (regress cell cycle) -> PCA -> UMAP -> cluster
pbmcs <- ScaleData(pbmcs, vars.to.regress = c("S.Score","G2M.Score"), verbose = FALSE)
pbmcs <- RunPCA(pbmcs, features = VariableFeatures(pbmcs), verbose = FALSE)
pbmcs <- RunUMAP(pbmcs, dims = 1:16, verbose = FALSE)
pbmcs <- FindNeighbors(pbmcs, dims = 1:16, verbose = FALSE)
pbmcs <- FindClusters(pbmcs, resolution = 0.4, verbose = FALSE)
msg("n clusters:", nlevels(Idents(pbmcs)))
print(table(Idents(pbmcs)))

# cell cycle figures
msg("Writing cell-cycle figures ...")
# UMAP colored by Phase
save_plot(DimPlot(pbmcs, group.by = "Phase") + ggtitle("Cell-cycle phase on UMAP"),
          "18_cellcycle_UMAP_phase.png", 8, 6)
# FeaturePlots of the two scores
save_plot(FeaturePlot(pbmcs, features = c("S.Score","G2M.Score")),
          "19_cellcycle_scores_featureplot.png", 11, 5)
# Ridge plot of the two scores by phase
rp <- RidgePlot(pbmcs, features = c("S.Score","G2M.Score"), group.by = "Phase", ncol = 2)
save_plot(rp, "20_cellcycle_scores_ridge.png", 11, 5)
# Stacked bar: phase proportion per cluster
df <- as.data.frame(prop.table(table(Cluster = Idents(pbmcs), Phase = pbmcs$Phase), 1))
p_bar <- ggplot(df, aes(Cluster, Freq, fill = Phase)) +
  geom_col() + labs(y = "Proportion of cells", title = "Cell-cycle phase composition per cluster") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_bar, "21_cellcycle_phase_by_cluster.png", 9, 6)

# markers (slow step)
msg("FindAllMarkers ... (this is the slow step)")
markers <- FindAllMarkers(pbmcs, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(plots_dir, "cluster_markers_all.csv"), row.names = FALSE)

top10 <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
top5  <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
write.csv(top10, file.path(plots_dir, "cluster_markers_top10.csv"), row.names = FALSE)
write.csv(top5,  file.path(plots_dir, "cluster_markers_top5.csv"),  row.names = FALSE)

cat("\n=================  TOP 5 MARKERS PER CLUSTER  =================\n")
for (cl in levels(Idents(pbmcs))) {
  g <- top5$gene[top5$cluster == cl]
  cat(sprintf("Cluster %s (n=%d): %s\n", cl, sum(Idents(pbmcs)==cl), paste(g, collapse=", ")))
}
cat("\n--- author Major.cell.type vs cluster (validation) ---\n")
print(table(Cluster = Idents(pbmcs), Major = pbmcs$Major.cell.type))

saveRDS(pbmcs, file.path(plots_dir, "pbmcs_clustered.rds"))
msg("DONE. elapsed:", round(difftime(Sys.time(), t0, units="mins"),1), "min")
