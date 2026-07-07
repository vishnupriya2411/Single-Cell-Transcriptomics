#!/usr/bin/env Rscript
# cellchat ligand-receptor analysis on the tumor cells.
# builds the object, runs the inference and writes all the figures to Plots/cellchat/.
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(patchwork); library(ggplot2); library(svglite)
})
options(stringsAsFactors = FALSE)
future::plan("sequential")   # keep deterministic / avoid multisession issues in Rscript

plots_dir <- "Plots"
cc_dir <- file.path(plots_dir, "cellchat"); dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)

# savers - write plots to svg (px sizes / res give inches)
cc_save_base <- function(filename, expr, width = 1600, height = 1200, res = 150) {
  filename <- sub("\\.png$", ".svg", filename)
  svglite::svglite(file.path(cc_dir, filename), width = width/res, height = height/res)
  on.exit(dev.off(), add = TRUE)
  tryCatch(force(expr), error = function(e) message("  [skip] ", filename, ": ", conditionMessage(e)))
}
cc_save_obj <- function(filename, plot_obj, width = 1600, height = 1200, res = 150) {
  filename <- sub("\\.png$", ".svg", filename)
  svglite::svglite(file.path(cc_dir, filename), width = width/res, height = height/res)
  on.exit(dev.off(), add = TRUE)
  tryCatch(print(plot_obj), error = function(e) message("  [skip] ", filename, ": ", conditionMessage(e)))
}
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# build the cellchat object (or reuse a saved one to just re-plot)
obj_path <- file.path(cc_dir, "cellchat_object.rds")
if (file.exists(obj_path)) {
  cellchat <- readRDS(obj_path)          # re-plot only: skip the (slow) inference
  cat("Loaded existing cellchat object; skipping inference.\n")
} else {
  pbmcs <- readRDS(file.path(plots_dir, "pbmcs_annotated.rds"))
  pbmcs$CellType <- as.character(Idents(pbmcs))
  tumor <- subset(pbmcs, subset = Tissue == "tumor")
  Idents(tumor) <- tumor$CellType
  cat("tumour cells:", ncol(tumor), " | cell types:", length(unique(tumor$CellType)), "\n")

  data_input <- GetAssayData(tumor, assay = "originalexp", layer = "data")   # log-normalized
  meta <- data.frame(labels = as.character(Idents(tumor)), row.names = colnames(tumor))

  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.mouse
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)  # no presto -> standard Wilcoxon
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
}

group_sizes  <- as.numeric(table(cellchat@idents))
pathways.all <- cellchat@netP$pathways
cat("Significant pathways (", length(pathways.all), "):",
    paste(pathways.all, collapse = ", "), "\n")
write.csv(subsetCommunication(cellchat), file.path(cc_dir, "significant_LR_interactions.csv"), row.names = FALSE)

# aggregated network
cc_save_base("18_aggregate_circle_count_weight.png", {
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count,  vertex.weight = group_sizes,
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = group_sizes,
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
}, width = 2000, height = 1000)

cc_save_base("19_aggregate_circle_per_source.png", {
  mat <- cellchat@net$weight; n <- nrow(mat)
  par(mfrow = c(ceiling(n / 3), 3), xpd = TRUE)
  for (i in seq_len(n)) {
    m <- matrix(0, n, n, dimnames = dimnames(mat)); m[i, ] <- mat[i, ]
    netVisual_circle(m, vertex.weight = group_sizes, weight.scale = TRUE,
                     edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}, width = 1800, height = 600 * ceiling(nrow(cellchat@net$weight) / 3))

cc_save_obj("20_aggregate_heatmap_count.png",
            netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds"))
cc_save_obj("21_aggregate_heatmap_weight.png",
            netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds"))

# global L-R views
cc_save_obj("22_bubble_all_LR.png",
            netVisual_bubble(cellchat, remove.isolate = FALSE), width = 2200, height = 2600)
cc_save_base("23_chord_gene_all.png",
             netVisual_chord_gene(cellchat, lab.cex = 0.6, legend.pos.y = 30),
             width = 1800, height = 1800)

# signaling roles
cc_save_obj("24_signalingRole_heatmap_outgoing.png",
            netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing"), width = 1200, height = 1600)
cc_save_obj("25_signalingRole_heatmap_incoming.png",
            netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming"), width = 1200, height = 1600)
cc_save_obj("26_signalingRole_scatter.png",
            netAnalysis_signalingRole_scatter(cellchat), width = 1200, height = 1000)

# per pathway plots
for (pw in pathways.all) {
  pwn <- safe_name(pw)
  cc_save_base(paste0("pathway_", pwn, "_a_circle.png"),
               netVisual_aggregate(cellchat, signaling = pw, layout = "circle"))
  cc_save_base(paste0("pathway_", pwn, "_b_chord.png"),
               netVisual_aggregate(cellchat, signaling = pw, layout = "chord"))
  cc_save_obj(paste0("pathway_", pwn, "_c_heatmap.png"),
              netVisual_heatmap(cellchat, signaling = pw, color.heatmap = "Reds"))
  cc_save_obj(paste0("pathway_", pwn, "_d_LR_contribution.png"),
              netAnalysis_contribution(cellchat, signaling = pw), width = 1200, height = 800)
  cc_save_base(paste0("pathway_", pwn, "_e_role_network.png"),
               netAnalysis_signalingRole_network(cellchat, signaling = pw,
                                                 width = 8, height = 2.5, font.size = 10),
               width = 1600, height = 700)
  cc_save_obj(paste0("pathway_", pwn, "_f_gene_expression.png"),
              plotGeneExpression(cellchat, signaling = pw), width = 1400, height = 1200)
}
cat("Per-pathway figures written.\n")

# NMF communication patterns (optional)
tryCatch({
  library(NMF); library(ggalluvial)
  cc_save_obj("27_selectK_outgoing.png", selectK(cellchat, pattern = "outgoing"), 1000, 800)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
  cc_save_obj("28_pattern_outgoing_river.png", netAnalysis_river(cellchat, pattern = "outgoing"), 1600, 1000)
  cc_save_obj("29_pattern_outgoing_dot.png",   netAnalysis_dot(cellchat, pattern = "outgoing"), 1400, 900)
  cc_save_obj("30_selectK_incoming.png", selectK(cellchat, pattern = "incoming"), 1000, 800)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 3)
  cc_save_obj("31_pattern_incoming_river.png", netAnalysis_river(cellchat, pattern = "incoming"), 1600, 1000)
  cc_save_obj("32_pattern_incoming_dot.png",   netAnalysis_dot(cellchat, pattern = "incoming"), 1400, 900)
}, error = function(e) message("Communication-pattern step skipped: ", conditionMessage(e)))

saveRDS(cellchat, file = file.path(cc_dir, "cellchat_object.rds"))
cat("DONE run_cellchat — figures in", cc_dir, "\n")
