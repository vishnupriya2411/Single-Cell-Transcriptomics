#!/usr/bin/env Rscript
# apply the cluster labels to the saved object and make the annotated umap +
# canonical marker dotplot to check them.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })
plots_dir <- "Plots"
save_plot <- function(p,f,w=10,h=7,dpi=200) ggsave(file.path(plots_dir,f),p,width=w,height=h,dpi=dpi)

pbmcs <- readRDS(file.path(plots_dir, "pbmcs_clustered.rds"))

# cluster -> label, from the top markers and the author Major.cell.type
labels <- c(
  "0"="Neutrophils (mature)", "1"="B cells", "2"="T cells",
  "3"="Neutrophils (Siglecf+ TAN)", "4"="Monocytes", "5"="NK cells",
  "6"="Neutrophils (immature)", "7"="Neutrophils (Ccl3+)",
  "8"="Alveolar macrophages", "9"="cDC2", "10"="mregDC",
  "11"="Macrophages (C1q+Trem2+ TAM)")
pbmcs <- RenameIdents(pbmcs, labels)
pbmcs$CellType_manual <- Idents(pbmcs)

# order identities lymphoid -> myeloid for tidy plots
lvl <- c("B cells","T cells","NK cells","Monocytes","Alveolar macrophages",
         "Macrophages (C1q+Trem2+ TAM)","cDC2","mregDC","Neutrophils (mature)",
         "Neutrophils (immature)","Neutrophils (Ccl3+)","Neutrophils (Siglecf+ TAN)")
Idents(pbmcs) <- factor(Idents(pbmcs), levels = lvl)

# Annotated UMAP
save_plot(DimPlot(pbmcs, label = TRUE, repel = TRUE, label.size = 3.2) +
            ggtitle("Manually annotated cell types (top-5 marker based)") +
            theme(legend.position = "right"),
          "22_UMAP_annotated_manual.png", 11, 7)

# dotplot of canonical markers as evidence for the labels
markers <- c("Cd79a","Ms4a1","Cd19",              # B
             "Cd3g","Lat","Thy1",                 # T
             "Ncr1","Gzma","Klra4",               # NK
             "Ace","F13a1","Clec4a1",             # Monocytes
             "Krt79","Cidec",                     # Alveolar mac
             "C1qa","C1qb","Trem2","Mmp12",       # TAM / interstitial mac
             "Cd209a","Clec10a","Ccl17",          # cDC2
             "Fscn1","Ccl22",                     # mregDC
             "S100a8","Ly6g","Retnlg","Ngp","Mmp8","Ccl3","Siglecf","Ffar2") # neutrophils
markers <- markers[markers %in% rownames(pbmcs)]
dp <- DotPlot(pbmcs, features = markers) + RotatedAxis() +
  ggtitle("Canonical marker expression per annotated cell type") +
  theme(axis.text.x = element_text(size = 8))
save_plot(dp, "23_canonical_marker_dotplot.png", 15, 7)

# annotation summary table
summ <- as.data.frame(table(CellType = Idents(pbmcs)))
write.csv(summ, file.path(plots_dir, "annotation_summary.csv"), row.names = FALSE)
cat("Annotated cell counts:\n"); print(summ)
saveRDS(pbmcs, file.path(plots_dir, "pbmcs_annotated.rds"))
cat("Phase 2 done.\n")
