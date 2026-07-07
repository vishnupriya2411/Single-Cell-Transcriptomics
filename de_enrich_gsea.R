#!/usr/bin/env Rscript
# DE (wilcox) + volcano, then GO/Reactome ORA and GSEA on the result.
# tumor vs healthy across all cell types, ident.1 = tumor so +log2FC is up in tumor.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2)
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
})
has_react <- requireNamespace("ReactomePA", quietly = TRUE)
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
plots_dir <- "Plots"
save_plot <- function(p,f,w=10,h=8,dpi=200){ f <- sub("\\.png$",".svg",f); ggsave(file.path(plots_dir,f),p,width=w,height=h,device=svglite::svglite); cat("saved",f,"\n") }
# shrink the y labels so long pathway names don't overlap
tidy_enr <- function(p, ts = 10) p +
  theme(axis.text.y = element_text(size = ts, lineheight = 0.9),
        strip.text  = element_text(size = 11, face = "bold"),
        plot.title  = element_text(size = 13))

pbmcs <- readRDS(file.path(plots_dir, "pbmcs_annotated.rds"))
Idents(pbmcs) <- pbmcs$Tissue    # compare tumor vs healthy across ALL cell types
cat("Tissue levels:", paste(unique(pbmcs$Tissue), collapse=", "), "\n")

# DE, wilcox, tumor vs healthy over all cells
de <- FindMarkers(pbmcs, ident.1 = "tumor", ident.2 = "healthy",
                  test.use = "wilcox", logfc.threshold = 0, min.pct = 0.1)
de$gene <- rownames(de)
write.csv(de, file.path(plots_dir, "DE_all_tumor_vs_healthy_wilcox.csv"), row.names = FALSE)
cat("DE genes tested:", nrow(de),
    " | padj<0.05 &|LFC|>0.25:", sum(de$p_val_adj<0.05 & abs(de$avg_log2FC)>0.25, na.rm=TRUE), "\n")

# volcano
lfc_cut <- 0.25; p_cut <- 0.05
de$dir <- "ns"
de$dir[de$p_val_adj < p_cut & de$avg_log2FC >  lfc_cut] <- "Up in tumor"
de$dir[de$p_val_adj < p_cut & de$avg_log2FC < -lfc_cut] <- "Down in tumor"
de$neglog10 <- -log10(pmax(de$p_val_adj, .Machine$double.xmin))
lab <- de %>% filter(dir != "ns") %>% arrange(desc(neglog10)) %>% head(20)
v <- ggplot(de, aes(avg_log2FC, neglog10, color = dir)) +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c("Up in tumor"="#c0392b","Down in tumor"="#2c6fbb","ns"="grey75"),
                     name = NULL) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey50") +
  labs(title = "All cell types: tumor (treatment) vs healthy (control) — Wilcoxon",
       x = "avg log2 fold change  (+ = up in tumor)", y = "-log10 adjusted p") +
  theme_bw(base_size = 13) + theme(legend.position = "top")
if (has_repel) v <- v + ggrepel::geom_text_repel(data = lab, aes(label = gene),
                    size = 3, max.overlaps = 30, show.legend = FALSE)
save_plot(v, "17_volcano_all_tumor_vs_healthy.png", 9, 8)

# gene sets for enrichment
sig_genes <- de$gene[de$p_val_adj < p_cut & abs(de$avg_log2FC) > lfc_cut]
cat("significant genes for ORA:", length(sig_genes), "\n")
geneList <- sort(setNames(de$avg_log2FC, de$gene), decreasing = TRUE)   # ranked, for GSEA
geneList <- geneList[!duplicated(names(geneList))]

# GO
# ORA dotplot
ego <- tryCatch(enrichGO(sig_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                         ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                error = function(e){message("enrichGO: ",conditionMessage(e)); NULL})
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  save_plot(tidy_enr(dotplot(ego, showCategory = 15, label_format = 45) +
                     ggtitle("GO BP over-representation (ORA)")),
            "24_GO_ORA_dotplot.png", 10, 10)
  save_plot(tidy_enr(barplot(ego, showCategory = 15, label_format = 45) +
                     ggtitle("GO BP over-representation (ORA)")),
            "16_GO_enrichment_barplot.png", 10, 10)
} else cat("GO ORA: no enriched terms\n")

# GSEA
gse_go <- tryCatch(gseGO(geneList, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
                         pvalueCutoff = 0.05, eps = 0, verbose = FALSE),
                   error = function(e){message("gseGO: ",conditionMessage(e)); NULL})
if (!is.null(gse_go) && nrow(as.data.frame(gse_go)) > 0) {
  save_plot(tidy_enr(dotplot(gse_go, showCategory = 10, split = ".sign", label_format = 40) +
              facet_grid(. ~ .sign) +
              ggtitle("GO BP GSEA (activated vs suppressed in tumor)")),
            "25_GO_GSEA_dotplot.png", 13, 10)
  save_plot(tidy_enr(ridgeplot(gse_go, showCategory = 15, label_format = 45) +
              ggtitle("GO BP GSEA — enrichment distributions")),
            "26_GO_GSEA_ridge.png", 11, 10)
} else cat("GO GSEA: no enriched terms\n")

# Reactome
if (has_react) {
  library(ReactomePA)
  em <- bitr(de$gene, "SYMBOL", "ENTREZID", org.Mm.eg.db)   # symbol -> entrez
  sig_entrez <- bitr(sig_genes, "SYMBOL", "ENTREZID", org.Mm.eg.db)$ENTREZID
  gl_e <- de %>% inner_join(em, by = c("gene"="SYMBOL")) %>%
          group_by(ENTREZID) %>% summarise(l = avg_log2FC[which.max(abs(avg_log2FC))]) %>%
          { setNames(.$l, .$ENTREZID) } %>% sort(decreasing = TRUE)

  er <- tryCatch(enrichPathway(sig_entrez, organism = "mouse", pvalueCutoff = 0.05, readable = TRUE),
                 error = function(e){message("enrichPathway: ",conditionMessage(e)); NULL})
  if (!is.null(er) && nrow(as.data.frame(er)) > 0)
    save_plot(tidy_enr(dotplot(er, showCategory = 15, label_format = 50) +
                       ggtitle("Reactome pathway ORA")),
              "27_reactome_ORA_dotplot.png", 11, 10)
  else cat("Reactome ORA: no enriched pathways\n")

  gse_r <- tryCatch(gsePathway(gl_e, organism = "mouse", pvalueCutoff = 0.05, eps = 0, verbose = FALSE),
                    error = function(e){message("gsePathway: ",conditionMessage(e)); NULL})
  if (!is.null(gse_r) && nrow(as.data.frame(gse_r)) > 0) {
    gse_r <- setReadable(gse_r, org.Mm.eg.db, keyType = "ENTREZID")
    save_plot(tidy_enr(dotplot(gse_r, showCategory = 10, split = ".sign", label_format = 45) +
                facet_grid(. ~ .sign) +
                ggtitle("Reactome GSEA (activated vs suppressed in tumor)")),
              "28_reactome_GSEA_dotplot.png", 13, 10)
  } else cat("Reactome GSEA: no enriched pathways\n")
} else cat("ReactomePA not installed — skipping Reactome\n")

cat("DONE de_enrich_gsea\n")
