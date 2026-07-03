#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc3 Script 10b: UCell Hallmark scoring on CD8 T cells
# Unbiased approach — all 50 Hallmark signatures tested
# Discriminant signatures identified by Wilcoxon MPR vs non-MPR
# Focuses on CD8.TEX and CD8.TPEX populations
# Input:  Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
# Output: Results/Figures/CD8/Preprint/Bloc3_UCell_Hallmark_*.png
#         Results/Tables/Bloc3_UCell_Hallmark_scores_summary.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(msigdbr)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(data.table)
  library(patchwork)
  library(BiocParallel)
})

# Paths
DATA_DIR        <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ          <- file.path(DATA_DIR, "Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/CD8/Preprint")
OUT_TAB         <- file.path(DATA_DIR, "Results/Tables")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)

# 1) Load CD8 ProjecTILs object
message("Loading CD8 ProjecTILs object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))
DefaultAssay(seu_CD8) <- "RNA"

# 2) Load ALL Hallmark signatures
message("Loading Hallmark signatures...")
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
message("Total Hallmark signatures: ", length(hallmark_sets))

# 3) Compute UCell scores for all Hallmark signatures
message("Computing UCell Hallmark scores...")
register(SnowParam(workers = 1))
seu_CD8 <- AddModuleScore_UCell(seu_CD8,
                                features = hallmark_sets,
                                name     = "")
message("Hallmark UCell scores added.")

# 4) Subset TEX and TPEX, MPR vs non-MPR only
seu_sub <- subset(seu_CD8,
                  subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX") &
                    pathological_response %in% c("MPR", "non-MPR"))

hallmark_cols <- names(hallmark_sets)

# 5) Wilcoxon test — identify discriminant signatures
message("Running Wilcoxon tests MPR vs non-MPR...")
wilcox_results <- lapply(hallmark_cols, function(sig) {
  mpr  <- seu_sub@meta.data[seu_sub@meta.data$pathological_response == "MPR",  sig]
  nmpr <- seu_sub@meta.data[seu_sub@meta.data$pathological_response == "non-MPR", sig]
  mpr  <- mpr[!is.na(mpr)]
  nmpr <- nmpr[!is.na(nmpr)]
  if (length(mpr) < 3 | length(nmpr) < 3) return(NULL)
  test <- wilcox.test(mpr, nmpr)
  data.frame(signature  = sig,
             p_value    = test$p.value,
             median_MPR = median(mpr),
             median_nonMPR = median(nmpr))
}) %>% bind_rows()

wilcox_results <- wilcox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

fwrite(wilcox_results,
       file.path(OUT_TAB, "Bloc3_GSE243013_UCell_Hallmark_CD8_wilcox_results.csv"))
message("Saved: Bloc3_GSE243013_UCell_Hallmark_CD8_wilcox_results.csv")

# 6) Select discriminant signatures p_adj < 0.05
sig_discriminant <- wilcox_results %>%
  filter(p_adj < 0.05) %>%
  pull(signature)

message("Discriminant Hallmark signatures (p_adj < 0.05): ",
        length(sig_discriminant))
message(paste(sig_discriminant, collapse = ", "))

# 7) Violin plots for discriminant signatures — TEX and TPEX MPR vs non-MPR
if (length(sig_discriminant) > 0) {
  meta_sub <- seu_sub@meta.data %>%
    select(response  = pathological_response,
           cd8_state = functional.cluster,
           all_of(sig_discriminant))
  
  plot_list <- lapply(sig_discriminant, function(sig) {
    ggplot(meta_sub, aes(x = response, y = .data[[sig]], fill = response)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
      scale_fill_manual(values = c("MPR" = "#4393C3", "non-MPR" = "#D73027")) +
      stat_compare_means(method = "wilcox.test", label = "p.format",
                         comparisons = list(c("MPR", "non-MPR"))) +
      facet_wrap(~cd8_state) +
      theme_bw() +
      theme(legend.position = "none",
            axis.title.x    = element_blank(),
            axis.text.x     = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(title = gsub("HALLMARK_", "", sig), y = "UCell score")
  })
  
  # Save in batches of 6 per figure
  n_batch <- ceiling(length(plot_list) / 6)
  for (i in seq_len(n_batch)) {
    idx   <- ((i - 1) * 6 + 1):min(i * 6, length(plot_list))
    p_out <- wrap_plots(plot_list[idx], ncol = 2) +
      plot_annotation(title = "Hallmark UCell scores — CD8 T cells GSE243013",
                      theme = theme(plot.title = element_text(size = 14, face = "bold")))
    ggsave(file.path(OUT_FIG_PREPRINT,
                     paste0("Bloc3_GSE243013_UCell_Hallmark_CD8_discriminant_batch", i, ".png")),
           p_out, width = 16, height = 12, dpi = 300, bg = "white")
    message("Saved batch ", i)
  }
}

# Save object with Hallmark scores
message("Saving CD8 object with Hallmark scores...")
saveRDS(seu_CD8, file.path(DATA_DIR, "Objects/Bloc3_GSE243013_12c_seu_CD8_Hallmark.rds"))
message("Saved: Objects/Bloc3_12c_seu_CD8_Hallmark.rds")

message("DONE Bloc3 Script 10b")