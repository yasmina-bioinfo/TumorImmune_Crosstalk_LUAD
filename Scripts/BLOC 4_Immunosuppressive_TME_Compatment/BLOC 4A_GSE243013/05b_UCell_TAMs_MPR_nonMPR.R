#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc4A Script 05b: UCell TAMs figures MPR vs non-MPR
# pCR excluded from main preprint narrative
# Input:  Objects/Bloc4A_05_seu_TAMs_UCell.rds
# Output: Results/Figures/BLOC4A_TAMs/Bloc4A_UCell_TAMs_barplot_MPR_nonMPR.png
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(ggpubr)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs/Preprint")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)

# 1) Load object with UCell scores already computed + subset MPR/non-MPR
seu_TAM <- readRDS(file.path(DATA_DIR, "Objects/Bloc4A_05_seu_TAMs_UCell.rds"))
seu_TAM <- subset(seu_TAM, subset = pathological_response %in% c("MPR", "non-MPR"))
message("Cells after pCR exclusion: ", ncol(seu_TAM))

score_cols <- c("M2_immunosuppressive_UCell", "M1_inflammatory_UCell",
                "SPP1_signature_UCell", "IFN_response_UCell")

# 2) Wilcoxon test MPR vs non-MPR
pvals <- sapply(score_cols, function(s) {
  wilcox.test(seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "MPR"],
              seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "non-MPR"])$p.value
})

message("Wilcoxon p-values MPR vs non-MPR:")
print(data.frame(signature = score_cols, p_value = pvals))

pval_df <- data.frame(
  signature = score_cols,
  label = ifelse(round(pvals, 4) == 0, "p < 0.0001",
                 paste0("p = ", round(pvals, 4)))
)

# 3) Barplot
score_summary_bar <- seu_TAM@meta.data %>%
  filter(!is.na(final_annotation), !is.na(pathological_response)) %>%
  select(pathological_response, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols),
               names_to = "signature", values_to = "score") %>%
  group_by(pathological_response, signature) %>%
  summarise(mean_score = mean(score), se = sd(score)/sqrt(n()), .groups = "drop")

p_bar <- ggplot(score_summary_bar,
                aes(x = pathological_response, y = mean_score, fill = pathological_response)) +
  geom_col(width = 0.6, color = "white") +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.2) +
  scale_fill_manual(values = c("MPR" = "#4393C3", "non-MPR" = "#D73027")) +
  facet_wrap(~signature, scales = "free_y", ncol = 2) +
  geom_text(data = pval_df,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 2, size = 3.5) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x    = element_blank(),
        strip.text      = element_text(size = 11, face = "bold")) +
  ylab("Mean UCell score")

ggsave(file.path(OUT_FIG_PREPRINT, "Bloc4A_UCell_TAMs_barplot_MPR_nonMPR.png"),
       p_bar, width = 10, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4A_UCell_TAMs_barplot_MPR_nonMPR.png")

# 4) Violin plots avec p-values — Preprint

meta_vln <- seu_TAM@meta.data %>%
  select(response = pathological_response, all_of(score_cols))

plot_list_vln <- lapply(score_cols, function(score) {
  ggplot(meta_vln, aes(x = response, y = .data[[score]], fill = response)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = c("MPR" = "#4393C3", "non-MPR" = "#D73027")) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       comparisons = list(c("MPR", "non-MPR"))) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title = gsub("_UCell", "", score), y = "UCell score")
})

p_vln <- wrap_plots(plot_list_vln, ncol = 2)

ggsave(file.path(OUT_FIG_PREPRINT, "Bloc4A_UCell_TAMs_violin_MPR_nonMPR_pvalues.png"),
       p_vln, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved: Bloc4A_UCell_TAMs_violin_MPR_nonMPR_pvalues.png")

message("DONE Bloc4A Script 05b")