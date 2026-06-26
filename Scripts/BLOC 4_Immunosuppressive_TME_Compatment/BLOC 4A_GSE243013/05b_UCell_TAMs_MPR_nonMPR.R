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
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")

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

ggsave(file.path(OUT_FIG, "Bloc4A_UCell_TAMs_barplot_MPR_nonMPR.png"),
       p_bar, width = 10, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4A_UCell_TAMs_barplot_MPR_nonMPR.png")
message("DONE Bloc4A Script 05b")