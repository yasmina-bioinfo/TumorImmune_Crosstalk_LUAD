#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5 Script 06: LIANA+ figures — main axis
# Produces comparative dotplots for CD8 <-> TAMs axis
# Three figure series:
#   - Preprint: MPR vs non-MPR
#   - Discussion: MPR vs pCR
# Input:  Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_MPR_aggregated.csv
#         Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_non_MPR_aggregated.csv
#         Results/Tables/BLOC5/GSE243013/Bloc5_03_LIANA_GSE243013_pCR_aggregated.csv
# Output: Results/Figures/BLOC5_Communication/GSE243013/
# Reference: Dimitrov et al., Nature Communications 2022 (LIANA+)
# ============================================================

suppressPackageStartupMessages({
  library(liana)
  library(tidyverse)
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE243013")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5/GSE243013")

# ============================================================
# 1) Load per condition results
# ============================================================
message("Loading LIANA per condition results...")
liana_MPR    <- fread(file.path(OUT_TAB, "Bloc5_03_LIANA_GSE243013_MPR_aggregated.csv"))
liana_nonMPR <- fread(file.path(OUT_TAB, "Bloc5_03_LIANA_GSE243013_non_MPR_aggregated.csv"))
liana_pCR    <- fread(file.path(OUT_TAB, "Bloc5_03_LIANA_GSE243013_pCR_aggregated.csv"))

# Define cell type groups
cd8_types <- c("CD8.TEX", "CD8.TPEX", "CD8.EM", "CD8.CM",
               "CD8.TEMRA", "CD8.NaiveLike")

short_tam_labels <- c(
  "Proliferating TAMs (cycling/MKI67+)"                                = "Proliferating",
  "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)"                 = "Inflammatory Mono-derived",
  "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)"  = "Resident M2",
  "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)"        = "Lipid-associated",
  "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)"      = "IFN-stimulated",
  "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)"    = "Stress-response",
  "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"           = "Classical Mono-derived"
)

# ============================================================
# 2) Comparison function
# ============================================================
plot_liana_comparison <- function(liana_cond1, liana_cond2,
                                  source_groups, target_groups,
                                  label_cond1 = "MPR", label_cond2 = "non-MPR",
                                  ntop = 5, title = "") {
  
  process_cond <- function(df, label) {
    df %>%
      filter(source %in% source_groups, target %in% target_groups,
             aggregate_rank <= 0.05) %>%
      group_by(target) %>%
      arrange(aggregate_rank) %>%
      slice_head(n = ntop) %>%
      ungroup() %>%
      mutate(condition = label,
             interaction = paste0(ligand.complex, " - ", receptor.complex))
  }
  
  df_all <- bind_rows(
    process_cond(liana_cond1, label_cond1),
    process_cond(liana_cond2, label_cond2)
  )
  
  # Apply short TAM labels
  df_all$source <- ifelse(df_all$source %in% names(short_tam_labels),
                          short_tam_labels[df_all$source], df_all$source)
  df_all$target <- ifelse(df_all$target %in% names(short_tam_labels),
                          short_tam_labels[df_all$target], df_all$target)
  
  df_all$condition <- factor(df_all$condition, levels = c(label_cond1, label_cond2))
  
  # Color mapping
  color_map <- c("#4393C3", "#D73027")
  names(color_map) <- c(label_cond1, label_cond2)
  
  ggplot(df_all, aes(x = interaction, y = target,
                     size = sca.LRscore, color = condition)) +
    geom_point() +
    scale_color_manual(values = color_map) +
    scale_size_continuous(range = c(3, 10)) +
    facet_wrap(~condition, ncol = 1) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
          axis.text.y  = element_text(size = 14, color = "black"),
          strip.text   = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 13),
          plot.title   = element_text(size = 15, face = "bold", hjust = 0.5)) +
    labs(title = title, x = "Ligand - Receptor", y = "Target",
         size = "LR Score", color = "Condition")
}

# ============================================================
# 3) Figures preprint — MPR vs non-MPR
# ============================================================
message("Generating preprint figures: MPR vs non-MPR...")

p1 <- plot_liana_comparison(
  liana_MPR, liana_nonMPR,
  source_groups = cd8_types,
  target_groups = names(short_tam_labels),
  label_cond1 = "MPR", label_cond2 = "non-MPR",
  title = "CD8 → TAMs interactions — MPR vs non-MPR (GSE243013)"
)
ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_CD8_TAMs_MPR_nonMPR.png"),
       p1, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs MPR vs non-MPR")

p2 <- plot_liana_comparison(
  liana_nonMPR, liana_MPR,
  source_groups = names(short_tam_labels),
  target_groups = cd8_types,
  label_cond1 = "non-MPR", label_cond2 = "MPR",
  title = "TAMs → CD8 interactions — MPR vs non-MPR (GSE243013)"
)
ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_TAMs_CD8_MPR_nonMPR.png"),
       p2, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> CD8 MPR vs non-MPR")

# ============================================================
# 4) Figures discussion — MPR vs pCR
# ============================================================
message("Generating discussion figures: MPR vs pCR...")

p3 <- plot_liana_comparison(
  liana_MPR, liana_pCR,
  source_groups = cd8_types,
  target_groups = names(short_tam_labels),
  label_cond1 = "MPR", label_cond2 = "pCR",
  title = "CD8 → TAMs interactions — MPR vs pCR (GSE243013)"
)
ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_CD8_TAMs_MPR_pCR.png"),
       p3, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs MPR vs pCR")

p4 <- plot_liana_comparison(
  liana_MPR, liana_pCR,
  source_groups = names(short_tam_labels),
  target_groups = cd8_types,
  label_cond1 = "MPR", label_cond2 = "pCR",
  title = "TAMs → CD8 interactions — MPR vs pCR (GSE243013)"
)
ggsave(file.path(OUT_FIG, "Bloc5_04_LIANA_TAMs_CD8_MPR_pCR.png"),
       p4, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> CD8 MPR vs pCR")

message("DONE LIANA+ GSE243013 figures")