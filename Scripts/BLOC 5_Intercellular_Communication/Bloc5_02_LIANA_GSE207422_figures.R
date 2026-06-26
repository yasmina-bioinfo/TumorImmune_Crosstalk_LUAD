#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 02: LIANA+ figures — main axes
# Produces comparative dotplots for the 3 biological axes:
# CD8 <-> TAMs, CD8 <-> Epithelial, TAMs <-> Epithelial
# Input:  Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv
#         Results/Tables/BLOC5/Bloc5_01_LIANA_GSE207422_NMPR_aggregated.csv
# Output: Results/Figures/BLOC5_Communication/GSE207422/
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
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

# ============================================================
# 1) Load per condition results
# ============================================================
message("Loading LIANA per condition results...")
liana_MPR  <- fread(file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_MPR_aggregated.csv"))
liana_NMPR <- fread(file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_NMPR_aggregated.csv"))

# Define cell type groups
cd8_types <- c("CD8.TEX", "CD8.TPEX", "CD8.EM", "CD8.CM", 
               "CD8.TEMRA", "CD8.NaiveLike")
short_tam_labels <- c(
  "TAM_like_MRC1"                                            = "MRC1+ M2-like",
  "TAM_like_SPP1"                                            = "SPP1+ immunosuppressive",
  "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)" = "Resident M2",
  "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)"                     = "IFN-stimulated",
  "TAM_like_monocyte (classical inflammatory)"               = "Monocyte-derived",
  "TAM_like_lipid (CCL18+/AKR+)"                            = "Lipid-associated",
  "TAM_like_stress (HSP-high/M1-like)"                      = "Stress-response",
  "TAM_like_regulatory (glucocorticoid-responsive)"          = "Regulatory",
  "TAM_like_M2 (SIGLEC8+/CCL18+)"                           = "M2-SIGLEC8+"
)
epi_types  <- c("tumor_MPR", "normal_MPR", "tumor_NMPR", "normal_NMPR", "Ciliated")

# Theme for all figures
theme_liana <- theme(
  axis.text.x  = element_text(angle = 55, hjust = 1, size = 14, color = "black"),
  axis.text.y  = element_text(size = 13, color = "black"),
  axis.title   = element_text(size = 14, color = "black"),
  strip.text   = element_text(size = 13, face = "bold"),
  legend.text  = element_text(size = 13),
  legend.title = element_text(size = 14),
  plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 13, hjust = 0.5)
)

# ============================================================
# 2) Figures 
# ============================================================
message("Generating Figures ...")

# Function to produce MPR/NMPR comparative dotplot
plot_liana_comparison <- function(liana_MPR, liana_NMPR, 
                                  source_groups, target_groups,
                                  ntop = 15, title = "") {
  
  # Filter by source and target groups
  df_MPR <- liana_MPR %>%
    filter(source %in% source_groups, target %in% target_groups,
           aggregate_rank <= 0.05) %>%
    group_by(target) %>%
    arrange(aggregate_rank) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(condition = "MPR",
           interaction = paste0(ligand.complex, " - ", receptor.complex))
  
  df_NMPR <- liana_NMPR %>%
    filter(source %in% source_groups, target %in% target_groups,
           aggregate_rank <= 0.05) %>%
    group_by(target) %>%
    arrange(aggregate_rank) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(condition = "NMPR",
           interaction = paste0(ligand.complex, " - ", receptor.complex))
  
  # Combine
  df_all <- bind_rows(df_MPR, df_NMPR)
  
  # Apply short labels to TAM types
  df_all$source <- ifelse(df_all$source %in% names(short_tam_labels),
                          short_tam_labels[df_all$source], df_all$source)
  df_all$target <- ifelse(df_all$target %in% names(short_tam_labels),
                          short_tam_labels[df_all$target], df_all$target)
  
  df_all$condition <- factor(df_all$condition, levels = c("MPR", "NMPR"))
  
  # Plot
  
  # Simplify epithelial labels
  df_all$target <- case_when(
    df_all$target %in% c("tumor_MPR", "tumor_NMPR") ~ "Tumor epithelial",
    df_all$target %in% c("normal_MPR", "normal_NMPR") ~ "Normal epithelial",
    df_all$target == "Ciliated" ~ "Ciliated",
    TRUE ~ df_all$target
  )
  
  df_all$source <- case_when(
    df_all$source %in% c("tumor_MPR", "tumor_NMPR") ~ "Tumor epithelial",
    df_all$source %in% c("normal_MPR", "normal_NMPR") ~ "Normal epithelial",
    df_all$source == "Ciliated" ~ "Ciliated",
    TRUE ~ df_all$source
  )
  
  ggplot(df_all, aes(x = interaction, y = target,
                     size = sca.LRscore, color = condition)) +
    geom_point() +
    scale_color_manual(values = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
    scale_size_continuous(range = c(3, 10)) +
    facet_wrap(~condition, ncol = 1) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
          axis.text.y  = element_text(size = 13, color = "black"),
          strip.text   = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 13),
          plot.title   = element_text(size = 15, face = "bold", hjust = 0.5)) +
    labs(title = title, x = "Ligand - Receptor", y = "Target",
         size = "LR Score", color = "Condition")
}

# Figure 1 — CD8 -> TAMs
p1 <- plot_liana_comparison(
  liana_MPR, liana_NMPR,
  source_groups = cd8_types,
  target_groups = tam_types,
  ntop = 15,
  title = "CD8 → TAMs interactions — MPR vs NMPR (GSE207422)"
)
ggsave(file.path(OUT_FIG, "Bloc5_02_LIANA_CD8_TAMs.png"),
       p1, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs")

# Figure 2 — CD8 -> Epithelial
p2 <- plot_liana_comparison(
  liana_MPR, liana_NMPR,
  source_groups = cd8_types,
  target_groups = epi_types,
  ntop = 15,
  title = "CD8 → Epithelial interactions — MPR vs NMPR (GSE207422)"
)
ggsave(file.path(OUT_FIG, "Bloc5_02_LIANA_CD8_Epithelial.png"),
       p2, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> Epithelial")

# Figure 3 — TAMs -> Epithelial
p3 <- plot_liana_comparison(
  liana_MPR, liana_NMPR,
  source_groups = tam_types,
  target_groups = epi_types,
  ntop = 15,
  title = "TAMs → Epithelial interactions — MPR vs NMPR (GSE207422)"
)
ggsave(file.path(OUT_FIG, "Bloc5_02_LIANA_TAMs_Epithelial.png"),
       p3, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> Epithelial")

message("DONE LIANA+ GSE207422 figures")