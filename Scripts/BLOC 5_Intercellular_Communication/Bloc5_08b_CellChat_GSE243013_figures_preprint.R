#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5 Script 08b: CellChat v2 — preprint figures
# Filtered on cells of interest guided by Blocs 3 and 4
# CD8 priority: CD8.TEX, CD8.TPEX
# TAMs priority: IFN-stimulated, Stress-response, Resident M2,
#                Lipid-associated, Inflammatory Mono-derived
# Two series: MPR vs non-MPR (preprint) + MPR vs pCR (discussion)
# Input:  Results/Tables/BLOC5/GSE243013/Bloc5_07_CellChat_GSE243013_*_interactions.csv
# Output: Results/Figures/BLOC5_Communication/GSE243013/Preprint/
# Reference: Jin et al., Nature Communications 2021 (CellChat v2)
# ============================================================

library(ggplot2)
library(dplyr)
library(data.table)

DATA_DIR    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG_PRE <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE243013/Preprint")
OUT_TAB     <- file.path(DATA_DIR, "Results/Tables/BLOC5/GSE243013")

dir.create(OUT_FIG_PRE, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load data
# ============================================================
df_MPR    <- fread(file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_MPR_interactions.csv"))
df_nonMPR <- fread(file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_nonMPR_interactions.csv"))
df_pCR    <- fread(file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_pCR_interactions.csv"))

df_MPR$condition    <- "MPR"
df_nonMPR$condition <- "non-MPR"
df_pCR$condition    <- "pCR"

df_main <- rbind(df_MPR, df_nonMPR)
df_pcr  <- rbind(df_MPR, df_pCR)

# ============================================================
# 2) Cells of interest
# ============================================================
cd8_poi <- c("CD8.TEX", "CD8.TPEX")

tam_poi <- c("IFN-stimulated", "Stress-response",
             "Resident M2", "Lipid-associated",
             "Monocyte FCN1+")

# ============================================================
# 3) Plot function
# ============================================================
plot_cellchat <- function(df, source_groups, target_groups,
                          cond1 = "MPR", cond2 = "non-MPR",
                          ntop = 5, title = "") {
  
  df_filt <- df %>%
    filter(source %in% source_groups,
           target %in% target_groups,
           pval < 0.05) %>%
    group_by(target, condition) %>%
    arrange(desc(prob)) %>%
    slice_head(n = ntop) %>%
    ungroup() %>%
    mutate(interaction = interaction_name_2,
           alpha_val = case_when(
             pval == 0   ~ 1.0,
             pval < 0.01 ~ 0.85,
             pval < 0.05 ~ 0.65,
             TRUE        ~ 0.3
           ))
  
  df_filt$condition <- factor(df_filt$condition, levels = c(cond1, cond2))
  
  color_map <- setNames(
    c("MPR" = "#4393C3", "non-MPR" = "#D73027", "pCR" = "#1A7A1A")[c(cond1, cond2)],
    c(cond1, cond2)
  )
  
  ggplot(df_filt, aes(x = interaction, y = target,
                      size = prob, color = condition,
                      alpha = alpha_val)) +
    geom_point(shape = 19) +
    scale_color_manual(values = color_map) +
    scale_size_continuous(range = c(3, 10)) +
    scale_alpha_identity() +
    facet_wrap(~condition, ncol = 1) +
    theme_bw() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
          axis.text.y  = element_text(size = 12, color = "black"),
          strip.text   = element_text(size = 13, face = "bold"),
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 13),
          plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5)) +
    labs(title    = title,
         subtitle = "Opacity: p < 0.001 = opaque | p < 0.01 = semi | p < 0.05 = light",
         x = "Ligand - Receptor", y = "Target",
         size = "Comm. Prob.", color = "Condition")
}

# ============================================================
# 4) Preprint figures — MPR vs non-MPR
# ============================================================
message("Generating preprint figures MPR vs non-MPR...")

p1 <- plot_cellchat(df_main, cd8_poi, tam_poi,
                    cond1 = "MPR", cond2 = "non-MPR", ntop = 5,
                    title = "CD8 → TAMs — MPR vs non-MPR (GSE243013)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_08b_CellChat_CD8_TAMs_MPR_nonMPR_preprint.png"),
       p1, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs MPR vs non-MPR preprint")

p2 <- plot_cellchat(df_main, tam_poi, cd8_poi,
                    cond1 = "MPR", cond2 = "non-MPR", ntop = 5,
                    title = "TAMs → CD8 — MPR vs non-MPR (GSE243013)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_08b_CellChat_TAMs_CD8_MPR_nonMPR_preprint.png"),
       p2, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> CD8 MPR vs non-MPR preprint")

# ============================================================
# 5) Discussion figures — MPR vs pCR
# ============================================================
message("Generating discussion figures MPR vs pCR...")

p3 <- plot_cellchat(df_pcr, cd8_poi, tam_poi,
                    cond1 = "MPR", cond2 = "pCR", ntop = 5,
                    title = "CD8 → TAMs — MPR vs pCR (GSE243013)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_08b_CellChat_CD8_TAMs_MPR_pCR_preprint.png"),
       p3, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs MPR vs pCR")

p4 <- plot_cellchat(df_pcr, tam_poi, cd8_poi,
                    cond1 = "MPR", cond2 = "pCR", ntop = 5,
                    title = "TAMs → CD8 — MPR vs pCR (GSE243013)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_08b_CellChat_TAMs_CD8_MPR_pCR_preprint.png"),
       p4, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> CD8 MPR vs pCR")

message("DONE CellChat GSE243013 — preprint figures")