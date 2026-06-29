#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5 Script 08: CellChat v2 — custom figures
# Bidirectional bubble plots from CellChat interaction tables
# CD8 <-> TAMs only (no epithelial compartment available)
# Two series: MPR vs non-MPR (preprint) + MPR vs pCR (discussion)
# Input:  Results/Tables/BLOC5/GSE243013/Bloc5_07_CellChat_GSE243013_*_interactions.csv
# Output: Results/Figures/BLOC5_Communication/GSE243013/
# Reference: Jin et al., Nature Communications 2021 (CellChat v2)
# ============================================================

library(ggplot2)
library(dplyr)
library(data.table)

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE243013")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5/GSE243013")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load data
# ============================================================
df_MPR    <- fread(file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_MPR_interactions.csv"))
df_nonMPR <- fread(file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_nonMPR_interactions.csv"))
df_pCR    <- fread(file.path(OUT_TAB, "Bloc5_07_CellChat_GSE243013_pCR_interactions.csv"))

df_MPR$condition    <- "MPR"
df_nonMPR$condition <- "non-MPR"
df_pCR$condition    <- "pCR"

# Combine MPR vs non-MPR
df_main <- rbind(df_MPR, df_nonMPR)

# Combine MPR vs pCR
df_pcr <- rbind(df_MPR, df_pCR)

# ============================================================
# 2) Cell type groups
# ============================================================
cd8_types <- c("CD8.TEX", "CD8.TPEX", "CD8.EM", "CD8.CM",
               "CD8.TEMRA", "CD8.NaiveLike", "CD8.MAIT")

tam_types <- c("Proliferating", "Inflammatory Mono-derived", "Resident M2",
               "Lipid-associated", "IFN-stimulated",
               "Stress-response", "Classical Mono-derived")

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
  
  color_map <- c("#4393C3", "#D73027")
  names(color_map) <- c(cond1, cond2)
  
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
# 4) Full figures — MPR vs non-MPR (supplementary)
# ============================================================
message("Generating full figures MPR vs non-MPR...")

p1 <- plot_cellchat(df_main, cd8_types, tam_types,
                    cond1 = "MPR", cond2 = "non-MPR", ntop = 5,
                    title = "CD8 → TAMs — MPR vs non-MPR (GSE243013)")
ggsave(file.path(OUT_FIG, "Bloc5_08_CellChat_CD8_TAMs_MPR_nonMPR.png"),
       p1, width = 16, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs MPR vs non-MPR")

p2 <- plot_cellchat(df_main, tam_types, cd8_types,
                    cond1 = "MPR", cond2 = "non-MPR", ntop = 5,
                    title = "TAMs → CD8 — MPR vs non-MPR (GSE243013)")
ggsave(file.path(OUT_FIG, "Bloc5_08_CellChat_TAMs_CD8_MPR_nonMPR.png"),
       p2, width = 16, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> CD8 MPR vs non-MPR")

# ============================================================
# 5) Full figures — MPR vs pCR (discussion)
# ============================================================
message("Generating full figures MPR vs pCR...")

p3 <- plot_cellchat(df_pcr, cd8_types, tam_types,
                    cond1 = "MPR", cond2 = "pCR", ntop = 5,
                    title = "CD8 → TAMs — MPR vs pCR (GSE243013)")
ggsave(file.path(OUT_FIG, "Bloc5_08_CellChat_CD8_TAMs_MPR_pCR.png"),
       p3, width = 16, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs MPR vs pCR")

p4 <- plot_cellchat(df_pcr, tam_types, cd8_types,
                    cond1 = "MPR", cond2 = "pCR", ntop = 5,
                    title = "TAMs → CD8 — MPR vs pCR (GSE243013)")
ggsave(file.path(OUT_FIG, "Bloc5_08_CellChat_TAMs_CD8_MPR_pCR.png"),
       p4, width = 16, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> CD8 MPR vs pCR")

message("DONE CellChat GSE243013 — full figures")