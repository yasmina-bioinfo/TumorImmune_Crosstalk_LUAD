#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 04b: CellChat v2 — preprint figures
# Bidirectional figures on cells of interest guided by Blocs 3 and 4
# CD8 priority: CD8.TEX, CD8.TPEX
# TAMs priority: IFN-stimulated, SPP1+, Lipid-associated, Resident M2
# Epithelial: Tumor + Normal + Ciliated (all kept — Discussion/Perspectives)
# Main preprint narrative: CD8 <-> TAMs only
# Input:  Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_MPR_interactions.csv
#         Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_NMPR_interactions.csv
# Output: Results/Figures/BLOC5_Communication/GSE207422/Preprint/
# Reference: Jin et al., Nature Communications 2021 (CellChat v2)
# ============================================================

library(ggplot2)
library(dplyr)
library(data.table)

DATA_DIR    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG_PRE <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422/Preprint")
dir.create(OUT_FIG_PRE, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load data
# ============================================================
df_MPR  <- fread(file.path(DATA_DIR, "Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_MPR_interactions.csv"))
df_NMPR <- fread(file.path(DATA_DIR, "Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_NMPR_interactions.csv"))

df_MPR$condition  <- "MPR"
df_NMPR$condition <- "NMPR"
df_all <- rbind(df_MPR, df_NMPR)

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

# ============================================================
# 2) Cells of interest
# ============================================================
cd8_poi <- c("CD8.TEX", "CD8.TPEX")

tam_poi <- c("IFN-stimulated", "SPP1+ immunosuppressive",
             "Lipid-associated", "Resident M2")

epi_types <- c("Tumor epithelial", "Normal epithelial", "Ciliated")

# ============================================================
# 3) Plot function
# ============================================================
plot_cellchat <- function(df, source_groups, target_groups, ntop = 5, title = "") {
  
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
  
  df_filt$condition <- factor(df_filt$condition, levels = c("MPR", "NMPR"))
  
  ggplot(df_filt, aes(x = interaction, y = target,
                      size = prob, color = condition,
                      alpha = alpha_val)) +
    geom_point(shape = 19) +
    scale_color_manual(values = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
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
# 4) Preprint figures : CD8 <-> TAMs (main narrative)
# ============================================================
message("Generating preprint figures CD8 <-> TAMs...")

p1 <- plot_cellchat(df_all, cd8_poi, tam_poi, ntop = 5,
                    title = "CD8 → TAMs — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_04b_CellChat_CD8_TAMs_preprint.png"),
       p1, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs preprint")

p2 <- plot_cellchat(df_all, tam_poi, cd8_poi, ntop = 5,
                    title = "TAMs → CD8 — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_04b_CellChat_TAMs_CD8_preprint.png"),
       p2, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: TAMs -> CD8 preprint")

# ============================================================
# 5) Discussion/Perspectives figures : Epithelial axes
# ============================================================
message("Generating discussion figures — Epithelial axes...")

p3 <- plot_cellchat(df_all, cd8_poi, epi_types, ntop = 5,
                    title = "CD8 → Epithelial — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_04b_CellChat_CD8_Epithelial_discussion.png"),
       p3, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: CD8 -> Epithelial discussion")

p4 <- plot_cellchat(df_all, epi_types, cd8_poi, ntop = 5,
                    title = "Epithelial → CD8 — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_04b_CellChat_Epithelial_CD8_discussion.png"),
       p4, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: Epithelial -> CD8 discussion")

p5 <- plot_cellchat(df_all, tam_poi, epi_types, ntop = 5,
                    title = "TAMs → Epithelial — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_04b_CellChat_TAMs_Epithelial_discussion.png"),
       p5, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: TAMs -> Epithelial discussion")

p6 <- plot_cellchat(df_all, epi_types, tam_poi, ntop = 5,
                    title = "Epithelial → TAMs — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG_PRE, "Bloc5_04b_CellChat_Epithelial_TAMs_discussion.png"),
       p6, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: Epithelial -> TAMs discussion")

message("DONE — CellChat GSE207422 preprint + discussion figures")