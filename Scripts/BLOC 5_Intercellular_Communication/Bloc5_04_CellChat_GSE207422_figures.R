#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 04: CellChat v2 — custom figures
# Produces readable bubble plots from CellChat interaction tables
# MPR vs NMPR — 3 biological axes: CD8-TAMs, CD8-Epithelial, TAMs-Epithelial
# Input:  Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_MPR_interactions.csv
#         Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_NMPR_interactions.csv
# Output: Results/Figures/BLOC5_Communication/GSE207422/
# Reference: Jin et al., Nature Communications 2021 (CellChat v2)
# ============================================================

library(ggplot2)
library(dplyr)
library(data.table)

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422")

# Load data
df_MPR  <- fread(file.path(DATA_DIR, "Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_MPR_interactions.csv"))
df_NMPR <- fread(file.path(DATA_DIR, "Results/Tables/BLOC5/Bloc5_03_CellChat_GSE207422_NMPR_interactions.csv"))

# Add condition column
df_MPR$condition  <- "MPR"
df_NMPR$condition <- "NMPR"

# Combine
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

# Define groups
cd8_types <- c("CD8.TEX", "CD8.TPEX", "CD8.EM", "CD8.CM", "CD8.TEMRA")
tam_types <- c("MRC1+ M2-like", "SPP1+ immunosuppressive", "IFN-stimulated",
               "M2-SIGLEC8+", "Resident M2", "Monocyte-derived", "Lipid-associated",
               "Stress-response", "Regulatory")
epi_types  <- c("Tumor epithelial", "Normal epithelial", "Ciliated")

# Function to plot
plot_cellchat <- function(df, source_groups, target_groups, ntop = 15, title = "") {
  
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
  
  p <- ggplot(df_filt, aes(x = interaction, y = target,
                           size = prob, color = condition,
                           alpha = alpha_val)) +
    geom_point(shape = 19) +
    scale_color_manual(values = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
    scale_size_continuous(range = c(3, 10)) +
    scale_alpha_identity() +
    facet_wrap(~condition, ncol = 1) +
    theme_bw() +
    theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
          axis.text.y   = element_text(size = 12, color = "black"),
          strip.text    = element_text(size = 13, face = "bold"),
          legend.text   = element_text(size = 12),
          legend.title  = element_text(size = 13),
          plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5)) +
    labs(title    = title,
         subtitle = "Opacity: p < 0.001 = opaque | p < 0.01 = semi | p < 0.05 = light",
         x = "Ligand - Receptor", y = "Target",
         size = "Comm. Prob.", color = "Condition")
  
  return(p)
}

# Figure 1 — CD8 -> TAMs
p1 <- plot_cellchat(df_all, cd8_types, tam_types, ntop = 5,
                    title = "CD8 → TAMs — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG, "Bloc5_04_CellChat_CD8_TAMs.png"),
       p1, width = 18, height = 10, dpi = 300, bg = "white")
message("Saved: CD8 -> TAMs")

# Figure 2 — CD8 -> Epithelial
p2 <- plot_cellchat(df_all, cd8_types, epi_types, ntop = 5,
                    title = "CD8 → Epithelial — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG, "Bloc5_04_CellChat_CD8_Epithelial.png"),
       p2, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: CD8 -> Epithelial")

# Figure 3 — TAMs -> Epithelial
p3 <- plot_cellchat(df_all, tam_types, epi_types, ntop = 5,
                    title = "TAMs → Epithelial — MPR vs NMPR (GSE207422)")
ggsave(file.path(OUT_FIG, "Bloc5_04_CellChat_TAMs_Epithelial.png"),
       p3, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: TAMs -> Epithelial")