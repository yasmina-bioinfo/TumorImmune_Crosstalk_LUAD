#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Bloc4B Script 05: UCell scoring on TAMs
# Computes functional gene module scores per TAM subtype
# Compares scores across MPR and NMPR
# Input:  Objects/Bloc4B_04_seu_TAMs_combined.rds (TumorImmune repo)
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_UCell_TAMs_combined_*.png
#         Results/Tables/BLOC4B/Bloc4B_UCell_TAMs_combined_scores_summary.csv
# Reference: Andreatta & Carmona 2021 (UCell)
#            Chen et al. 2021 (M1/M2 signatures — PMC8053174)
#            Italiani & Boraschi 2019 (M1/M2 markers — PMC6543837)
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(patchwork)
  library(ggpubr)
  library(tidyr)
})

# Paths
DATA_DIR_PORTFOLIO <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy"
DATA_DIR_OUTPUT    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ  <- file.path(DATA_DIR_OUTPUT, "Objects/Bloc4B_04_seu_TAMs_combined.rds")
OUT_FIG <- file.path(DATA_DIR_OUTPUT, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB <- file.path(DATA_DIR_OUTPUT, "Results/Tables/BLOC4B")

# 1) Load TAMs object
message("Loading TAMs object...")
seu_TAM <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_TAM))

# 2) Define gene signatures
# NOTE: M1/M2 signatures from Chen et al. 2021 (PMC8053174) and
#       Italiani & Boraschi 2019 (PMC6543837)
#       SPP1 and IFN signatures from Cheng et al. 2021 (Cell)
signatures <- list(
  M2_immunosuppressive = c("MRC1", "CD163", "TGFB1", "IL10", "VEGFA",
                           "CD274", "IDO1", "CSF1R"),
  M1_inflammatory      = c("TNF", "IL1B", "IL6", "CXCL10", "NOS2"),
  SPP1_signature       = c("SPP1", "GPNMB", "APOE", "TREM2"),
  IFN_response         = c("ISG15", "IFIT1", "IFIT3", "CXCL9", "CXCL10")
)

# 3) Compute UCell scores
message("Computing UCell scores...")
DefaultAssay(seu_TAM) <- "RNA"
seu_TAM <- JoinLayers(seu_TAM)
seu_TAM <- AddModuleScore_UCell(seu_TAM,
                                features = signatures,
                                name     = "_UCell")

message("UCell scores added:")
print(head(seu_TAM@meta.data[, grep("_UCell", colnames(seu_TAM@meta.data))]))

# 4) Violin plots: scores by pathological response (MPR/NMPR)
message("Generating violin plots by response...")

score_cols <- paste0(names(signatures), "_UCell")

plot_list <- lapply(score_cols, function(score) {
  VlnPlot(seu_TAM,
          features = score,
          group.by = "PathResponse",
          pt.size  = 0,
          cols     = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank())
})

p_vln_response <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4B_UCell_TAMs_combined_violin_response.png"),
       p_vln_response, width = 12, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4B_UCell_TAMs_combined_violin_response.png")

# 5) Violin plots : scores by TAM subtype (9 in total)
message("Generating violin plots by TAM subtype...")

plot_list2 <- lapply(score_cols, function(score) {
  VlnPlot(seu_TAM,
          features = score,
          group.by = "combined_annotation",
          pt.size  = 0,
          cols = c(
            "TAM_like_MRC1" = "#A50026",
            "TAM_like_SPP1" = "#1A7A1A",
            "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)" = "#D73027",
            "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)" = "#4575B4",
            "TAM_like_monocyte (classical inflammatory)" = "#FDAE61",
            "TAM_like_lipid (CCL18+/AKR+)" = "#FEE090",
            "TAM_like_stress (HSP-high/M1-like)" = "#F46D43",
            "TAM_like_regulatory (glucocorticoid-responsive)" = "#ABD9E9",
            "TAM_like_M2 (SIGLEC8+/CCL18+)" = "#74ADD1"
          )) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1))
})

p_vln_state <- wrap_plots(plot_list2, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4B_UCell_TAMs_combined_violin_subtype.png"),
       p_vln_state, width = 12, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4B_UCell_TAMs_combined_violin_subtype.png")

# 6) Summary table : mean scores per response group per TAM subtype
message("Generating summary table...")

summary_scores <- seu_TAM@meta.data %>%
  filter(!is.na(combined_annotation), !is.na(PathResponse)) %>%
  group_by(PathResponse, combined_annotation) %>%
  summarise(across(all_of(score_cols), mean, .names = "mean_{.col}"),
            n_cells = n(),
            .groups = "drop")

print(summary_scores)
fwrite(summary_scores, file.path(OUT_TAB, "Bloc4B_UCell_TAMs_combined_scores_summary.csv"))
message("Saved: Bloc4B_UCell_TAMs_combined_scores_summary.csv")

# 7) Barplot : mean UCell scores by response group with statistical test
score_summary_bar <- seu_TAM@meta.data %>%
  filter(!is.na(combined_annotation), !is.na(PathResponse)) %>%
  select(PathResponse, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols),
               names_to = "signature",
               values_to = "score") %>%
  group_by(PathResponse, signature) %>%
  summarise(mean_score = mean(score),
            se = sd(score)/sqrt(n()),
            .groups = "drop")
# HOW TO COMPUTE A WILCOXON P-VALUE MANUALLY IN R:
# wilcox.test(group1_values, group2_values)$p.value
# Example for M2 signature:
# wilcox.test(
#   seu_TAM@meta.data$M2_immunosuppressive_UCell[seu_TAM@meta.data$PathResponse == "MPR"],
#   seu_TAM@meta.data$M2_immunosuppressive_UCell[seu_TAM@meta.data$PathResponse == "NMPR"]
# )$p.value
# sapply() applies this to all signatures at once and returns a named vector

# Calcul of p-values : Wilcoxon test cell-level
pvals <- sapply(score_cols, function(s) {
  wilcox.test(seu_TAM@meta.data[[s]][seu_TAM@meta.data$PathResponse == "MPR"],
              seu_TAM@meta.data[[s]][seu_TAM@meta.data$PathResponse == "NMPR"])$p.value
})
print(round(pvals, 4))

# Manual p-value annotations dataframe
# NOTE: round(..., 4) returns 0 for very small p-values (e.g. p < 0.0001)
# which is misleading. ifelse() checks if the rounded value equals 0 and
# replaces it with "p < 0.0001" for clarity in the figure.
# The original paste0("p = ", round(..., 4)) is kept as comment above
# to show the base logic before this correction was applied.
pval_df <- data.frame(
  signature = score_cols,
  label = c(
   # paste0("p = ", round(pvals["M2_immunosuppressive_UCell"], 4)),
    #paste0("p = ", round(pvals["M1_inflammatory_UCell"], 4)),
    #paste0("p = ", round(pvals["SPP1_signature_UCell"], 4)),
    #paste0("p = ", round(pvals["IFN_response_UCell"], 4))
    
    ifelse(round(pvals["M2_immunosuppressive_UCell"], 4) == 0, "p < 0.0001", paste0("p = ", round(pvals["M2_immunosuppressive_UCell"], 4))),
    paste0("p = ", round(pvals["M1_inflammatory_UCell"], 4)),
    ifelse(round(pvals["SPP1_signature_UCell"], 4) == 0, "p < 0.0001", paste0("p = ", round(pvals["SPP1_signature_UCell"], 4))),
    paste0("p = ", round(pvals["IFN_response_UCell"], 4))
  )
)

p_bar_ucell <- ggplot(score_summary_bar,
                      aes(x = PathResponse, y = mean_score, fill = PathResponse)) +
  geom_col(width = 0.6, color = "white") +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                width = 0.2) +
  scale_fill_manual(values = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
  facet_wrap(~signature, scales = "free_y", ncol = 2) +
  geom_text(data = pval_df,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 2, size = 3.5) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x    = element_blank(),
        strip.text      = element_text(size = 11, face = "bold")) +
  ylab("Mean UCell score")

ggsave(file.path(OUT_FIG, "Bloc4B_UCell_TAMs_combined_barplot_scores.png"),
       p_bar_ucell, width = 10, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4B_UCell_TAMs_combined_barplot_scores.png")

# 8) Save updated object
message("Saving updated TAMs object with UCell scores...")
saveRDS(seu_TAM, file.path(DATA_DIR_OUTPUT, "Objects/Bloc4B_05_seu_TAMs_combined_UCell.rds"))
message("Saved: Objects/Bloc4B_05_seu_TAMs_combined_UCell.rds")

# 9) Dotplot by TAM subtype and response
p_dot <- DotPlot(seu_TAM,
                 features = score_cols,
                 group.by = "combined_annotation",
                 split.by = "PathResponse",
                 cols = c("MPR" = "#E63946", "NMPR" = "#457B9D"),
                 dot.scale = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title  = element_text(size = 13, color = "black", face = "bold"),
        panel.grid.major = element_line(color = "grey90")) +
  ggtitle("UCell scores by TAM subtype and response — GSE207422") +
  labs(subtitle = "Red = MPR | Blue = NMPR") +
  xlab("") + ylab("")

ggsave(file.path(OUT_FIG, "Bloc4B_UCell_Dotplot_TAMs_split.png"),
       p_dot, width = 14, height = 10, dpi = 300, bg = "white")
message("Saved: Bloc4B_UCell_Dotplot_TAMs_split.png")

message("DONE Bloc4B Script 05")