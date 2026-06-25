#!/usr/bin/env Rscript
# ============================================================
# GSE207422  Bloc4A Script 05: UCell scoring on TAMs
# Computes functional gene module scores per TAM subtype
# Compares scores across MPR and NMPR
# Input:  Objects/Bloc4A_04_seu_TAMs.rds (TumorImmune repo)
# Output: Results/Figures/Bloc4A_Epithelial_TAMs/Bloc4A_UCell_TAMs_*.png
#         Results/Tables/Bloc4A/Bloc4A_UCell_TAMs_scores_summary.csv
# Reference: Andreatta & Carmona 2021 (UCell)
#            Chen et al. 2021 (M1/M2 signatures  PMC8053174)
#            Italiani & Boraschi 2019 (M1/M2 markers  PMC6543837)
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
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4A_04_seu_TAMs_annotated.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4A_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC4A")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")

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
          group.by = "pathological_response",
          pt.size  = 0,
          cols     = c("MPR" = "#4393C3", "non-MPR" = "#D73027", "pCR" = "#1A7A1A")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank())
})

p_vln_response <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4A_UCell_TAMs_violin_response.png"),
       p_vln_response, width = 12, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4A_UCell_TAMs_violin_response.png")

# 5) Violin plots : scores by TAM subtypes (7 in total)
message("Generating violin plots by TAM subtype...")

# Short labels for readability
short_labels <- c(
  "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)" = "Resident M2",
  "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)"        = "LAMs",
  "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)"                 = "Monocyte FCN1+",
  "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)"    = "Stress-response",
  "Proliferating TAMs (cycling/MKI67+)"                                = "Proliferating",
  "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)"      = "IFN-stimulated",
  "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"           = "Classical monocyte"
)

plot_list2 <- lapply(score_cols, function(score) {
  VlnPlot(seu_TAM,
          features = score,
          group.by = "final_annotation",
          pt.size  = 0,
          cols     = c(
            "Tissue-resident immunosuppressive TAMs (anti-inflammatory/M2-like)" = "#D73027",
            "TREM2+/APOE+ lipid-associated immunosuppressive TAMs (LAMs)"        = "#A50026",
            "Inflammatory monocyte-derived TAMs (FCN1+/S100A8+)"                 = "#FDAE61",
            "Stress-response immunosuppressive TAMs (MARCO+/PPARG+/HSP-high)"    = "#F46D43",
            "Proliferating TAMs (cycling/MKI67+)"                                = "#ABD9E9",
            "IFN-stimulated immunomodulatory TAMs (ISG-high/PD-L1+/IDO1+)"      = "#4575B4",
            "Classical monocyte-derived TAMs (S100A8+/S100A9+/CCR2+)"           = "#D6604D"
          )) +
    scale_x_discrete(labels = short_labels) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1, size = 9)) +
    labs(title = score)
})

p_vln_state <- wrap_plots(plot_list2, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4A_UCell_TAMs_violin_subtype.png"),
       p_vln_state, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4A_UCell_TAMs_violin_subtype.png")

# 6) Summary table : mean scores per response group per TAM subtype
message("Generating summary table...")

summary_scores <- seu_TAM@meta.data %>%
  filter(!is.na(final_annotation), !is.na(pathological_response)) %>%
  group_by(pathological_response, final_annotation) %>%
  summarise(across(all_of(score_cols), mean, .names = "mean_{.col}"),
            n_cells = n(),
            .groups = "drop")

print(summary_scores)
fwrite(summary_scores, file.path(OUT_TAB, "Bloc4A_UCell_TAMs_scores_summary.csv"))
message("Saved: Bloc4A_UCell_TAMs_scores_summary.csv")

# 7) Barplot : mean UCell scores by response group with statistical test
score_summary_bar <- seu_TAM@meta.data %>%
  filter(!is.na(final_annotation), !is.na(pathological_response)) %>%
  select(pathological_response, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols),
               names_to = "signature",
               values_to = "score") %>%
  group_by(pathological_response, signature) %>%
  summarise(mean_score = mean(score),
            se = sd(score)/sqrt(n()),
            .groups = "drop")

# Pairwise Wilcoxon tests with Bonferroni correction (3 comparisons)
# Reference: Wilcoxon 1945; Jaakkola et al. 2017 (PMC6979262)
# Bonferroni threshold: p < 0.05/3 = 0.0167
pvals_MPR_nonMPR <- sapply(score_cols, function(s) {
  wilcox.test(seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "MPR"],
              seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "non-MPR"])$p.value
})
pvals_MPR_pCR <- sapply(score_cols, function(s) {
  wilcox.test(seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "MPR"],
              seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "pCR"])$p.value
})
pvals_nonMPR_pCR <- sapply(score_cols, function(s) {
  wilcox.test(seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "non-MPR"],
              seu_TAM@meta.data[[s]][seu_TAM@meta.data$pathological_response == "pCR"])$p.value
})
print("MPR vs non-MPR:"); print(round(pvals_MPR_nonMPR, 4))
print("MPR vs pCR:"); print(round(pvals_MPR_pCR, 4))
print("non-MPR vs pCR:"); print(round(pvals_nonMPR_pCR, 4))

# Manual p-value annotations dataframe : 3 pairwise comparisons per signature
# NOTE: ifelse() replaces rounded 0 with "p < 0.0001" for clarity
# Format: "MPR/non-MPR | MPR/pCR | non-MPR/pCR"
format_p <- function(p) {
  ifelse(round(p, 4) == 0, "< 0.0001", as.character(round(p, 4)))
}

pval_df <- data.frame(
  signature = score_cols,
  label = paste0(
    "MPR/non-MPR: p", format_p(pvals_MPR_nonMPR), "\n",
    "MPR/pCR: p", format_p(pvals_MPR_pCR), "\n",
    "non-MPR/pCR: p", format_p(pvals_nonMPR_pCR)
  )
)

p_bar_ucell <- ggplot(score_summary_bar,
                      aes(x = pathological_response, y = mean_score, fill = pathological_response)) +
  geom_col(width = 0.6, color = "white") +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                width = 0.2) +
  scale_fill_manual(values = c("MPR" = "#4393C3", "non-MPR" = "#D73027", "pCR" = "#1A7A1A"))+
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

ggsave(file.path(OUT_FIG, "Bloc4A_UCell_TAMs_barplot_scores.png"),
       p_bar_ucell, width = 10, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4A_UCell_TAMs_barplot_scores.png")


# 8) Save updated object
message("Saving updated TAMs object with UCell scores...")
saveRDS(seu_TAM, file.path(DATA_DIR_OUTPUT, "Objects/Bloc4A_05_seu_TAMs_UCell.rds"))
message("Saved: Objects/Bloc4A_05_seu_TAMs_UCell.rds")

# 9) Dotplot by TAM subtype and response
p_dot <- DotPlot(seu_TAM,
                 features = score_cols,
                 group.by = "final_annotation",
                 split.by = "pathological_response",
                 cols = c("MPR" = "#E63946", "non-MPR" = "#457B9D", "pCR" = "#2A9D8F"),
                 dot.scale = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        plot.title  = element_text(size = 13, color = "black", face = "bold"),
        panel.grid.major = element_line(color = "grey90")) +
  ggtitle("UCell scores by TAM subtype and response — GSE243013") +
  labs(subtitle = "Red = MPR | Blue = non-MPR | Green = pCR") +
  xlab("") + ylab("")

ggsave(file.path(OUT_FIG, "Bloc4A_UCell_Dotplot_TAMs_split.png"),
       p_dot, width = 14, height = 10, dpi = 300, bg = "white")

message("Saved: Bloc4A_UCell_Dotplot_TAMs_split.png")
message("DONE Bloc4A Script 05")