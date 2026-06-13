#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Bloc4B Script 02: UCell scoring on TAMs
# Computes functional gene module scores per TAM subtype
# Compares scores across MPR and NMPR
# Input:  Objects/Bloc4B_01_seu_TAMs_annotated.rds (TumorImmune repo)
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_UCell_TAMs_*.png
#         Results/Tables/BLOC4B/Bloc4B_UCell_TAMs_scores_summary.csv
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
})

# Paths
DATA_DIR_PORTFOLIO <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy"
DATA_DIR_OUTPUT    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ  <- file.path(DATA_DIR_OUTPUT, "Objects/Bloc4B_01_seu_TAMs_annotated.rds")
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

# 4) Violin plots — scores by pathological response
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
ggsave(file.path(OUT_FIG, "Bloc4B_UCell_TAMs_violin_response.png"),
       p_vln_response, width = 12, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4B_UCell_TAMs_violin_response.png")

# 5) Violin plots — scores by TAM subtype
message("Generating violin plots by TAM subtype...")

plot_list2 <- lapply(score_cols, function(score) {
  VlnPlot(seu_TAM,
          features = score,
          group.by = "final_annotation",
          pt.size  = 0,
          cols     = c("TAM_like"      = "#D73027",
                       "TAM_like_MRC1" = "#A50026",
                       "TAM_like_SPP1" = "#4575B4")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1))
})

p_vln_state <- wrap_plots(plot_list2, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4B_UCell_TAMs_violin_subtype.png"),
       p_vln_state, width = 12, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4B_UCell_TAMs_violin_subtype.png")

# 6) Focus: TAM_like_MRC1 and TAM_like_SPP1 — compare scores across response
# NOTE: TAM_like_SPP1 small (n=362) — interpret with caution
message("Focus analysis: TAM_like_MRC1 and TAM_like_SPP1...")

seu_focus <- subset(seu_TAM,
                    subset = final_annotation %in% c("TAM_like_MRC1", "TAM_like_SPP1"))

plot_list3 <- lapply(score_cols, function(score) {
  VlnPlot(seu_focus,
          features = score,
          group.by = "PathResponse",
          split.by = "final_annotation",
          pt.size  = 0) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1))
})

p_vln_focus <- wrap_plots(plot_list3, ncol = 2)
ggsave(file.path(OUT_FIG, "Bloc4B_UCell_TAMs_MRC1_SPP1_response.png"),
       p_vln_focus, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc4B_UCell_TAMs_MRC1_SPP1_response.png")

# 7) Summary table — mean scores per response group per TAM subtype
message("Generating summary table...")

summary_scores <- seu_TAM@meta.data %>%
  filter(!is.na(final_annotation), !is.na(PathResponse)) %>%
  group_by(PathResponse, final_annotation) %>%
  summarise(across(all_of(score_cols), mean, .names = "mean_{.col}"),
            n_cells = n(),
            .groups = "drop")

print(summary_scores)
fwrite(summary_scores, file.path(OUT_TAB, "Bloc4B_UCell_TAMs_scores_summary.csv"))
message("Saved: Bloc4B_UCell_TAMs_scores_summary.csv")

# 8) Save updated object
message("Saving updated TAMs object with UCell scores...")
saveRDS(seu_TAM, file.path(DATA_DIR_OUTPUT, "Objects/Bloc4B_02_seu_TAMs_UCell.rds"))
message("Saved: Objects/Bloc4B_02_seu_TAMs_UCell.rds")
message("DONE Bloc4B Script 02")