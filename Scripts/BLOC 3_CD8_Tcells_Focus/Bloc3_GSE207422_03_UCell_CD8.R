#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc3 Script 03: UCell scoring on CD8 T cells
# Computes functional gene module scores per cell
# Compares scores across MPR, NMPR
# Focuses on CD8.TEX and CD8.TPEX populations
# Input:  Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds
# Output: Results/Figures/CD8/Bloc3_GSE207422_UCell_*.png
#         Results/Tables/Bloc3_GSE207422_UCell_scores_summary.csv
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
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ <- file.path(DATA_DIR, "Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load CD8 ProjecTILs object
message("Loading CD8 ProjecTILs object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))

# 2) Define gene signatures
# NOTE: exhaustion signature same as portfolio (GSE207422) for cross-dataset comparison
signatures <- list(
  Exhaustion   = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "CTLA4"),
  Cytotoxicity = c("GZMB", "GZMA", "PRF1", "IFNG", "FASLG"),
  TPEX         = c("TCF7", "SLAMF6", "TOX", "PDCD1"),
  Memory       = c("CCR7", "SELL", "IL7R", "TCF7")
)

# 3) Compute UCell scores
# NOTE: UCell computes a per-cell score based on gene set enrichment
# Score range: 0 to 1; higher = stronger signature activity
message("Computing UCell scores...")
# NOTE: JoinLayers required for Seurat v5 + set DefaultAssay to RNA

DefaultAssay(seu_CD8) <- "RNA"
seu_CD8 <- AddModuleScore_UCell(seu_CD8,
                                features = signatures,
                                name     = "_UCell")

message("UCell scores added:")
print(head(seu_CD8@meta.data[, grep("_UCell", colnames(seu_CD8@meta.data))]))

# 4) Violin plots — scores by pathological response
message("Generating violin plots by response...")

score_cols <- paste0(names(signatures), "_UCell")

plot_list <- lapply(score_cols, function(score) {
  VlnPlot(seu_CD8,
          features  = score,
          group.by  = "PathResponse",
          pt.size   = 0,
          cols      = c("MPR" = "#4393C3", "NMPR" = "#D73027")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank())
})

p_vln_response <- wrap_plots(plot_list, ncol = 2)

ggsave(file.path(OUT_FIG, "Bloc3_GSE207422_UCell_violin_response.png"),
       p_vln_response, width = 12, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc3_GSE207422_UCell_violin_response.png")

# 5) Violin plots : scores by ProjecTILs state
message("Generating violin plots by CD8 state...")

plot_list2 <- lapply(score_cols, function(score) {
  VlnPlot(seu_CD8,
          features = score,
          group.by = "functional.cluster",
          pt.size  = 0) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x    = element_blank(),
          axis.text.x     = element_text(angle = 45, hjust = 1))
})

p_vln_state <- wrap_plots(plot_list2, ncol = 2)

ggsave(file.path(OUT_FIG, "Bloc3_GSE207422_UCell_violin_CD8state.png"),
       p_vln_state, width = 12, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc3_GSE207422_UCell_violin_CD8state.png")

# 6) Focus: TEX and TPEX : compare scores across response groups
message("Focus analysis: CD8.TEX and CD8.TPEX...")

seu_TEX_TPEX <- subset(seu_CD8,
                       subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX"))

plot_list3 <- lapply(score_cols, function(score) {
  VlnPlot(seu_TEX_TPEX,
          features = score,
          group.by = "PathResponse",
          split.by = "functional.cluster",
          pt.size  = 0) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_text(angle = 45, hjust = 1))
})

p_vln_focus <- wrap_plots(plot_list3, ncol = 2)

ggsave(file.path(OUT_FIG, "Bloc3_GSE207422_UCell_TEX_TPEX_response.png"),
       p_vln_focus, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved: Bloc3_GSE207422_UCell_TEX_TPEX_response.png")

# 7) Summary table : mean scores per response group per CD8 state
message("Generating summary table...")

summary_scores <- seu_CD8@meta.data %>%
  filter(!is.na(functional.cluster), !is.na(PathResponse)) %>%
  group_by(PathResponse, functional.cluster) %>%
  summarise(across(all_of(score_cols), mean, .names = "mean_{.col}"),
            n_cells = n(),
            .groups = "drop")

print(summary_scores)
fwrite(summary_scores, file.path(OUT_TAB, "Bloc3_GSE207422_UCell_scores_summary.csv"))
message("Saved: Bloc3_GSE207422_UCell_scores_summary.csv")

# 8) Save CD8 object with UCell scores
saveRDS(seu_CD8, file.path(DATA_DIR, "Objects/Bloc3_GSE207422_02_seu_CD8_UCell.rds"))
message("Saved: Objects/Bloc3_GSE207422_02_seu_CD8_UCell.rds")

# 9) Dotplot: all signatures x all CD8 states
Idents(seu_CD8) <- "functional.cluster"

p_dot <- DotPlot(seu_CD8,
                 features = score_cols,
                 group.by = "functional.cluster",
                 cols = c("lightblue", "red"),
                 dot.scale = 8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        plot.title  = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_line(color = "grey90")) +
  ggtitle("UCell scores — CD8 T cells GSE207422") +
  xlab("") + ylab("")

ggsave(file.path(OUT_FIG, "Bloc3_GSE207422_UCell_Dotplot_CD8.png"),
       p_dot, width = 10, height = 6, dpi = 300, bg = "white")
message("Saved: Bloc3_GSE207422_UCell_Dotplot_CD8.png")

# 10) Dotplot x condition
Idents(seu_CD8) <- "PathResponse"

p_dot_response <- DotPlot(seu_CD8,
                          features = score_cols,
                          group.by = "PathResponse",
                          cols = c("lightblue", "red"),
                          dot.scale = 8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        plot.title  = element_text(size = 14, color = "black", face = "bold"),
        panel.grid.major = element_line(color = "grey90")) +
  ggtitle("UCell scores by response — CD8 T cells GSE207422") +
  xlab("") + ylab("")

ggsave(file.path(OUT_FIG, "Bloc3_GSE207422_UCell_Dotplot_CD8_response.png"),
       p_dot_response, width = 8, height = 5, dpi = 300, bg = "white")

# 11) Dotplot x CD8 states x condition
p_dot_split <- DotPlot(seu_CD8,
                       features = score_cols,
                       group.by = "functional.cluster",
                       split.by = "PathResponse",
                       cols = c("MPR" = "#E63946", "NMPR" = "#457B9D"),
                       dot.scale = 6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        plot.title  = element_text(size = 13, color = "black", face = "bold"),
        panel.grid.major = element_line(color = "grey90")) +
  ggtitle("UCell scores by CD8 state and response — GSE207422") +
  labs(subtitle = "Red = MPR | Blue = NMPR") +
  xlab("") + ylab("")

ggsave(file.path(OUT_FIG, "Bloc3_GSE207422_UCell_Dotplot_CD8_split.png"),
       p_dot_split, width = 14, height = 7, dpi = 300, bg = "white")

message("DONE Bloc3_GSE207422 Script 03")