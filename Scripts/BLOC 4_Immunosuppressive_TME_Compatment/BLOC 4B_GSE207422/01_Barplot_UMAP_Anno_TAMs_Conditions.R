#!/usr/bin/env Rscript
# ============================================================
# # GSE207422 Bloc4B Script 01: TAMs annotation + UMAP + Barplot
# NOTE: TAM annotation already performed during global TME annotation
# (TME_cell_type column in 04_TME_MPR_NMPR.rds)
# Three TAM subtypes identified: TAM_like, TAM_like_MRC1, TAM_like_SPP1
# TAM_like_SPP1 small population (n=362), interpret with caution
# Input:  objects/04_TME_MPR_NMPR.rds (portfolio repo)
# Output: Results/Figures/BLOC4B_TAMs/
#         Results/Tables/BLOC4B/
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(scales)
})

# Paths
DATA_DIR_PORTFOLIO <- "C:/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Immunotherapy"
DATA_DIR_OUTPUT    <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"

IN_OBJ  <- file.path(DATA_DIR_PORTFOLIO, "objects/04_TME_MPR_NMPR.rds")
OUT_OBJ <- file.path(DATA_DIR_OUTPUT, "Objects")
OUT_FIG <- file.path(DATA_DIR_OUTPUT, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB <- file.path(DATA_DIR_OUTPUT, "Results/Tables/BLOC4B")

# 1) Load Seurat object
message("Loading Seurat object...")
seu <- readRDS(IN_OBJ)
seu_TAM <- subset(seu, subset = TME_cell_type %in% c("TAM_like", "TAM_like_MRC1", "TAM_like_SPP1"))
rm(seu); gc()
message("TAMs: ", ncol(seu_TAM))

# 2) NOTE: annotation already in TME_cell_type column, use directly
seu_TAM$final_annotation <- seu_TAM$TME_cell_type
message("Annotation distribution:")
print(table(seu_TAM$final_annotation))


# 3) Define color palette
# Colors organized by cell lineage for visual coherence
tam_colors <- c(
  "TAM_like"       = "#D73027",
  "TAM_like_MRC1"  = "#A50026",
  "TAM_like_SPP1"  = "#4575B4"
)

# Ensure factor order matches palette
seu_TAM$final_annotation <- factor(seu_TAM$final_annotation, levels = names(tam_colors))

# 4) UMAP : final annotation
message("Generating UMAP final annotation...")
png(file.path(OUT_FIG, "Bloc4B_UMAP_TAMs_annotation.png"),
    width = 12, height = 8, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              group.by  = "final_annotation",
              cols      = tam_colors,
              label     = FALSE,
              pt.size   = 0.1,
              raster    = FALSE) +
        theme_classic() +
        labs(title = "TAM subtypes — GSE207422 NSCLC") +
        theme(plot.title       = element_text(size = 14, face = "bold"),
              legend.text      = element_text(size = 14),
              legend.key.size  = unit(0.6, "cm"),
              legend.title     = element_blank()))
dev.off()
message("Saved: Bloc4B_UMAP_TAMs_annotation.png")

# UMAP split by pathological response
png(file.path(OUT_FIG, "Bloc4B_UMAP_TAMs_split_response.png"),
    width = 18, height = 10, units = "in", res = 300)
print(DimPlot(seu_TAM,
              reduction = "umap",
              group.by  = "final_annotation",
              split.by  = "PathResponse",
              cols      = tam_colors,
              label     = FALSE,
              pt.size   = 0.1,
              raster    = FALSE) +
        theme_classic() +
        labs(title = "TAM subtypes split by pathological response") +
        theme(plot.title      = element_text(size = 13, face = "bold"),
              strip.text      = element_text(size = 12, face = "bold"),
              legend.text     = element_text(size = 14),
              legend.key.size = unit(0.7, "cm"),
              legend.title    = element_blank()))
dev.off()
message("Saved: Bloc4B_UMAP_TAMs_split_response.png")

# 5) Barplot : proportions by response group
message("Computing proportions...")

df <- seu_TAM@meta.data %>%
  filter(!is.na(final_annotation), 
         !is.na(PathResponse)) %>%
  transmute(response  = PathResponse,
            cell_type = final_annotation)

df_prop <- df %>%
  count(response, cell_type, name = "n") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar <- ggplot(df_prop,
                aes(x = response, y = prop, fill = cell_type)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = tam_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.1)) +
  annotate("text",
           x        = 2,
           y        = 1.08,
           label    = "Chi-2 p < 2.2e-16",
           size     = 3.5,
           fontface = "italic") +
  ylab("TAM subtype proportion") +
  theme_classic() +
  theme(axis.title.x  = element_blank(),
        axis.text.x   = element_text(size = 13, face = "bold"),
        axis.text.y   = element_text(size = 11),
        axis.title.y  = element_text(size = 12),
        legend.position = "right",
        legend.text   = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"),
        legend.title  = element_blank(),
        plot.caption  = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Bloc4B_Barplot_TAMs_proportions.png"),
       p_bar, width = 8, height = 7, dpi = 300, bg = "white")
message("Saved.")
message("Saved: Bloc4B_Barplot_TAMs_proportions.png")

# 6) Statistical tests
message("Running statistical tests...")

# Contingency table: cell types x response groups
cont_table <- table(seu_TAM$final_annotation, seu_TAM$PathResponse)

# Chi-2 global test
chisq_res <- chisq.test(cont_table)
message("Chi-2 test p-value: ", chisq_res$p.value)

chisq_summary <- data.frame(
  statistic = chisq_res$statistic,
  df        = chisq_res$parameter,
  p.value   = chisq_res$p.value
)
fwrite(chisq_summary, file.path(OUT_TAB, "Bloc4B_TAMs_chisq_test.csv"))
message("Saved: Bloc4B_TAMs_chisq_test.csv")

# Fisher post-hoc — GSE207422 has only MPR and NMPR (no pCR)
comparisons <- list(
  c("MPR", "NMPR")
)

fisher_results <- lapply(comparisons, function(pair) {
  sub_table <- cont_table[, pair]
  test <- fisher.test(sub_table, simulate.p.value = TRUE, B = 10000)
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    p.value    = test$p.value,
    p.bonf     = min(test$p.value * length(comparisons), 1)
  )
})

fisher_df <- do.call(rbind, fisher_results)
print(fisher_df)
fwrite(fisher_df, file.path(OUT_TAB, "Bloc4B_TAMs_fisher_posthoc.csv"))
message("Saved: Bloc4B_TAMs_fisher_posthoc.csv")

# 7) Save updated Seurat object with final_annotation
message("Saving updated Seurat object...")
saveRDS(seu_TAM, file.path(OUT_OBJ, "Bloc4B_01_seu_TAMs_annotated.rds"))
message("Saved: Objects/Bloc4B_01_seu_TAMs_annotated.rds")
message("DONE Bloc2 Script 01")