#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 09b: LAG3 expression, CD8 side
# Purpose: verify whether LAG3 (CD8 side of the HLA-D-LAG3
# communication signal) shows a significant MPR vs NMPR
# difference, to determine whether the NMPR-enriched direction
# of this ligand-receptor pair is driven by LAG3 rather than
# by TAM-side HLA-D expression (see Script 09)
# Input:  Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds
# Output: Results/Tables/BLOC5/Bloc5_09b_LAG3_CD8_wilcox.csv
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

# 1) Load CD8 object, restrict to TEX/TPEX
seu_CD8 <- readRDS(file.path(DATA_DIR, "Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds"))
seu_CD8 <- subset(seu_CD8, subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX"))
seu_CD8$cell_type <- seu_CD8$functional.cluster

# 2) Extract LAG3 expression
DefaultAssay(seu_CD8) <- "RNA"
gene_of_interest <- "LAG3"
stopifnot(gene_of_interest %in% rownames(seu_CD8))

expr_vec <- GetAssayData(seu_CD8, layer = "data")[gene_of_interest, ]

meta <- seu_CD8@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, cell_type, Sample, PathResponse)

expr_df <- data.frame(cell_id = names(expr_vec), LAG3 = expr_vec) %>%
  left_join(meta, by = "cell_id")

# 3) Patient-level pseudobulk mean, per CD8 state x patient
pseudobulk <- expr_df %>%
  group_by(cell_type, Sample, PathResponse) %>%
  summarise(LAG3 = mean(LAG3), .groups = "drop")

# 4) Wilcoxon MPR vs NMPR, within each CD8 state
results <- lapply(c("CD8.TEX", "CD8.TPEX"), function(state) {
  mpr_vals  <- pseudobulk %>% filter(cell_type == state, PathResponse == "MPR") %>% pull(LAG3)
  nmpr_vals <- pseudobulk %>% filter(cell_type == state, PathResponse == "NMPR") %>% pull(LAG3)
  if (length(mpr_vals) < 3 | length(nmpr_vals) < 3) {
    return(data.frame(cell_type = state, n_MPR = length(mpr_vals), n_NMPR = length(nmpr_vals),
                      median_MPR = ifelse(length(mpr_vals)>0, median(mpr_vals), NA),
                      median_NMPR = ifelse(length(nmpr_vals)>0, median(nmpr_vals), NA),
                      p_value = NA))
  }
  test <- wilcox.test(mpr_vals, nmpr_vals)
  data.frame(cell_type = state, n_MPR = length(mpr_vals), n_NMPR = length(nmpr_vals),
             median_MPR = median(mpr_vals), median_NMPR = median(nmpr_vals),
             p_value = test$p.value)
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

fwrite(results, file.path(OUT_TAB, "Bloc5_09b_LAG3_CD8_wilcox.csv"))
message("Saved: Bloc5_09b_LAG3_CD8_wilcox.csv")
print(results)