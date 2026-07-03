#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 05b: UCell Hallmark scoring on TAMs
# Unbiased approach — all 50 Hallmark signatures tested
# Discriminant signatures identified by Wilcoxon MPR vs NMPR
# Input:  Objects/Bloc4B_05_seu_TAMs_combined_UCell.rds
# Output: Results/Tables/GSE207422_UCell_Hallmark_TAMs_wilcox_results.csv
#         Objects/GSE207422_seu_TAMs_Hallmark.rds
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(msigdbr)
  library(dplyr)
  library(data.table)
  library(BiocParallel)
})

register(SnowParam(workers = 1))

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_05_seu_TAMs_combined_UCell.rds")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load TAMs object
message("Loading TAMs object...")
seu_TAM <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_TAM))
DefaultAssay(seu_TAM) <- "RNA"

# 2) Load ALL Hallmark signatures
message("Loading Hallmark signatures...")
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)
hallmark_cols <- names(hallmark_sets)
message("Total Hallmark signatures: ", length(hallmark_sets))

# 3) Compute UCell scores
message("Computing UCell Hallmark scores...")
seu_TAM <- AddModuleScore_UCell(seu_TAM,
                                features = hallmark_sets,
                                name     = "")
message("Hallmark UCell scores added.")

# 4) Wilcoxon test — MPR vs NMPR
message("Running Wilcoxon tests MPR vs NMPR...")

wilcox_results <- lapply(hallmark_cols, function(sig) {
  mpr  <- seu_TAM@meta.data[seu_TAM@meta.data$PathResponse == "MPR",  sig]
  nmpr <- seu_TAM@meta.data[seu_TAM@meta.data$PathResponse == "NMPR", sig]
  mpr  <- mpr[!is.na(mpr)]
  nmpr <- nmpr[!is.na(nmpr)]
  if (length(mpr) < 3 | length(nmpr) < 3) return(NULL)
  test <- wilcox.test(mpr, nmpr)
  data.frame(signature     = sig,
             p_value       = test$p.value,
             median_MPR    = median(mpr),
             median_NMPR   = median(nmpr))
})

wilcox_results <- do.call(rbind, wilcox_results)

wilcox_results <- wilcox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

fwrite(wilcox_results,
       file.path(OUT_TAB, "GSE207422_UCell_Hallmark_TAMs_wilcox_results.csv"))
message("Saved: GSE207422_UCell_Hallmark_TAMs_wilcox_results.csv")

# 5) Discriminant signatures
sig_discriminant <- wilcox_results %>%
  filter(p_adj < 0.05) %>%
  pull(signature)

message("Discriminant Hallmark signatures (p_adj < 0.05): ",
        length(sig_discriminant))
message(paste(sig_discriminant, collapse = ", "))

# 6) Save object
message("Saving TAMs object with Hallmark scores...")
saveRDS(seu_TAM,
        file.path(DATA_DIR, "Objects/GSE207422_seu_TAMs_Hallmark.rds"))
message("Saved: Objects/GSE207422_seu_TAMs_Hallmark.rds")

message("DONE Script 05b")