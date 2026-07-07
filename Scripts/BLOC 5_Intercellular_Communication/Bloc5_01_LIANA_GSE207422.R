#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 01: LIANA+ intercellular communication
# Infers ligand-receptor interactions across CD8, TAMs, and Epithelial compartments
# Compares interactions across MPR and NMPR conditions
# Methods: sca, natmi, connectome (aggregated consensus)
# NOTE: aggregate_rank is a rank-aggregation consensus score across methods,
# NOT a hypothesis-test p-value. Differential significance between conditions
# is assessed separately (Step 8) via patient-level pseudobulk Wilcoxon testing,
# consistent with the pseudobulk methodology used elsewhere in this project
# (avoids pseudo-replication at the single-cell level).
# Epithelial compartment included here for completeness (available in GSE207422
# only); CD8<->TAM interactions are the main-text focus, CD8/TAM<->Epithelial
# results are reserved for supplementary figures / perspectives.
# Input:  Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds
#         Objects/Bloc4B_04_seu_TAMs_combined.rds
#         Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds
# Output: Results/Figures/BLOC5_Communication/GSE207422/
#         Results/Tables/BLOC5/
# Reference: Dimitrov et al., Nature Communications 2022 (LIANA+)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(liana)
  library(SingleCellExperiment)
  library(tidyverse)
  library(ggplot2)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE207422")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load objects
# ============================================================
message("Loading objects...")
seu_CD8 <- readRDS(file.path(DATA_DIR, "Objects/Bloc3_GSE207422_01_seu_CD8_ProjecTILs.rds"))
seu_TAM <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_04_seu_TAMs_combined.rds"))
seu_Epi <- readRDS(file.path(DATA_DIR, "Objects/Bloc4B_09c_seu_Epithelial_FinalGroups.rds"))

# ============================================================
# 2) Unify cell_type column — TAM subtypes use short labels,
#    consistent with tam_short convention used in Blocs 3/4
# ============================================================
seu_CD8$cell_type <- seu_CD8$functional.cluster

tam_short_labels <- c(
  "TAM_like_MRC1"                                            = "MRC1+ M2-like",
  "TAM_like_SPP1"                                            = "SPP1+ immunosuppressive",
  "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)" = "Resident M2",
  "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)"                     = "IFN-stimulated",
  "TAM_like_monocyte (classical inflammatory)"               = "Monocyte-derived",
  "TAM_like_lipid (CCL18+/AKR+)"                            = "Lipid-associated",
  "TAM_like_stress (HSP-high/M1-like)"                      = "Stress-response",
  "TAM_like_regulatory (glucocorticoid-responsive)"          = "Regulatory",
  "TAM_like_M2 (SIGLEC8+/CCL18+)"                           = "M2-SIGLEC8+"
)
seu_TAM$cell_type <- unname(tam_short_labels[as.character(seu_TAM$combined_annotation)])
if (any(is.na(seu_TAM$cell_type))) {
  stop("Unmapped TAM subtype label(s): ",
       paste(unique(seu_TAM$combined_annotation[is.na(seu_TAM$cell_type)]), collapse = ", "))
}

seu_Epi$cell_type <- seu_Epi$final_group

# ============================================================
# 3) Merge objects
# ============================================================
message("Merging objects...")
seu_TME <- merge(seu_CD8, y = list(seu_TAM, seu_Epi),
                 add.cell.ids = c("CD8", "TAM", "Epi"),
                 merge.data = TRUE)

# ============================================================
# 4) Remove NA cell types
# ============================================================
seu_TME <- subset(seu_TME, cells = colnames(seu_TME)[!is.na(seu_TME$cell_type)])
message("Total cells after QC: ", ncol(seu_TME))
message("Cell types: ", length(unique(seu_TME$cell_type)))
message("Conditions: ", paste(unique(seu_TME$PathResponse), collapse=", "))

# ============================================================
# 5) Normalize and prepare object
# ============================================================
DefaultAssay(seu_TME) <- "RNA"
seu_TME <- JoinLayers(seu_TME)
seu_TME <- NormalizeData(seu_TME)
Idents(seu_TME) <- "cell_type"

# ============================================================
# 6) Convert to SCE and run LIANA+ - 3 methods aggregated (pooled, both conditions)
# ============================================================
sce_TME <- as.SingleCellExperiment(seu_TME, assay = "RNA")

message("Running LIANA+ (sca + natmi + connectome) — pooled...")
liana_results <- liana_wrap(sce_TME,
                            method = c("sca", "natmi", "connectome"),
                            resource = "Consensus",
                            idents_col = "cell_type")

liana_agg <- liana_aggregate(liana_results)
message("LIANA aggregated — ", nrow(liana_agg), " interactions")

fwrite(as.data.frame(liana_agg),
       file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_aggregated.csv"))
message("Saved: CSV aggregated interactions (pooled)")

# ============================================================
# 7) Run LIANA per condition - MPR and NMPR separately
#    (descriptive/prioritization use only — NOT a statistical test,
#    see Step 8 for the formal differential comparison)
# ============================================================
message("Running LIANA per condition...")

liana_per_cond <- list()
for (cond in c("MPR", "NMPR")) {
  sce_sub <- sce_TME[, sce_TME$PathResponse == cond]
  message("Running LIANA for ", cond, " — ", ncol(sce_sub), " cells")
  
  liana_cond <- liana_wrap(sce_sub,
                           method = c("sca", "natmi", "connectome"),
                           resource = "Consensus",
                           idents_col = "cell_type")
  
  liana_cond_agg <- liana_aggregate(liana_cond)
  liana_per_cond[[cond]] <- liana_cond_agg
  
  fwrite(as.data.frame(liana_cond_agg),
         file.path(OUT_TAB, paste0("Bloc5_01_LIANA_GSE207422_", cond, "_aggregated.csv")))
  message("Saved: LIANA ", cond, " (descriptive)")
}

# 7b) Additional Wilcoxon tests: tumor vs normal, within each response condition
#     (malignancy effect, tested separately per condition rather than pooled,
#     since tumor vs normal differences may not be uniform across MPR/NMPR)
message("Running Wilcoxon tests: tumor vs normal, within each condition...")

run_malignancy_comparison <- function(scores, meta, level_a, level_b, label) {
  full <- scores %>% left_join(meta, by = "condition") %>%
    filter(final_group %in% c(level_a, level_b), source %in% top20_tfs)
  results <- lapply(top20_tfs, function(tf) {
    a <- full$score[full$final_group == level_a & full$source == tf]
    b <- full$score[full$final_group == level_b & full$source == tf]
    if (length(a) < 3 | length(b) < 3) return(NULL)
    test <- wilcox.test(a, b)
    data.frame(comparison = label, source = tf, p_value = test$p.value,
               median_A = median(a), median_B = median(b))
  }) %>% bind_rows()
  names(results)[names(results) == "median_A"] <- paste0("median_", level_a)
  names(results)[names(results) == "median_B"] <- paste0("median_", level_b)
  results
}

res_MPR   <- run_malignancy_comparison(tf_scores, meta, "tumor_MPR",  "normal_MPR",  "tumor_vs_normal_MPR")
res_NMPR  <- run_malignancy_comparison(tf_scores, meta, "tumor_NMPR", "normal_NMPR", "tumor_vs_normal_NMPR")

malignancy_results <- bind_rows(res_MPR, res_NMPR) %>%
  group_by(comparison) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(comparison, p_adj)

fwrite(malignancy_results,
       file.path(OUT_TAB, "Bloc4B_11_CollecTRI_Epithelial_wilcox_malignancy.csv"))
message("Saved: Bloc4B_11_CollecTRI_Epithelial_wilcox_malignancy.csv")

for (comp in unique(malignancy_results$comparison)) {
  n_sig <- sum(malignancy_results$comparison == comp & malignancy_results$p_adj < 0.05, na.rm = TRUE)
  message(comp, ": ", n_sig, "/", length(top20_tfs), " TFs significant (p_adj < 0.05)")
}

# ============================================================
# 8) Formal differential test: MPR vs NMPR, patient-level pseudobulk
#    For each prioritized interaction (aggregate_rank <= 0.05 in at least one
#    condition), compute a per-patient interaction score as the product of
#    mean ligand expression (source cell type) and mean receptor expression
#    (target cell type), then Wilcoxon rank-sum test across patients,
#    BH-corrected across all tested interactions.
# ============================================================
message("Running patient-level differential test (MPR vs NMPR)...")

prioritized <- bind_rows(liana_per_cond, .id = "condition") %>%
  filter(aggregate_rank <= 0.05) %>%
  distinct(source, target, ligand.complex, receptor.complex)

message("Prioritized interactions to test: ", nrow(prioritized))

expr_mat <- GetAssayData(seu_TME, layer = "data")
meta_TME <- seu_TME@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, cell_type, sampleID, PathResponse)

get_patient_score <- function(gene, ctype) {
  # Mean expression of `gene` in cells of type `ctype`, per patient
  cells_ct <- meta_TME %>% filter(cell_type == ctype)
  if (!(gene %in% rownames(expr_mat)) || nrow(cells_ct) == 0) return(NULL)
  vals <- expr_mat[gene, cells_ct$cell_id]
  data.frame(cell_id = cells_ct$cell_id, expr = vals) %>%
    left_join(cells_ct, by = "cell_id") %>%
    group_by(sampleID, PathResponse) %>%
    summarise(mean_expr = mean(expr), .groups = "drop")
}

diff_results <- lapply(seq_len(nrow(prioritized)), function(i) {
  row <- prioritized[i, ]
  # simple single-subunit case; multi-subunit complexes (e.g. "A_B") use the
  # first subunit as a practical proxy score, consistent with common LIANA
  # downstream conventions
  lig_gene <- strsplit(row$ligand.complex, "_")[[1]][1]
  rec_gene <- strsplit(row$receptor.complex, "_")[[1]][1]
  
  lig_score <- get_patient_score(lig_gene, row$source)
  rec_score <- get_patient_score(rec_gene, row$target)
  if (is.null(lig_score) || is.null(rec_score)) return(NULL)
  
  merged <- inner_join(lig_score, rec_score, by = c("sampleID", "PathResponse"),
                       suffix = c("_lig", "_rec")) %>%
    mutate(interaction_score = mean_expr_lig * mean_expr_rec)
  
  mpr  <- merged$interaction_score[merged$PathResponse == "MPR"]
  nmpr <- merged$interaction_score[merged$PathResponse == "NMPR"]
  if (length(mpr) < 3 | length(nmpr) < 3) return(NULL)
  
  test <- wilcox.test(mpr, nmpr)
  data.frame(source = row$source, target = row$target,
             ligand = row$ligand.complex, receptor = row$receptor.complex,
             p_value = test$p.value,
             median_MPR = median(mpr), median_NMPR = median(nmpr),
             n_MPR = length(mpr), n_NMPR = length(nmpr))
}) %>% bind_rows()

diff_results <- diff_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

fwrite(diff_results,
       file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_differential_wilcox.csv"))
message("Saved: Bloc5_01_LIANA_GSE207422_differential_wilcox.csv (",
        sum(diff_results$p_adj < 0.05, na.rm = TRUE), " significant of ",
        nrow(diff_results), " tested)")

# ============================================================
# 9) Visualizations - dot plot top interactions (descriptive, prioritization-based)
#    NOTE: figure generation split by target audience —
#    CD8<->TAM only for main text (this plot); CD8/TAM<->Epithelial
#    reserved for supplementary (separate script)
# ============================================================
message("Generating main-text figure (CD8 <-> TAM only)...")

tam_subtypes <- unname(tam_short_labels)

p1 <- liana_agg %>%
  filter(aggregate_rank <= 0.05) %>%
  liana_dotplot(source_groups = c("CD8.TEX", "CD8.TPEX"),
                target_groups = tam_subtypes,
                ntop = 20) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5)) +
  ggtitle("Top CD8 <-> TAM interactions — GSE207422 (descriptive prioritization)")

ggsave(file.path(OUT_FIG, "Bloc5_01_LIANA_dotplot_CD8_to_TAM_main.png"),
       p1, width = 18, height = 11, dpi = 300, bg = "white")
message("Saved: Dotplot CD8 <-> TAM (main text)")

# ============================================================
# 10) Save TME object
# ============================================================
saveRDS(seu_TME, file.path(DATA_DIR, "Objects/Bloc5_01_seu_TME_GSE207422.rds"))
message("Saved: Bloc5_01_seu_TME_GSE207422.rds")

message("DONE LIANA+ GSE207422")