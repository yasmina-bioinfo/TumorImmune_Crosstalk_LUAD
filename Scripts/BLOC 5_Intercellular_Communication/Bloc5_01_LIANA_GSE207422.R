#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc5 Script 01: LIANA+ intercellular communication
# Infers ligand-receptor interactions across CD8, TAMs, and Epithelial compartments
# Compares interactions across MPR and NMPR conditions
# Methods: sca, natmi, connectome (aggregated consensus)
# NOTE: aggregate_rank is a rank-aggregation consensus score across methods,
# NOT a hypothesis-test p-value (Dimitrov et al., 2022; RRA method, Kolde et al., 2012).
# Threshold 0.05 (more permissive than the tutorial's example 0.01), chosen to
# avoid over-restricting discovery in a multi-compartment TME.
# CD8 restricted to CD8.TEX/CD8.TPEX only.
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
  library(Matrix)
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

# Restrict CD8 to TEX/TPEX only
seu_CD8 <- subset(seu_CD8, subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX"))
message("CD8 restricted to TEX/TPEX: ", ncol(seu_CD8), " cells")

# ============================================================
# 2) Unify cell_type column — TAM subtypes use short labels
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
# 7) Run LIANA per condition - MPR and NMPR separately (descriptive)
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

# ============================================================
# 8) Build the two interaction lists first (CD8<->TAM/Epithelial,
#    and TAM/CD8<->tumor/normal), THEN compute expression means ONCE
#    for all genes needed by both, THEN run both statistical tests.
# ============================================================

cd8_states_of_interest <- c("CD8.TEX", "CD8.TPEX")
tam_states_of_interest <- unname(tam_short_labels)
epi_states_of_interest <- c("tumor_MPR", "normal_MPR", "tumor_NMPR", "normal_NMPR", "Ciliated")
valid_partners <- c(tam_states_of_interest, epi_states_of_interest)

# List 1: CD8 <-> TAM / Ciliated (tested directly MPR vs NMPR)
prioritized <- bind_rows(liana_per_cond, .id = "condition") %>%
  filter(aggregate_rank <= 0.05) %>%
  distinct(source, target, ligand.complex, receptor.complex) %>%
  filter((source %in% cd8_states_of_interest & target %in% valid_partners) |
           (target %in% cd8_states_of_interest & source %in% valid_partners))
message("List 1 (CD8<->TAM/Epithelial): ", nrow(prioritized), " interactions")

# List 2: TAM/CD8 <-> tumor/normal (tested as two arms of the same axis)
epi_pairs <- list(
  tumor  = c(mpr = "tumor_MPR",  nmpr = "tumor_NMPR"),
  normal = c(mpr = "normal_MPR", nmpr = "normal_NMPR")
)
epi_partners <- liana_agg %>%
  filter(aggregate_rank <= 0.05) %>%
  filter((source %in% c(tam_states_of_interest, cd8_states_of_interest) &
            target %in% unlist(epi_pairs)) |
           (target %in% c(tam_states_of_interest, cd8_states_of_interest) &
              source %in% unlist(epi_pairs))) %>%
  distinct(source, target, ligand.complex, receptor.complex)
message("List 2 (TAM/CD8<->tumor/normal): ", nrow(epi_partners), " interactions")

# Genes needed for BOTH lists combined
all_complexes <- c(prioritized$ligand.complex, prioritized$receptor.complex,
                   epi_partners$ligand.complex, epi_partners$receptor.complex)
needed_genes <- unique(unlist(strsplit(all_complexes, "_")))

expr_mat <- GetAssayData(seu_TME, layer = "data")
needed_genes <- intersect(needed_genes, rownames(expr_mat))
message("Unique genes needed (both lists): ", length(needed_genes))

expr_sub <- expr_mat[needed_genes, , drop = FALSE]

meta_TME <- seu_TME@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, cell_type, Sample, PathResponse) %>%
  mutate(group_key = paste(cell_type, Sample, sep = "__"))

stopifnot(identical(colnames(expr_sub), meta_TME$cell_id))

groups <- unique(meta_TME$group_key)
group_idx <- match(meta_TME$group_key, groups)
indicator <- Matrix::sparseMatrix(i = seq_along(group_idx), j = group_idx,
                                  x = 1, dims = c(nrow(meta_TME), length(groups)))
colnames(indicator) <- groups

group_sums   <- expr_sub %*% indicator
group_counts <- Matrix::colSums(indicator)
group_means  <- sweep(as.matrix(group_sums), 2, group_counts, "/")

group_response <- meta_TME %>% distinct(group_key, cell_type, Sample, PathResponse)

# get_scores: works for single genes AND multi-subunit complexes (geometric mean)
# FIXED: group_response contains one row per (cell_type, Sample) combination
# across ALL cell types in the merged object (CD8 + TAM + Epithelial), not one
# row per patient. Using it directly duplicated each patient once per cell type
# present elsewhere in the object, inflating n_MPR/n_NMPR far beyond the true
# patient count and producing statistically invalid, artificially tiny p-values.
# Corrected to use a deduplicated list of unique patients before building keys.
get_scores <- function(gene_complex, ctype) {
  genes <- strsplit(gene_complex, "_")[[1]]
  genes <- intersect(genes, rownames(expr_sub))
  if (length(genes) == 0) return(NULL)
  
  unique_samples <- group_response %>% distinct(Sample, PathResponse)
  keys <- paste(ctype, unique_samples$Sample, sep = "__")
  present <- keys %in% colnames(group_means)
  
  if (length(genes) == 1) {
    vals <- rep(NA_real_, length(keys))
    vals[present] <- group_means[genes, keys[present]]
  } else {
    sub_vals <- sapply(genes, function(g) {
      v <- rep(NA_real_, length(keys))
      v[present] <- group_means[g, keys[present]]
      v
    })
    vals <- exp(rowMeans(log1p(sub_vals))) - 1
  }
  
  data.frame(Sample = unique_samples$Sample, PathResponse = unique_samples$PathResponse,
             mean_expr = vals)
}

# ---- Test 1: CD8 <-> TAM/Ciliated, MPR vs NMPR ----
message("Testing List 1 (CD8<->TAM/Epithelial)...")
diff_results <- lapply(seq_len(nrow(prioritized)), function(i) {
  row <- prioritized[i, ]
  lig_score <- get_scores(row$ligand.complex, row$source)
  rec_score <- get_scores(row$receptor.complex, row$target)
  if (is.null(lig_score) || is.null(rec_score)) return(NULL)
  
  merged <- inner_join(lig_score, rec_score, by = c("Sample", "PathResponse"),
                       suffix = c("_lig", "_rec")) %>%
    filter(!is.na(mean_expr_lig), !is.na(mean_expr_rec)) %>%
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
message("Saved List 1 results (", sum(diff_results$p_adj < 0.05, na.rm = TRUE),
        " significant of ", nrow(diff_results), ")")

# ---- Test 2: TAM/CD8 <-> tumor/normal, comparing the two arms ----
message("Testing List 2 (TAM/CD8<->tumor/normal)...")
epi_axis_results <- lapply(seq_len(nrow(epi_partners)), function(i) {
  row <- epi_partners[i, ]
  epi_side <- if (row$source %in% unlist(epi_pairs)) row$source else row$target
  epi_base <- names(epi_pairs)[sapply(epi_pairs, function(x) epi_side %in% x)]
  if (length(epi_base) == 0) return(NULL)
  
  epi_mpr  <- epi_pairs[[epi_base]]["mpr"]
  epi_nmpr <- epi_pairs[[epi_base]]["nmpr"]
  partner  <- if (row$source %in% unlist(epi_pairs)) row$target else row$source
  partner_is_source <- !(row$source %in% unlist(epi_pairs))
  
  score_for <- function(epi_group) {
    lig_ct <- if (partner_is_source) partner else epi_group
    rec_ct <- if (partner_is_source) epi_group else partner
    lig_score <- get_scores(row$ligand.complex, lig_ct)
    rec_score <- get_scores(row$receptor.complex, rec_ct)
    if (is.null(lig_score) || is.null(rec_score)) return(NULL)
    inner_join(lig_score, rec_score, by = c("Sample", "PathResponse"),
               suffix = c("_lig", "_rec")) %>%
      filter(!is.na(mean_expr_lig), !is.na(mean_expr_rec)) %>%
      mutate(interaction_score = mean_expr_lig * mean_expr_rec)
  }
  
  mpr_data  <- score_for(epi_mpr)
  nmpr_data <- score_for(epi_nmpr)
  if (is.null(mpr_data) || is.null(nmpr_data)) return(NULL)
  mpr_data  <- mpr_data  %>% filter(PathResponse == "MPR")
  nmpr_data <- nmpr_data %>% filter(PathResponse == "NMPR")
  if (nrow(mpr_data) < 3 | nrow(nmpr_data) < 3) return(NULL)
  
  test <- wilcox.test(mpr_data$interaction_score, nmpr_data$interaction_score)
  data.frame(partner = partner, epithelial_axis = epi_base,
             ligand = row$ligand.complex, receptor = row$receptor.complex,
             p_value = test$p.value,
             median_MPR_arm = median(mpr_data$interaction_score),
             median_NMPR_arm = median(nmpr_data$interaction_score),
             n_MPR = nrow(mpr_data), n_NMPR = nrow(nmpr_data))
}) %>% bind_rows()

epi_axis_results <- epi_axis_results %>%
  distinct(partner, epithelial_axis, ligand, receptor, .keep_all = TRUE) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

fwrite(epi_axis_results,
       file.path(OUT_TAB, "Bloc5_01_LIANA_GSE207422_TAM_CD8_Epithelial_wilcox.csv"))
message("Saved List 2 results (", sum(epi_axis_results$p_adj < 0.05, na.rm = TRUE),
        " significant of ", nrow(epi_axis_results), ")")

# ============================================================
# 9) Visualizations - dot plot top interactions (descriptive, prioritization-based)
# ============================================================
message("Generating main-text figure (CD8 <-> TAM only)...")

p1 <- liana_agg %>%
  filter(aggregate_rank <= 0.05) %>%
  liana_dotplot(source_groups = cd8_states_of_interest,
                target_groups = tam_states_of_interest,
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