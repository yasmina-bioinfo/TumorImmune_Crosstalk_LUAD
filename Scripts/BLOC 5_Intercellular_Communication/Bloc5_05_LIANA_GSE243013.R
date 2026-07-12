#!/usr/bin/env Rscript
# ============================================================
# GSE243013 : Bloc5 Script 05: LIANA+ intercellular communication
# Infers ligand-receptor interactions across CD8 and TAMs compartments
# Compares interactions across MPR, non-MPR, and pCR conditions
# Methods: sca, natmi, connectome (aggregated consensus)
# CD8 restricted to CD8.TEX/CD8.TPEX only, consistent with Bloc 5
# methodology (GSE207422).
# NOTE: aggregate_rank is a rank-aggregation consensus score across
# methods, NOT a hypothesis-test p-value (Dimitrov et al., 2022).
# A patient-level pseudobulk Wilcoxon test was added on top of
# LIANA's own prioritization score for the MPR vs non-MPR comparison,
# treating each patient as an independent observation, consistent
# with the pseudobulk methodology used throughout this project.
# Input:  Objects/Bloc3_03_seu_Tcells_ProjecTILs.rds
#         Objects/Bloc4A_04_seu_TAMs_annotated.rds
# Output: Results/Figures/BLOC5_Communication/GSE243013/
#         Results/Tables/BLOC5/GSE243013/
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
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC5_Communication/GSE243013")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables/BLOC5/GSE243013")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1) Load objects
# ============================================================
message("Loading objects...")
seu_CD8 <- readRDS(file.path(DATA_DIR, "Objects/Bloc3_03_seu_Tcells_ProjecTILs.rds"))
seu_TAM <- readRDS(file.path(DATA_DIR, "Objects/Bloc4A_04_seu_TAMs_annotated.rds"))

# Restrict CD8 to TEX/TPEX only
seu_CD8 <- subset(seu_CD8, subset = functional.cluster %in% c("CD8.TEX", "CD8.TPEX"))
message("CD8 restricted to TEX/TPEX: ", ncol(seu_CD8), " cells")

# ============================================================
# 2) Unify cell_type column
# ============================================================
seu_CD8$cell_type <- seu_CD8$functional.cluster
seu_TAM$cell_type <- seu_TAM$final_annotation

# ============================================================
# 3) Merge objects
# ============================================================
message("Merging objects...")
seu_TME <- merge(seu_CD8, y = seu_TAM,
                 add.cell.ids = c("CD8", "TAM"),
                 merge.data = TRUE)

# ============================================================
# 4) Remove NA cell types
# ============================================================
seu_TME <- subset(seu_TME, cells = colnames(seu_TME)[!is.na(seu_TME$cell_type)])
message("Total cells after QC: ", ncol(seu_TME))
message("Cell types: ", length(unique(seu_TME$cell_type)))
message("Conditions: ", paste(unique(seu_TME$pathological_response), collapse=", "))

# ============================================================
# 5) Normalize and prepare object
# ============================================================
DefaultAssay(seu_TME) <- "RNA"
seu_TME <- JoinLayers(seu_TME)
seu_TME <- NormalizeData(seu_TME)
Idents(seu_TME) <- "cell_type"

# ============================================================
# 6) Convert to SCE and run LIANA+ - full object (MPR + non-MPR + pCR)
# ============================================================
sce_TME <- as.SingleCellExperiment(seu_TME, assay = "RNA")

message("Running LIANA+ full object (sca + natmi + connectome)...")
liana_results <- liana_wrap(sce_TME,
                            method = c("sca", "natmi", "connectome"),
                            resource = "Consensus",
                            idents_col = "cell_type")

liana_agg <- liana_aggregate(liana_results)
message("LIANA aggregated — ", nrow(liana_agg), " interactions")

fwrite(as.data.frame(liana_agg),
       file.path(OUT_TAB, "Bloc5_03_LIANA_GSE243013_aggregated.csv"))
message("Saved: CSV aggregated interactions — full object")

# ============================================================
# 7) Run LIANA per condition — MPR, non-MPR, pCR separately
# ============================================================
message("Running LIANA per condition...")

liana_per_cond <- list()
for (cond in c("MPR", "non-MPR", "pCR")) {
  sce_sub <- sce_TME[, sce_TME$pathological_response == cond]
  message("Running LIANA for ", cond, " — ", ncol(sce_sub), " cells")
  
  liana_cond <- liana_wrap(sce_sub,
                           method = c("sca", "natmi", "connectome"),
                           resource = "Consensus",
                           idents_col = "cell_type")
  
  liana_cond_agg <- liana_aggregate(liana_cond)
  liana_per_cond[[cond]] <- liana_cond_agg
  
  fwrite(as.data.frame(liana_cond_agg),
         file.path(OUT_TAB, paste0("Bloc5_03_LIANA_GSE243013_",
                                   gsub("-", "_", cond), "_aggregated.csv")))
  message("Saved: LIANA ", cond)
}

# ============================================================
# 8) Formal differential test: CD8(TEX/TPEX) <-> TAM, MPR vs non-MPR,
#    patient-level pseudobulk (optimized: sparse matrix multiplication)
# ============================================================
message("Running patient-level differential test: CD8 <-> TAM, MPR vs non-MPR...")

cd8_states <- c("CD8.TEX", "CD8.TPEX")
tam_states <- unique(seu_TAM$final_annotation)

prioritized <- bind_rows(liana_per_cond[c("MPR", "non-MPR")], .id = "condition") %>%
  filter(aggregate_rank <= 0.05) %>%
  distinct(source, target, ligand.complex, receptor.complex) %>%
  filter((source %in% cd8_states & target %in% tam_states) |
           (target %in% cd8_states & source %in% tam_states))

message("Prioritized CD8<->TAM interactions to test: ", nrow(prioritized))

expr_mat <- GetAssayData(seu_TME, layer = "data")

needed_genes <- unique(unlist(strsplit(
  c(prioritized$ligand.complex, prioritized$receptor.complex), "_")))
needed_genes <- intersect(needed_genes, rownames(expr_mat))
message("Unique genes needed: ", length(needed_genes))

expr_sub <- expr_mat[needed_genes, , drop = FALSE]

meta_TME <- seu_TME@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, cell_type, sampleID, pathological_response) %>%
  filter(pathological_response %in% c("MPR", "non-MPR")) %>%
  mutate(group_key = paste(cell_type, sampleID, sep = "__"))

expr_sub <- expr_sub[, meta_TME$cell_id, drop = FALSE]
stopifnot(identical(colnames(expr_sub), meta_TME$cell_id))

groups <- unique(meta_TME$group_key)
group_idx <- match(meta_TME$group_key, groups)
indicator <- Matrix::sparseMatrix(i = seq_along(group_idx), j = group_idx,
                                  x = 1, dims = c(nrow(meta_TME), length(groups)))
colnames(indicator) <- groups

group_sums   <- expr_sub %*% indicator
group_counts <- Matrix::colSums(indicator)
group_means  <- sweep(as.matrix(group_sums), 2, group_counts, "/")

group_response <- meta_TME %>% distinct(group_key, cell_type, sampleID, pathological_response)

get_scores <- function(gene_complex, ctype) {
  genes <- strsplit(gene_complex, "_")[[1]]
  genes <- intersect(genes, rownames(expr_sub))
  if (length(genes) == 0) return(NULL)
  
  # Use only unique sampleIDs, not the full group_response (which has
  # one row per cell_type x sampleID combination across ALL cell types)
  unique_samples <- group_response %>% distinct(sampleID, pathological_response)
  keys <- paste(ctype, unique_samples$sampleID, sep = "__")
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
  data.frame(sampleID = unique_samples$sampleID, pathological_response = unique_samples$pathological_response,
             mean_expr = vals)
}

message("Testing prioritized CD8<->TAM interactions...")
diff_results <- lapply(seq_len(nrow(prioritized)), function(i) {
  row <- prioritized[i, ]
  lig_score <- get_scores(row$ligand.complex, row$source)
  rec_score <- get_scores(row$receptor.complex, row$target)
  if (is.null(lig_score) || is.null(rec_score)) return(NULL)
  
  merged <- inner_join(lig_score, rec_score, by = c("sampleID", "pathological_response"),
                       suffix = c("_lig", "_rec")) %>%
    filter(!is.na(mean_expr_lig), !is.na(mean_expr_rec)) %>%
    mutate(interaction_score = mean_expr_lig * mean_expr_rec)
  
  mpr    <- merged$interaction_score[merged$pathological_response == "MPR"]
  nonmpr <- merged$interaction_score[merged$pathological_response == "non-MPR"]
  if (length(mpr) < 3 | length(nonmpr) < 3) return(NULL)
  
  test <- wilcox.test(mpr, nonmpr)
  data.frame(source = row$source, target = row$target,
             ligand = row$ligand.complex, receptor = row$receptor.complex,
             p_value = test$p.value,
             median_MPR = median(mpr), median_nonMPR = median(nonmpr),
             n_MPR = length(mpr), n_nonMPR = length(nonmpr))
}) %>% bind_rows()

diff_results <- diff_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

fwrite(diff_results,
       file.path(OUT_TAB, "Bloc5_05_LIANA_GSE243013_differential_wilcox.csv"))
message("Saved: Bloc5_05_LIANA_GSE243013_differential_wilcox.csv (",
        sum(diff_results$p_adj < 0.05, na.rm = TRUE), " significant of ",
        nrow(diff_results), " tested)")

# ============================================================
# 9) Visualization — dot plot top interactions (full object, exploratory)
# ============================================================
message("Generating figures...")

p1 <- liana_agg %>%
  filter(aggregate_rank <= 0.05) %>%
  liana_dotplot(source_groups = c("CD8.TEX", "CD8.TPEX"),
                target_groups = unique(seu_TAM$final_annotation),
                ntop = 20) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 16, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5)) +
  ggtitle("Top CD8 → TAMs interactions — GSE243013 (descriptive prioritization)")
ggsave(file.path(OUT_FIG, "Bloc5_03_LIANA_dotplot_CD8_to_TAMs.png"),
       p1, width = 18, height = 11, dpi = 300, bg = "white")
message("Saved: Dotplot CD8 to TAMs")

# ============================================================
# 10) Save TME object
# ============================================================
saveRDS(seu_TME, file.path(DATA_DIR, "Objects/Bloc5_03_seu_TME_GSE243013.rds"))
message("Saved: Bloc5_03_seu_TME_GSE243013.rds")

message("DONE LIANA+ GSE243013")