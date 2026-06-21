#!/usr/bin/env Rscript
# ============================================================
# GSE207422 : Bloc4B Script 08: SCEVAN , CNV inference on epithelial cells
# MPR and NMPR analyzed separately to reduce RAM requirements
# Ciliated_epithelial used as normal reference in each run
# Input:  Objects/Bloc4B_07_seu_Epithelial.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_08_SCEVAN_*.png
#         Results/Tables/Bloc4B_08_SCEVAN_predictions.csv
#         Objects/Bloc4B_08_seu_Epithelial_SCEVAN.rds
# Reference: De Falco et al., Nature Communications 2023
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SCEVAN)
  library(ggplot2)
  library(dplyr)
  library(data.table)
})

# Paths
DATA_DIR <- "/mnt/c/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc4B_07_seu_Epithelial.rds")
OUT_OBJ  <- file.path(DATA_DIR, "Objects")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load epithelial object
message("Loading epithelial object...")
seu_Epi <- readRDS(IN_OBJ)
message("Total cells: ", ncol(seu_Epi))
DefaultAssay(seu_Epi) <- "RNA"

# 2) Run SCEVAN separately for MPR and NMPR
results_list <- list()

for (resp in c("MPR", "NMPR")) {
  
  message("\n========== Running SCEVAN for ", resp, " ==========")
  
  # Subset cells for this condition (tumor + Ciliated reference)
  cells_resp <- colnames(seu_Epi)[seu_Epi$PathResponse == resp | 
                                   seu_Epi$TME_cell_type == "Ciliated_epithelial"]
  seu_sub <- subset(seu_Epi, cells = cells_resp)
  message("Cells in ", resp, " run: ", ncol(seu_sub))
  message("  Tumor cells: ", sum(seu_sub$PathResponse == resp))
  message("  Ciliated reference: ", sum(seu_sub$TME_cell_type == "Ciliated_epithelial"))
  
  # Count matrix
  count_matrix <- GetAssayData(seu_sub, layer = "counts")
  
  # Normal reference
  normal_cells <- colnames(seu_sub)[seu_sub$TME_cell_type == "Ciliated_epithelial"]
  message("Normal reference cells: ", length(normal_cells))
  
  # Run SCEVAN
  setwd(OUT_FIG)
  results <- pipelineCNA(
    count_matrix,
    norm_cell  = normal_cells,
    par_cores  = 1,
    SUBCLONES  = FALSE,
    plotTree   = FALSE,
    ClonalCN   = FALSE,
    output_dir = file.path(OUT_FIG, "SCEVAN_output")
  )
  
  results$response <- resp
  results_list[[resp]] <- results
  message("SCEVAN done for ", resp)
}

# 3) Combine predictions
message("\nCombining predictions...")
predictions_combined <- do.call(rbind, lapply(names(results_list), function(resp) {
  df <- results_list[[resp]]
  df$cell <- rownames(df)
  df$response <- resp
  df
}))

message("Combined predictions:")
print(table(predictions_combined$class, predictions_combined$response))

# 4) Add to Seurat object
seu_Epi$scevan_class <- predictions_combined$class[match(
  colnames(seu_Epi), predictions_combined$cell)]

message("Predictions by cell type:")
print(table(seu_Epi$scevan_class, seu_Epi$TME_cell_type))

# 5) Save predictions
pred_df <- data.frame(
  cell         = colnames(seu_Epi),
  scevan_class = seu_Epi$scevan_class,
  cell_type    = seu_Epi$TME_cell_type,
  response     = seu_Epi$PathResponse
)
fwrite(pred_df, file.path(OUT_TAB, "Bloc4B_08_SCEVAN_predictions.csv"))
message("Saved: Bloc4B_08_SCEVAN_predictions.csv")

# 6) UMAP by SCEVAN prediction
png(file.path(OUT_FIG, "Bloc4B_08_UMAP_SCEVAN_class.png"),
    width = 10, height = 7, units = "in", res = 300)
print(DimPlot(seu_Epi,
              reduction = "umap",
              group.by  = "scevan_class",
              raster    = FALSE) +
        scale_color_manual(values = c("tumor"  = "#D73027",
                                      "normal" = "#4393C3")) +
        theme_bw() +
        theme(legend.text = element_text(size = 12)) +
        labs(title = "SCEVAN predictions — Epithelial cells GSE207422"))
dev.off()

# 7) Barplot by response
df_bar <- pred_df %>%
  filter(!is.na(scevan_class), !is.na(response)) %>%
  group_by(response, scevan_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n))

png(file.path(OUT_FIG, "Bloc4B_08_Barplot_SCEVAN_response.png"),
    width = 8, height = 6, units = "in", res = 300)
print(ggplot(df_bar, aes(x = response, y = prop, fill = scevan_class)) +
        geom_col(width = 0.8, color = "white") +
        scale_fill_manual(values = c("tumor"  = "#D73027",
                                     "normal" = "#4393C3")) +
        scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
        theme_classic() +
        theme(axis.title.x = element_blank(),
              axis.text.x  = element_text(size = 13, face = "bold"),
              legend.title = element_blank()) +
        ylab("Proportion") +
        labs(title = "Malignant vs normal epithelial cells by response — GSE207422"))
dev.off()

# 8) Save updated object
saveRDS(seu_Epi, file.path(OUT_OBJ, "Bloc4B_08_seu_Epithelial_SCEVAN.rds"))
message("Saved: Objects/Bloc4B_08_seu_Epithelial_SCEVAN.rds")

# 9) Session info
writeLines(capture.output(sessionInfo()),
           file.path(DATA_DIR, "session_info_linux.txt"))
message("DONE Bloc4B Script 08 SCEVAN")