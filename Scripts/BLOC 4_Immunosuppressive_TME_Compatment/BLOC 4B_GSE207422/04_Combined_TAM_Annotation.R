#!/usr/bin/env Rscript
# ============================================================
# GSE207422 — Bloc4B Script 04: Combined TAM annotation + UMAP + Barplot
# Merges TAM_like subclusters (Script 03) with TAM_like_MRC1 and TAM_like_SPP1
# Final annotation: 9 TAM subtypes (cluster 7 TAM_like excluded — T cell contamination)
# Input:  Objects/Bloc4B_01_seu_TAMs_annotated.rds
#         Objects/Bloc4B_03_seu_TAMlike_annotated.rds
# Output: Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_UMAP_TAMs_combined.png
#         Results/Figures/BLOC4B_Epithelial_TAMs/Bloc4B_Barplot_TAMs_combined.png
#         Results/Tables/BLOC4B/Bloc4B_TAMs_combined_chisq.csv
#         Objects/Bloc4B_04_seu_TAMs_combined.rds
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(scales)
})

# Paths
DATA_DIR_OUTPUT <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ_TAM     <- file.path(DATA_DIR_OUTPUT, "Objects/Bloc4B_01_seu_TAMs_annotated.rds")
IN_OBJ_TAMLIKE <- file.path(DATA_DIR_OUTPUT, "Objects/Bloc4B_03_seu_TAMlike_annotated.rds")
OUT_FIG <- file.path(DATA_DIR_OUTPUT, "Results/Figures/BLOC4B_Epithelial_TAMs")
OUT_TAB <- file.path(DATA_DIR_OUTPUT, "Results/Tables/BLOC4B")
OUT_OBJ <- file.path(DATA_DIR_OUTPUT, "Objects")

# 1) Load objects
message("Loading objects...")
seu_TAM     <- readRDS(IN_OBJ_TAM)
seu_TAMlike <- readRDS(IN_OBJ_TAMLIKE)
message("TAMs total: ", ncol(seu_TAM))
message("TAM_like reclustered: ", ncol(seu_TAMlike))

# 2) Transfer detailed annotation from TAM_like to main object
# NOTE: TAM_like cells get their detailed annotation from Script 03
#       TAM_like_MRC1 and TAM_like_SPP1 keep their original annotation
message("Transferring TAM_like detailed annotations...")

# Get barcodes of TAM_like cells
tamlike_barcodes <- colnames(seu_TAMlike)

# Get detailed annotations from TAM_like reclustering
tamlike_annot <- seu_TAMlike@meta.data %>%
  select(final_annotation_detailed) %>%
  tibble::rownames_to_column("barcode")

# Build combined annotation vector
combined_annot <- seu_TAM@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  select(barcode, final_annotation) %>%
  left_join(tamlike_annot, by = "barcode") %>%
  mutate(combined_annotation = case_when(
    !is.na(final_annotation_detailed) ~ final_annotation_detailed,
    TRUE ~ final_annotation
  )) %>%
  pull(combined_annotation)

names(combined_annot) <- colnames(seu_TAM)
seu_TAM$combined_annotation <- combined_annot

message("Combined annotation distribution:")
print(table(seu_TAM$combined_annotation))

# 3) Exclude T cell contamination
seu_plot <- subset(seu_TAM,
                   subset = combined_annotation != "EXCLUDE — T cell contamination")
message("Cells after exclusion: ", ncol(seu_plot))

# 4) Color palette  9 TAM subtypes
tam_colors <- c(
  "TAM_like_MRC1"                                    = "#A50026",
  "TAM_like_SPP1"                                    = "#1A7A1A",
  "TAM_like_resident_M2 (iron metabolism/anti-inflammatory)" = "#D73027",
  "TAM_like_IFN (PD-L1+/IDO1+/CXCL9+)"             = "#4575B4",
  "TAM_like_monocyte (classical inflammatory)"        = "#FDAE61",
  "TAM_like_lipid (CCL18+/AKR+)"                    = "#FEE090",
  "TAM_like_stress (HSP-high/M1-like)"               = "#F46D43",
  "TAM_like_regulatory (glucocorticoid-responsive)"  = "#ABD9E9",
  "TAM_like_M2 (SIGLEC8+/CCL18+)"                   = "#74ADD1"
)

# 5) UMAP  combined annotation
message("Generating UMAP combined annotation...")
png(file.path(OUT_FIG, "Bloc4B_UMAP_TAMs_combined.png"),
    width = 16, height = 10, units = "in", res = 300)
print(DimPlot(seu_plot,
              group.by   = "combined_annotation",
              cols       = tam_colors,
              label      = FALSE,
              pt.size    = 0.3,
              reduction  = "umap") +
        labs(title = "TAM subtypes combined — GSE207422 NSCLC") +
        theme(plot.title      = element_text(size = 14, face = "bold"),
              legend.text     = element_text(size = 12),
              legend.key.size = unit(0.6, "cm"),
              legend.title    = element_blank()))

dev.off()
message("Saved: Bloc4B_UMAP_TAMs_combined.png")

# 6) UMAP split by response
message("Generating UMAP split by response...")
png(file.path(OUT_FIG, "Bloc4B_UMAP_TAMs_combined_split.png"),
    width = 18, height = 8, units = "in", res = 300)
print(DimPlot(seu_plot,
              group.by   = "combined_annotation",
              split.by   = "PathResponse",
              cols       = tam_colors,
              label      = FALSE,
              pt.size    = 0.3,
              reduction  = "umap") +
        labs(title = "TAM subtypes split by pathological response — GSE207422 NSCLC") +
        theme(plot.title      = element_text(size = 14, face = "bold"),
              legend.text     = element_text(size = 12),
              legend.key.size = unit(0.6, "cm"),
              legend.title    = element_blank(),
              strip.text      = element_text(size = 13, face = "bold")))

dev.off()
message("Saved: Bloc4B_UMAP_TAMs_combined_split.png")

# 7) Barplot proportions
message("Computing proportions...")
df <- seu_plot@meta.data %>%
  filter(!is.na(combined_annotation), !is.na(PathResponse)) %>%
  transmute(response  = PathResponse,
            cell_type = combined_annotation)

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
                     limits = c(0, 1.05)) +
  ylab("TAM subtype proportion") +
  theme_classic() +
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_text(size = 13, face = "bold"),
        axis.text.y     = element_text(size = 11),
        axis.title.y    = element_text(size = 12),
        legend.position = "right",
        legend.text     = element_text(size = 11),
        legend.key.size = unit(0.5, "cm"),
        legend.title    = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Bloc4B_Barplot_TAMs_combined.png"),
       p_bar, width = 10, height = 6, dpi = 300, bg = "white")
message("Saved: Bloc4B_Barplot_TAMs_combined.png")

# 8) Chi-2 test
message("Running statistical tests...")
cont_table <- table(seu_plot$combined_annotation, seu_plot$PathResponse)
chisq <- chisq.test(cont_table)
message("Chi-2 p-value: ", chisq$p.value)

chisq_summary <- data.frame(
  test    = "Chi-squared",
  p.value = chisq$p.value,
  df      = chisq$parameter
)
fwrite(chisq_summary, file.path(OUT_TAB, "Bloc4B_TAMs_combined_chisq.csv"))
message("Saved: Bloc4B_TAMs_combined_chisq.csv")

# 9) Save combined object
message("Saving combined object...")
saveRDS(seu_plot,
        file.path(OUT_OBJ, "Bloc4B_04_seu_TAMs_combined.rds"))
message("Saved: Objects/Bloc4B_04_seu_TAMs_combined.rds")
message("DONE Bloc4B Script 04")