#!/usr/bin/env Rscript
# ============================================================
# GSE243013: Bloc3 Script 09b: Barplot CD8 ProjecTILs states
# MPR vs non-MPR only — pCR excluded from main preprint narrative
# Input:  Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
# Output: Results/Figures/CD8/Bloc3_Barplot_CD8_ProjecTILs_MPR_nonMPR.png
#         Results/Tables/Bloc3_CD8_chisq_MPR_nonMPR.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(data.table)
})

DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load and subset MPR + non-MPR only
message("Loading CD8 ProjecTILs object...")
seu_CD8 <- readRDS(IN_OBJ)
seu_CD8 <- subset(seu_CD8, subset = pathological_response %in% c("MPR", "non-MPR"))
message("Cells after pCR exclusion: ", ncol(seu_CD8))

# 2) Color palette
cd8_colors <- c(
  "CD8.TEX"       = "#D73027",
  "CD8.TPEX"      = "#FC8D59",
  "CD8.EM"        = "#4393C3",
  "CD8.CM"        = "#2166AC",
  "CD8.NaiveLike" = "#74ADD1",
  "CD8.TEMRA"     = "#762A83",
  "CD8.MAIT"      = "#A6D96A"
)

# 3) Compute proportions
df_prop <- seu_CD8@meta.data %>%
  filter(!is.na(functional.cluster), !is.na(pathological_response)) %>%
  group_by(pathological_response, functional.cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pathological_response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 4) Chi-2 test
cont_table <- table(seu_CD8$functional.cluster, seu_CD8$pathological_response)
chisq_res  <- chisq.test(cont_table)
message("Chi-2 p-value: ", chisq_res$p.value)

fwrite(data.frame(statistic = chisq_res$statistic,
                  df        = chisq_res$parameter,
                  p.value   = chisq_res$p.value),
       file.path(OUT_TAB, "Bloc3_CD8_chisq_MPR_nonMPR.csv"))
message("Saved: Chi-2 table")

# 5) Barplot
p_bar <- ggplot(df_prop,
                aes(x = pathological_response, y = prop, fill = functional.cluster)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = cd8_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.1)) +
  annotate("text",
           x = 1.5, y = 1.08,
           label    = paste0("Chi-2 p = ", signif(chisq_res$p.value, 3)),
           size     = 3.5,
           fontface = "italic") +
  ylab("CD8 state proportion") +
  labs(caption = "Note: proportions based on ProjecTILs functional.cluster annotation (Script 08)") +
  theme_classic() +
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_text(size = 13, face = "bold"),
        axis.text.y     = element_text(size = 11),
        axis.title.y    = element_text(size = 12),
        legend.position = "right",
        legend.text     = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.title    = element_blank(),
        plot.caption    = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Bloc3_Barplot_CD8_ProjecTILs_MPR_nonMPR.png"),
       p_bar, width = 7, height = 6, dpi = 300, bg = "white")
message("Saved: Bloc3_Barplot_CD8_ProjecTILs_MPR_nonMPR.png")

message("DONE Bloc3 Script 09b")