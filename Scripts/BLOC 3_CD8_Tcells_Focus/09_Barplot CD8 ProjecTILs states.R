#!/usr/bin/env Rscript
# ============================================================
# GSE243013: Bloc3 Script 09: Barplot CD8 ProjecTILs states
# Proportions of CD8 functional states by pathological response
# Input:  Objects/Bloc3_08_seu_CD8_ProjecTILs.rds
# Output: Results/Figures/CD8/Bloc3_Barplot_CD8_ProjecTILs.png
#         Results/Tables/Bloc3_CD8_chisq_test.csv
#         Results/Tables/Bloc3_CD8_fisher_posthoc.csv
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(data.table)
})

# Paths
DATA_DIR <- "C:/Users/yasmi/OneDrive/Desktop/Mini-Projets/TumorImmune_Crosstalk_LUAD"
IN_OBJ   <- file.path(DATA_DIR, "Objects/Bloc3_08_seu_CD8_ProjecTILs.rds")
OUT_FIG  <- file.path(DATA_DIR, "Results/Figures/CD8")
OUT_TAB  <- file.path(DATA_DIR, "Results/Tables")

# 1) Load CD8 ProjecTILs object
message("Loading CD8 ProjecTILs object...")
seu_CD8 <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu_CD8))

# 2) Color palette for CD8 states
# Colors chosen to reflect biological meaning:
# TEX/TPEX : warm (red/pink): exhaustion spectrum
# EM/CM : cool (blue/teal): memory/functional spectrum
# Others : neutral
cd8_colors <- c(
  "CD8.TEX"      = "#D73027",  # red — terminal exhaustion
  "CD8.TPEX"     = "#FC8D59",  # orange — precursor exhaustion
  "CD8.EM"       = "#4393C3",  # blue — effector memory
  "CD8.CM"       = "#2166AC",  # dark blue — central memory
  "CD8.NaiveLike"= "#74ADD1",  # light blue — naive-like
  "CD8.TEMRA"    = "#762A83",  # purple — terminally differentiated
  "CD8.MAIT"     = "#A6D96A"   # green — MAIT cells
)

# 3) Compute proportions
message("Computing proportions...")

df <- seu_CD8@meta.data %>%
  filter(!is.na(functional.cluster), !is.na(pathological_response)) %>%
  transmute(response  = pathological_response,
            cd8_state = functional.cluster)

df_prop <- df %>%
  count(response, cd8_state, name = "n") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 4) Barplot
message("Generating barplot...")

p_bar <- ggplot(df_prop,
                aes(x = response, y = prop, fill = cd8_state)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = cd8_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.1)) +
  annotate("text",
           x        = 2,
           y        = 1.08,
           label    = "Chi-2 p < 2.2e-16",
           size     = 3.5,
           fontface = "italic") +
  ylab("CD8 state proportion") +
  ggtitle("CD8 T cell state proportions — GSE243013") +
  theme_classic() +
  theme(axis.title.x   = element_blank(),
        axis.text.x    = element_text(size = 13, face = "bold"),
        axis.text.y    = element_text(size = 11),
        axis.title.y   = element_text(size = 12),
        legend.position = "right",
        legend.text    = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.title   = element_blank(),
        plot.caption   = element_text(size = 8, face = "italic", hjust = 0)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG, "Bloc3_Barplot_CD8_ProjecTILs.png"),
       p_bar, width = 7, height = 6, dpi = 300, bg = "white")
message("Saved: Bloc3_Barplot_CD8_ProjecTILs.png")

# 5) Statistical tests
message("Running statistical tests...")

cont_table <- table(seu_CD8$functional.cluster, seu_CD8$pathological_response)

# Chi-2 global
chisq_res <- chisq.test(cont_table)
message("Chi-2 p-value: ", chisq_res$p.value)

chisq_summary <- data.frame(
  statistic = chisq_res$statistic,
  df        = chisq_res$parameter,
  p.value   = chisq_res$p.value
)
fwrite(chisq_summary, file.path(OUT_TAB, "Bloc3_CD8_chisq_test.csv"))
message("Saved: Bloc3_CD8_chisq_test.csv")

# Fisher post-hoc pairwise
comparisons <- list(
  c("MPR", "non-MPR"),
  c("pCR", "non-MPR"),
  c("MPR", "pCR")
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
fwrite(fisher_df, file.path(OUT_TAB, "Bloc3_CD8_fisher_posthoc.csv"))
message("Saved: Bloc3_CD8_fisher_posthoc.csv")

# 6) Preprint figure — MPR vs non-MPR only
message("Generating preprint barplot MPR vs non-MPR...")

OUT_FIG_PREPRINT <- file.path(DATA_DIR, "Results/Figures/CD8/Preprint")
dir.create(OUT_FIG_PREPRINT, recursive = TRUE, showWarnings = FALSE)

# Chi-2 sur MPR vs non-MPR uniquement
cont_table_2cond <- table(
  seu_CD8$functional.cluster,
  seu_CD8$pathological_response)[, c("MPR", "non-MPR")]
chisq_2cond <- chisq.test(cont_table_2cond)
message("Chi-2 MPR vs non-MPR p-value: ", chisq_2cond$p.value)

df_prop_preprint <- seu_CD8@meta.data %>%
  filter(!is.na(functional.cluster),
         pathological_response %in% c("MPR", "non-MPR")) %>%
  transmute(response  = pathological_response,
            cd8_state = functional.cluster) %>%
  count(response, cd8_state, name = "n") %>%
  group_by(response) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_bar_preprint <- ggplot(df_prop_preprint,
                         aes(x = response, y = prop, fill = cd8_state)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = cd8_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.1)) +
  annotate("text",
           x        = 1.5,
           y        = 1.08,
           label    = paste0("Chi-2 p = ", signif(chisq_2cond$p.value, 3)),
           size     = 3.5,
           fontface = "italic") +
  ggtitle("CD8 T cell state proportions — GSE243013") +
  ylab("CD8 state proportion") +
  theme_classic() +
  theme(axis.title.x    = element_blank(),
        axis.text.x     = element_text(size = 13, face = "bold"),
        axis.text.y     = element_text(size = 11),
        axis.title.y    = element_text(size = 12),
        legend.position = "right",
        legend.text     = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.title    = element_blank()) +
  guides(fill = guide_legend(ncol = 1))

ggsave(file.path(OUT_FIG_PREPRINT, "Bloc3_Barplot_CD8_ProjecTILs_MPR_nonMPR.png"),
       p_bar_preprint, width = 7, height = 6, dpi = 300, bg = "white")
message("Saved: Bloc3_Barplot_CD8_ProjecTILs_MPR_nonMPR.png — Preprint")

message("DONE Bloc3 Script 09")