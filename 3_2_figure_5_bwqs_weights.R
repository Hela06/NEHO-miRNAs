# =============================================================================
# Figure 5 — BWQS Mixture Weights Plot
# -----------------------------------------------------------------------------
# Description:
#   For each analysis (EEs, OCs, TOT) and each side (MATERNAL, CORD):
#     - Extracts BWQS weights (mix_name rows) for significant miRNAs
#     - Computes summary statistics per chemical (mean, median, IQR)
#     - Computes support ratio: proportion of miRNAs with weight > 1/k threshold
#     - Plots weight distributions as boxplot + jitter + mean diamond
#   Maternal panels use scale_y_continuous (left-to-right);
#   cordonal panels use scale_y_continuous (mirrored layout via ggarrange).
#   Plots are combined into a 3x3 grid (rows = EEs/POPs/TOT, cols = M/label/C).
#
# Input:
#   - output_bWQS/FINAL_MODELS_significant_TOT.xlsx     (output of 2_1_BWQS.R)
#   - output_bWQS/FINAL_MODELS_significant_divided.xlsx (output of 2_1_BWQS.R)
#   - output_bWQS/FINAL_MODELS_nonsignificant.xlsx      (output of 2_1_BWQS.R)
#   - output_bWQS/WQS_all_results.xlsx                  (output of 2_1_BWQS.R)
#
# Output:
#   - Figure5_bwqs_weights.tiff
#   - Figure5_bwqs_weights_input.xlsx  (one sheet per panel)
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
library(ggpubr)


# --- Load data ---------------------------------------------------------------
bwqs_TOT  <- as.data.frame(readxl::read_excel("./output_bWQS/FINAL_MODELS_significant_TOT.xlsx"))
bwqs_div  <- as.data.frame(readxl::read_excel("./output_bWQS/FINAL_MODELS_significant_divided.xlsx"))
bwqs_insg <- as.data.frame(readxl::read_excel("./output_bWQS/FINAL_MODELS_nonsignificant.xlsx"))
WQS_all   <- as.data.frame(readxl::read_excel("./output_bWQS/WQS_all_results.xlsx"))

bwqs_TOT$ANALYSIS <- "TOT"


# =============================================================================
# SECTION 1: ASSEMBLE WEIGHT DATASET
# =============================================================================

# Combine all BWQS results and keep only weight rows (not beta1/intercept)
WGT_all <- bind_rows(bwqs_TOT, bwqs_div, bwqs_insg)

# Pollutant columns are identified by mix_name starting with M_ or C_
WGT_all <- WGT_all[grepl("^M_|^C_", WGT_all$mix_name), ]

# Keep only miRNAs that had a significant beta1 (credible interval excludes zero)
sig_mirna <- unique(WQS_all$miRNA[!is.na(WQS_all$mean_flagged)])
WGT_plot  <- WGT_all[WGT_all$miRNA %in% sig_mirna, ]

# Add side and clean chemical name
WGT_plot$SIDE     <- ifelse(grepl("_CD$", WGT_plot$miRNA), "CHILD", "MOTHER")
WGT_plot$met_name <- gsub("^M_|^C_", "", WGT_plot$mix_name)
WGT_plot$met_name <- factor(WGT_plot$met_name,
                             levels = c("Cu", "Se", "Zn", "DDE", "HCB",
                                        "PCB74", "PCB118", "PCB138", "PCB153", "PCB180"))


# =============================================================================
# SECTION 2: HELPER FUNCTIONS
# =============================================================================

# Compute per-chemical summary and support ratio
make_summary <- function(df, threshold) {
  df %>%
    mutate(mean_perc = mean * 100) %>%
    group_by(met_name) %>%
    summarise(
      MeanWeight   = mean(mean, na.rm = TRUE),
      MedianWeight = median(mean, na.rm = TRUE),
      P25          = quantile(mean, 0.25, na.rm = TRUE),
      P75          = quantile(mean, 0.75, na.rm = TRUE),
      P10          = quantile(mean, 0.10, na.rm = TRUE),
      P90          = quantile(mean, 0.90, na.rm = TRUE),
      n_miRNAs     = sum(mean * 100 > threshold, na.rm = TRUE),
      total_miRNAs = n(),
      support_ratio = n_miRNAs / total_miRNAs,
      .groups = "drop"
    )
}

# Build one weight plot panel
make_weight_plot <- function(dfj, summary_df, color, threshold, y_reverse = FALSE) {
  p <- ggplot(dfj, aes(x = met_name, y = mean_perc)) +
    geom_col(data = summary_df,
             aes(x = met_name, y = support_ratio * 100),
             fill = color, alpha = 0.1, width = 0.9) +
    geom_jitter(shape = 21, color = "black", size = 3,
                fill = color, alpha = 0.2, width = 0.25) +
    geom_boxplot(aes(group = met_name), outlier.shape = NA,
                 alpha = 0.3, fill = color, width = 0.8) +
    geom_point(aes(y = MeanWeight * 100),
               color = ifelse(color == "maroon2", "violetred", "royalblue3"),
               size = 10, shape = 18) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
    theme_minimal() +
    xlab(" ") +
    theme(
      axis.text.x      = element_text(size = 20, face = "bold"),
      panel.grid.major.x = element_blank(),
      axis.line        = element_line(linewidth = 1, colour = "gray")
    ) +
    coord_flip()

  if (y_reverse) {
    p <- p + scale_y_continuous(name = " ", limits = c(0, 100),
                                sec.axis = sec_axis(~ . / 100)) +
             scale_x_discrete(position = "top")
  } else {
    p <- p + scale_y_continuous(name = " ", limits = c(0, 100),
                                sec.axis = sec_axis(~ . / 100))
  }
  p
}


# =============================================================================
# SECTION 3: BUILD PANELS
# =============================================================================

# --- EEs maternal (threshold = 1/3) -----------------------------------------
t1    <- (1 / 3) * 100
EE_M  <- WGT_plot[WGT_plot$ANALYSIS == "EEs" & WGT_plot$SIDE == "MOTHER", ]
EE_Mj <- left_join(EE_M, make_summary(EE_M, t1), by = "met_name") %>%
           mutate(mean_perc = mean * 100)
chemical_supp_EE_M <- make_summary(EE_M, t1)

EE_M_plot <- make_weight_plot(EE_Mj, chemical_supp_EE_M, "maroon2", t1, y_reverse = TRUE)


# --- EEs cordonal ------------------------------------------------------------
EE_C  <- WGT_plot[WGT_plot$ANALYSIS == "EEs" & WGT_plot$SIDE == "CHILD", ]
EE_Cj <- left_join(EE_C, make_summary(EE_C, t1), by = "met_name") %>%
           mutate(mean_perc = mean * 100)
chemical_supp_EE_C <- make_summary(EE_C, t1)

EE_C_plot <- make_weight_plot(EE_Cj, chemical_supp_EE_C, "royalblue1", t1)


# --- POPs maternal (threshold = 1/7) ----------------------------------------
t2     <- (1 / 7) * 100
POP_M  <- WGT_plot[WGT_plot$ANALYSIS == "OCs" & WGT_plot$SIDE == "MOTHER", ]
POP_Mj <- left_join(POP_M, make_summary(POP_M, t2), by = "met_name") %>%
            mutate(mean_perc = mean * 100)
chemical_supp_POP_M <- make_summary(POP_M, t2)

POP_M_plot <- make_weight_plot(POP_Mj, chemical_supp_POP_M, "maroon2", t2, y_reverse = TRUE)


# --- POPs cordonal (threshold = 1/5, add NA rows for PCB74/PCB118) ----------
t2a    <- (1 / 5) * 100
POP_C  <- WGT_plot[WGT_plot$ANALYSIS == "OCs" & WGT_plot$SIDE == "CHILD", ]

# PCB74 and PCB118 not present in cordonal mixture — add as NA placeholders
df_miss <- POP_C[1:2, ]
df_miss[, setdiff(colnames(df_miss), c("mix_name", "met_name"))] <- NA
df_miss$mix_name <- c("C_PCB74", "C_PCB118")
df_miss$met_name <- factor(c("PCB74", "PCB118"), levels = levels(POP_C$met_name))
POP_C  <- bind_rows(POP_C, df_miss)

POP_Cj <- left_join(POP_C, make_summary(POP_C, t2a), by = "met_name") %>%
            mutate(mean_perc = mean * 100)
chemical_supp_POP_C <- make_summary(POP_C, t2a)

POP_C_plot <- make_weight_plot(POP_Cj, chemical_supp_POP_C, "royalblue1", t2a)


# --- TOT maternal (threshold = 1/10) ----------------------------------------
t3     <- (1 / 10) * 100
TOT_M  <- WGT_plot[WGT_plot$ANALYSIS == "TOT" & WGT_plot$SIDE == "MOTHER", ]
TOT_Mj <- left_join(TOT_M, make_summary(TOT_M, t3), by = "met_name") %>%
            mutate(mean_perc = mean * 100)
chemical_supp_TOT_M <- make_summary(TOT_M, t3)

TOT_M_plot <- make_weight_plot(TOT_Mj, chemical_supp_TOT_M, "maroon2", t3, y_reverse = TRUE)


# --- TOT cordonal (threshold = 1/8, reuse df_miss) --------------------------
t3a    <- (1 / 8) * 100
TOT_C  <- WGT_plot[WGT_plot$ANALYSIS == "TOT" & WGT_plot$SIDE == "CHILD", ]
TOT_C  <- bind_rows(TOT_C, df_miss)

TOT_Cj <- left_join(TOT_C, make_summary(TOT_C, t3a), by = "met_name") %>%
            mutate(mean_perc = mean * 100)
chemical_supp_TOT_C <- make_summary(TOT_C, t3a)

TOT_C_plot <- make_weight_plot(TOT_Cj, chemical_supp_TOT_C, "royalblue1", t3a)


# =============================================================================
# SECTION 4: EXPORT INPUT DATA (one sheet per panel)
# =============================================================================

writexl::write_xlsx(
  list(
    EE_Mj  = as.data.frame(EE_Mj),
    EE_Cj  = as.data.frame(EE_Cj),
    POP_Mj = as.data.frame(POP_Mj),
    POP_Cj = as.data.frame(POP_Cj),
    TOT_Mj = as.data.frame(TOT_Mj),
    TOT_Cj = as.data.frame(TOT_Cj)
  ),
  "Figure5_bwqs_weights_input.xlsx"
)

cat("Input data saved to Figure5_bwqs_weights_input.xlsx\n")


# =============================================================================
# SECTION 5: COMBINE AND SAVE
# =============================================================================

y_label_EE  <- text_grob(c("Zn", "Se", "Cu"),
                          color = "black", face = "bold", size = 20,
                          y = c(0.7, 0.5, 0.3))
y_label_POP <- text_grob(c("PCB180","PCB153","PCB138","PCB118","PCB74","HCB","DDE"),
                          color = "black", face = "bold", size = 20,
                          y = c(0.88, 0.76, 0.64, 0.52, 0.40, 0.28, 0.16))
y_label_TOT <- text_grob(c("PCB180","PCB153","PCB138","PCB118","PCB74","HCB","DDE","Zn","Se","Cu"),
                          color = "black", face = "bold", size = 20,
                          y = c(0.9, 0.82, 0.73, 0.64, 0.55, 0.47, 0.39, 0.30, 0.21, 0.12))

tot_plt <- ggarrange(
  EE_M_plot  + rremove("y.text"), y_label_EE,  EE_C_plot  + rremove("y.text"),
  POP_M_plot + rremove("y.text"), y_label_POP, POP_C_plot + rremove("y.text"),
  TOT_M_plot + rremove("y.text"), y_label_TOT, TOT_C_plot + rremove("y.text"),
  ncol = 3, nrow = 3, widths = c(1, 0.1, 1), heights = c(10, 20, 27)
)

tot_plt_annt <- annotate_figure(
  tot_plt,
  top    = text_grob(c("MATERNAL", "WEIGHT", "CORD"),
                     color = "black", face = c("bold","plain","bold"),
                     size = c(25, 20, 25), x = c(0.2, 0.5, 0.8)),
  bottom = text_grob("%", color = "black", size = 20),
  left   = text_grob("CHEMICALS", color = "black", rot = 90, face = "bold", size = 20),
  right  = text_grob(c("EEs", "POPs", "TOTAL"), color = "black",
                     rot = 90, face = "bold", size = 20, y = c(0.9, 0.65, 0.25))
)

tiff("Figure5_bwqs_weights.tiff", res = 300, width = 27, height = 15, units = "in")
print(tot_plt_annt)
dev.off()

cat("Figure saved to Figure5_bwqs_weights.tiff\n")
