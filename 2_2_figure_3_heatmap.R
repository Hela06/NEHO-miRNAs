# =============================================================================
# Figure 3 — Heatmap: miRNA associations with pollutants
# -----------------------------------------------------------------------------
# Description:
#   Combines results from linear regression and BWQS into a single heatmap
#   showing effect sizes (estimates / beta1) for each miRNA-pollutant pair.
#   Rows = pollutants and mixture scores; columns = miRNAs (maternal + cordonal).
#   Only FDR-significant pairs are shown for linear regression;
#   only credible-interval-significant pairs for BWQS (mean_flagged).
#
# Input:
#   - regression_MM.xlsx                          (output of 2_0_linearRegression_analysis.R)
#   - regression_CC.xlsx                          (output of 2_0_linearRegression_analysis.R)
#   - output_bWQS/WQS_all_results.xlsx            (output of 2_1_BWQS.R)
#
# Output:
#   - Figure3_heatmap.tiff
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(data.table)
library(plyr)


# --- Load data ---------------------------------------------------------------
df_MM   <- as.data.frame(readxl::read_excel("./regression_MM.xlsx"))
df_CC   <- as.data.frame(readxl::read_excel("./regression_CC.xlsx"))
WQS_all <- as.data.frame(readxl::read_excel("./output_bWQS/WQS_all_results.xlsx"))


# =============================================================================
# SECTION 1: BWQS — wide matrix (rows = mixture labels, cols = miRNA)
# =============================================================================

input_mat_WQS <- WQS_all %>%
  select(MET_LIP, mean_flagged, miRNA_fin) %>%
  pivot_wider(names_from = miRNA_fin, values_from = mean_flagged) %>%
  as.data.frame()

row.names(input_mat_WQS) <- input_mat_WQS$MET_LIP
input_mat_WQS$MET_LIP    <- NULL


# =============================================================================
# SECTION 2: LINEAR REGRESSION — wide matrix (FDR < 0.05 only)
# =============================================================================

# Combine MM and CC results, keep only significant pairs
df_reg <- rbind(df_MM, df_CC)

df_reg_sig <- df_reg[which(df_reg$p.adjust.FDR < 0.05), c("pollutant", "estimate", "miRNA")]
colnames(df_reg_sig)[1] <- "MET_LIP"
df_reg_sig$miRNA_fin    <- gsub("_CD$", "", df_reg_sig$miRNA)

input_mat_reg <- df_reg_sig %>%
  select(MET_LIP, estimate, miRNA_fin) %>%
  pivot_wider(names_from = miRNA_fin, values_from = estimate) %>%
  as.data.frame()

row.names(input_mat_reg) <- input_mat_reg$MET_LIP
input_mat_reg$MET_LIP    <- NULL

# Remove derived/low-signal pollutant rows not shown in figure
rows_to_drop <- c("RAT_Cu_Zn", "M_PCB156", "M_PCB170", "M_PCB183",
                  "M_PCB187", "SUM_PCB", "RAT_C_Cu_Zn")
input_mat_reg <- input_mat_reg[!row.names(input_mat_reg) %in% rows_to_drop, ]

# Order rows by biological grouping (maternal first, then cordonal)
row_order_reg <- c(
  "M_Cu", "M_Se", "M_Zn", "M_DDE", "M_HCB",
  "M_PCB74", "M_PCB118", "M_PCB138", "M_PCB153", "M_PCB180", "SUM_PCBind",
  "C_Cu", "C_Se", "C_Zn", "C_HCB",
  "C_PCB138", "C_PCB153", "C_PCB180", "SUM_C_PCBind"
)
input_mat_reg <- input_mat_reg[
  order(factor(row.names(input_mat_reg), levels = row_order_reg)), ,
  drop = FALSE
]


# =============================================================================
# SECTION 3: COMBINE REGRESSION + BWQS INTO ONE MATRIX
# =============================================================================

input_mat_combined <- plyr::rbind.fill(input_mat_reg, input_mat_WQS)
row.names(input_mat_combined) <- c(row.names(input_mat_reg), row.names(input_mat_WQS))

# Remove columns that are entirely NA (miRNAs with no associations)
input_mat_combined <- Filter(function(x) !all(is.na(x)), input_mat_combined)

# Add empty separator rows between maternal and cordonal blocks
sep_rows           <- as.data.frame(matrix(NA, nrow = 4, ncol = ncol(input_mat_combined)))
colnames(sep_rows) <- colnames(input_mat_combined)
row.names(sep_rows) <- c("M_Cu", "M_Se", "C_DDE", "C_HCB")
input_mat_combined  <- rbind(input_mat_combined, sep_rows)

# Final row order (regression rows + WQS rows)
row_order_final <- c(
  "M_Cu", "M_Se", "M_Zn", "M_DDE", "M_HCB",
  "M_PCB74", "M_PCB118", "M_PCB138", "M_PCB153", "M_PCB180", "SUM_PCBind",
  "M_WQS_MIX_EEs", "M_WQS_MIX_OCs", "M_WQS_MIX_TOT",
  "C_Cu", "C_Se", "C_Zn", "C_DDE", "C_HCB",
  "C_PCB138", "C_PCB153", "C_PCB180", "SUM_C_PCBind",
  "C_WQS_MIX_EEs", "C_WQS_MIX_OCs", "C_WQS_MIX_TOT"
)
input_mat_combined <- input_mat_combined[
  order(factor(row.names(input_mat_combined), levels = row_order_final)), ,
  drop = FALSE
]


# =============================================================================
# SECTION 4: ROW ANNOTATIONS
# =============================================================================

annot_row       <- data.frame(row.names = row.names(input_mat_combined))
annot_row$SIDE  <- ifelse(grepl("^M_|^SUM_P", row.names(annot_row)), "MATERNAL", "CORD")
annot_row$SIDE[row.names(annot_row) == "SUM_C_PCBind"] <- "CORD"

# Clean display labels
annot_row$NAME <- gsub("^M_|^C_", "", row.names(annot_row))
annot_row$NAME[row.names(annot_row) %in% c("SUM_PCBind", "SUM_C_PCBind")]        <- "PCBs-SUM"
annot_row$NAME[row.names(annot_row) %in% c("M_WQS_MIX_EEs", "C_WQS_MIX_EEs")]  <- "EEs"
annot_row$NAME[row.names(annot_row) %in% c("M_WQS_MIX_OCs", "C_WQS_MIX_OCs")]  <- "POPs"
annot_row$NAME[row.names(annot_row) %in% c("M_WQS_MIX_TOT", "C_WQS_MIX_TOT")]  <- "TOTAL"

row_labels     <- annot_row$NAME
annot_row$NAME <- NULL


# =============================================================================
# SECTION 5: PLOT AND SAVE
# =============================================================================

ph <- pheatmap(
  as.matrix(input_mat_combined),
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  main          = "Associations miRNA vs Pollutant",
  legend        = TRUE,
  fontsize_row  = 17,
  fontsize_col  = 17,
  angle_col     = 45,
  cellwidth     = 30,
  cellheight    = 30,
  border_color  = "gray",
  fontsize      = 15,
  annotation_row = annot_row,
  labels_row    = row_labels,
  gaps_row      = c(11, 14, 23),
  gaps_col      = 22
)

tiff("./Figure3_heatmap.tiff", res = 300, width = 27, height = 15, units = "in")
print(ph)
dev.off()

cat("Figure saved to Figure3_heatmap.tiff\n")
