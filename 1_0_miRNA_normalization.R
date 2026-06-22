# =============================================================================
# miRNA Preprocessing and Normalization
# -----------------------------------------------------------------------------
# Description:
#   This script performs:
#   1. Loading of calibrated miRNA Ct values (maternal and cordonal serum)
#   2. Filtering of miRNAs with >= 30% undetermined values (Ct = 40)
#   3. Gene normalization (delta-Ct using housekeeping genes)
#   4. Loading and LOQ-based filtering of pollutant data (metals and POPs)
#   5. Lipid adjustment of pollutant concentrations
#   6. Placental transfer ratio computation (cordonal / maternal)
#   7. Clinical data integration
#   8. Final merged dataset export
#
# Input files:
#   - data/example/synthetic_miRNA.xlsx       (calibrated Ct values, maternal and cordonal serum)
#   - data/example/synthetic_poll.xlsx     (chemical concentrations, maternal and cordonal)
#   - data/example/synthetic_lipids.xlsx   (total lipids, maternal and cordonal)
#   - data/example/synthetic_clinical.xlsx (clinical and demographic data)
#
# NOTE: Original patient data cannot be shared due to ethical restrictions.
#       Synthetic datasets with identical structure are provided for reproducibility.
#       For data access requests, contact:
#         - gaspare.drago@irib.cnr.it
#         - silvia.ruggieri@irib.cnr.it
#
# Output:
#   - input_dataset.xlsx   (merged, normalized, ready-for-analysis dataset)
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(writexl)
library(data.table)
library(dplyr)
library(tidyr)


# --- File paths --------------------------------------------------------------
# Replace with your actual data paths if you have access to the original data
mirna.input        <- "./data/example/synthetic_miRNA.xlsx"
pollutants.input   <- "./data/example/synthetic_poll.xlsx"
lipidomics.input   <- "./data/example/synthetic_lipids.xlsx"
clinical.input     <- "./data/example/synthetic_clinical.xlsx"


# =============================================================================
# SECTION 1: miRNA PREPROCESSING
# =============================================================================

ds <- readxl::read_xlsx(mirna.input)
ds[, 2:177] <- as.data.frame(lapply(ds[, 2:177], as.numeric))


# --- 1.1 Filter miRNAs with >= 30% undetermined values (Ct = 40) -------------
cnt_mR <- data.frame(MIRNA = colnames(ds)[2:177], COUNT = NA, PERC = NA)
for (r in colnames(ds)[2:177]) {
  cnt_mR$COUNT[cnt_mR$MIRNA == r] <- sum(ds[, r] >= 40, na.rm = TRUE)
  cnt_mR$PERC[cnt_mR$MIRNA == r]  <- (cnt_mR$COUNT[cnt_mR$MIRNA == r] * 100) / nrow(ds)
}
miR_rm <- cnt_mR$MIRNA[cnt_mR$PERC >= 30]


# --- 1.2 Gene normalization (delta-Ct) ---------------------------------------
# Columns 85-88  = maternal housekeeping genes
# Columns 173-176 = cordonal housekeeping genes
ds_st_genN <- ds
row.names(ds_st_genN) <- ds_st_genN$Row.names
ds_st_genN$Row.names  <- NULL

for (rr in 1:nrow(ds_st_genN)) {
  hk_maternal <- as.numeric(ds_st_genN[rr, 85:88])
  hk_cordonal <- as.numeric(ds_st_genN[rr, 173:176])

  # Subtract mean Ct of housekeeping genes (delta-Ct normalization)
  for (cc1 in 1:84)    ds_st_genN[rr, cc1] <- ds_st_genN[rr, cc1] - mean(hk_maternal)
  for (cc2 in 89:172)  ds_st_genN[rr, cc2] <- ds_st_genN[rr, cc2] - mean(hk_cordonal)
}

# Remove housekeeping gene columns and low-quality miRNAs
ds_st_genN[, c(85:88, 173:176)] <- NULL
ds_st_genN <- ds_st_genN[, !colnames(ds_st_genN) %in% miR_rm]


# =============================================================================
# SECTION 2: POLLUTANT PREPROCESSING
# =============================================================================

# --- 2.1 Load pollutant data -------------------------------------------------
poll <- as.data.frame(readxl::read_xlsx(pollutants.input))
poll[, 2:ncol(poll)] <- lapply(poll[, 2:ncol(poll)], as.numeric)


# --- 2.2 LOQ-based filtering -------------------------------------------------
# LOQ replacement values: 2.5 for most analytes, 0.2 for Hg, 0.5 for As
pollutant_cols <- colnames(poll)[2:ncol(poll)]

cnt <- data.frame(INQUINANTI = pollutant_cols, COUNT = NA, PERC = NA)
for (r in pollutant_cols) {
  cnt$COUNT[cnt$INQUINANTI == r] <- sum(poll[, r] == 2.5, na.rm = TRUE)
  cnt$PERC[cnt$INQUINANTI == r]  <- (cnt$COUNT[cnt$INQUINANTI == r] * 100) / nrow(poll)
}
# Special LOQ thresholds for Hg and As
for (analyte in c("M_Hg", "C_Hg")) {
  cnt$COUNT[cnt$INQUINANTI == analyte] <- sum(poll[, analyte] == 0.2, na.rm = TRUE)
  cnt$PERC[cnt$INQUINANTI == analyte]  <- (cnt$COUNT[cnt$INQUINANTI == analyte] * 100) / nrow(poll)
}
for (analyte in c("M_As", "C_As")) {
  cnt$COUNT[cnt$INQUINANTI == analyte] <- sum(poll[, analyte] == 0.5, na.rm = TRUE)
  cnt$PERC[cnt$INQUINANTI == analyte]  <- (cnt$COUNT[cnt$INQUINANTI == analyte] * 100) / nrow(poll)
}

pll_rm <- cnt$INQUINANTI[cnt$PERC >= 30]
poll[, colnames(poll) %in% pll_rm] <- NULL

row.names(poll) <- poll$shortCode
poll$shortCode  <- NULL
pll_st          <- na.omit(poll)


# --- 2.3 Lipid adjustment of pollutant concentrations ------------------------
lipids <- as.data.frame(readxl::read_xlsx(lipidomics.input))
pll_st <- merge(pll_st, lipids, by.x = 0, by.y = "shortCode", all.x = TRUE)
row.names(pll_st) <- pll_st$Row.names
pll_st$Row.names  <- NULL

# Identify maternal (M_) and cordonal (C_) pollutant columns, excluding TL
m_poll_cols <- setdiff(grep("^M_", colnames(pll_st), value = TRUE), "M_TL_1")
c_poll_cols <- setdiff(grep("^C_", colnames(pll_st), value = TRUE), "C_TL_1")

pll_st[, m_poll_cols] <- lapply(pll_st[, m_poll_cols], function(x) x / pll_st$M_TL_1)
pll_st[, c_poll_cols] <- lapply(pll_st[, c_poll_cols], function(x) x / pll_st$C_TL_1)
pll_st[, c("M_TL_1", "C_TL_1")] <- NULL


# --- 2.4 Derived pollutant variables -----------------------------------------
pcb_cols_m <- intersect(c("M_PCB74","M_PCB118","M_PCB138","M_PCB153",
                           "M_PCB156","M_PCB170","M_PCB180","M_PCB187"),
                        colnames(pll_st))
if (length(pcb_cols_m) > 0)
  pll_st$SUM_PCB <- rowSums(pll_st[, pcb_cols_m], na.rm = TRUE)


# --- 2.5 Placental transfer ratios (cordonal / maternal) ---------------------
# Computed for pollutants present in both maternal and cordonal serum
shared_analytes <- gsub("^C_", "", intersect(
  gsub("^M_", "", m_poll_cols),
  gsub("^C_", "", c_poll_cols)
))

for (analyte in shared_analytes) {
  m_col <- paste0("M_", analyte)
  c_col <- paste0("C_", analyte)
  if (m_col %in% colnames(pll_st) && c_col %in% colnames(pll_st)) {
    pll_st[[paste0("PLAC_TR_", analyte)]] <- pll_st[, c_col] / pll_st[, m_col]
  }
}


# =============================================================================
# SECTION 3: CLINICAL DATA INTEGRATION
# =============================================================================

info_tot <- as.data.frame(readxl::read_xlsx(clinical.input))


# =============================================================================
# SECTION 4: FINAL DATASET ASSEMBLY
# =============================================================================

# Z-score standardization of miRNA values
mir_scaled <- as.data.frame(scale(ds_st_genN))

# Merge miRNA + pollutants + clinical
DS_TOT <- merge(mir_scaled, pll_st,   by = 0, all.x = TRUE)
row.names(DS_TOT) <- DS_TOT$Row.names
DS_TOT$Row.names  <- NULL

DS_TOT <- merge(DS_TOT, info_tot, by.x = 0, by.y = "shortCode", all.x = TRUE)
row.names(DS_TOT) <- DS_TOT$Row.names
DS_TOT$Row.names  <- NULL

writexl::write_xlsx(tibble::rownames_to_column(DS_TOT, "shortCode"), "./input_dataset.xlsx")

cat("Done. Final dataset saved to input_dataset.xlsx\n")
cat(sprintf("  Subjects : %d\n", nrow(DS_TOT)))
cat(sprintf("  Variables: %d\n", ncol(DS_TOT)))
