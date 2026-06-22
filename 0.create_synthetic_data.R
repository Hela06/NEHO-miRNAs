# =============================================================================
# Synthetic Dataset Generator
# -----------------------------------------------------------------------------
# Generates synthetic versions of the three datasets used in the analysis:
#   1. synthetic_miRNA.xlsx       — miRNA Ct values (maternal + cordonal)
#   2. synthetic_poll.xlsx     — pollutant concentrations (maternal + cordonal)
#   3. synthetic_lipids.xlsx   — total lipids (maternal + cordonal)
#   4. synthetic_clinical.xlsx
#
# All values are randomly generated and do NOT represent real patients.
# Row identifier (shortCode) is consistent across the four files.
#
# Usage: source this script once to populate data/example/
# =============================================================================

set.seed(42)
library(writexl)

n <- 35

# --- Subject codes -----------------------------------------------------------
subject_ids <- c(
  paste0("SIR", sprintf("%02d", 1:28)),
  paste0("LEN", sprintf("%03d", 101:107))
)


# =============================================================================
# 1. miRNA DATASET
# =============================================================================

mirna_maternal <- c(
  "hsa-let-7a-5p", "hsa-miR-1-3p", "hsa-miR-100-5p", "hsa-miR-106b-5p",
  "hsa-miR-10b-5p", "hsa-miR-122-5p", "hsa-miR-124-3p", "hsa-miR-125b-5p",
  "hsa-miR-126-3p", "hsa-miR-133a-3p", "hsa-miR-133b", "hsa-miR-134-5p",
  "hsa-miR-141-3p", "hsa-miR-143-3p", "hsa-miR-146a-5p", "hsa-miR-150-5p",
  "hsa-miR-155-5p", "hsa-miR-17-5p", "hsa-miR-17-3p", "hsa-miR-18a-5p",
  "hsa-miR-192-5p", "hsa-miR-195-5p", "hsa-miR-196a-5p", "hsa-miR-19a-3p",
  "hsa-miR-19b-3p", "hsa-miR-200a-3p", "hsa-miR-200b-3p", "hsa-miR-200c-3p",
  "hsa-miR-203a-3p", "hsa-miR-205-5p", "hsa-miR-208a-3p", "hsa-miR-20a-5p",
  "hsa-miR-21-5p", "hsa-miR-210-3p", "hsa-miR-214-3p", "hsa-miR-215-5p",
  "hsa-miR-221-3p", "hsa-miR-222-3p", "hsa-miR-223-3p", "hsa-miR-224-5p",
  "hsa-miR-23a-3p", "hsa-miR-25-3p", "hsa-miR-27a-3p", "hsa-miR-296-5p",
  "hsa-miR-29a-3p", "hsa-miR-30d-5p", "hsa-miR-34a-5p", "hsa-miR-375-3p",
  "hsa-miR-423-5p", "hsa-miR-499a-5p", "hsa-miR-574-3p", "hsa-miR-885-5p",
  "hsa-miR-9-5p", "hsa-miR-92a-3p", "hsa-miR-93-5p", "hsa-let-7c-5p",
  "hsa-miR-107", "hsa-miR-10a-5p", "hsa-miR-128-3p", "hsa-miR-130b-3p",
  "hsa-miR-145-5p", "hsa-miR-148a-3p", "hsa-miR-15a-5p", "hsa-miR-184",
  "hsa-miR-193a-5p", "hsa-miR-204-5p", "hsa-miR-206", "hsa-miR-211-5p",
  "hsa-miR-26b-5p", "hsa-miR-30e-5p", "hsa-miR-372-3p", "hsa-miR-373-3p",
  "hsa-miR-374a-5p", "hsa-miR-376c-3p", "hsa-miR-7-5p", "hsa-miR-96-5p",
  "hsa-miR-103a-3p", "hsa-miR-15b-5p", "hsa-miR-16-5p", "hsa-miR-191-5p",
  "hsa-miR-22-3p", "hsa-miR-24-3p", "hsa-miR-26a-5p", "hsa-miR-31-5p",
  "hsa-miR-30c-5p", "hsa-miR-103a-3p", "hsa-miR-451a", "hsa-miR-23a-3p"
)
mirna_cordonal <- paste0(mirna_maternal, "_CD")

gen_ct <- function(n_subjects, n_mirna, pct_undet = 0.08) {
  mat <- matrix(round(rnorm(n_subjects * n_mirna, mean = 29, sd = 4), 2),
                nrow = n_subjects)
  mat <- pmin(pmax(mat, 17), 39.99)
  mat[sample(length(mat), floor(length(mat) * pct_undet))] <- 40
  as.data.frame(mat)
}

ct_mat <- gen_ct(n, length(mirna_maternal)); colnames(ct_mat) <- mirna_maternal
ct_cd  <- gen_ct(n, length(mirna_cordonal)); colnames(ct_cd)  <- mirna_cordonal

ds_synthetic <- cbind(data.frame(Row.names = subject_ids), ct_mat, ct_cd)


# =============================================================================
# 2. POLLUTANT DATASET
# =============================================================================
# Values generated within the observed range from summary statistics.
# No NAs included (synthetic data for code testing only).

# Helper: truncated normal within [lo, hi]
rtrunc <- function(n, mean, sd, lo, hi) {
  x <- rnorm(n, mean, sd)
  x <- pmax(pmin(x, hi), lo)
  round(x, 3)
}

poll_synthetic <- data.frame(
  shortCode = subject_ids,

  # --- Maternal metals -------------------------------------------------------
  M_Se  = rtrunc(n, mean = 79,   sd = 14,   lo = 51,    hi = 113),
  M_Hg  = rtrunc(n, mean = 0.84, sd = 0.80, lo = 0.2,   hi = 4.56),
  M_As  = rtrunc(n, mean = 1.52, sd = 3.0,  lo = 0.5,   hi = 25.3),
  M_Zn  = rtrunc(n, mean = 647,  sd = 120,  lo = 395,   hi = 1064),
  M_Cu  = rtrunc(n, mean = 2145, sd = 450,  lo = 967,   hi = 3431),

  # --- Maternal POPs ---------------------------------------------------------
  M_HCB  = rtrunc(n, mean = 56.7,  sd = 28,   lo = 15.9,  hi = 149.9),
  M_TNC  = rtrunc(n, mean = 6.94,  sd = 7.5,  lo = 2.5,   hi = 58.1),
  M_DDE  = rtrunc(n, mean = 424.8, sd = 450,  lo = 50.4,  hi = 3853),
  M_PCB74  = rtrunc(n, mean = 11.2,  sd = 7,    lo = 2.5,   hi = 38.5),
  M_PCB118 = rtrunc(n, mean = 22.8,  sd = 15,   lo = 2.5,   hi = 93.9),
  M_PCB138 = rtrunc(n, mean = 67.4,  sd = 45,   lo = 11.2,  hi = 285),
  M_PCB153 = rtrunc(n, mean = 118.5, sd = 80,   lo = 18.2,  hi = 519),
  M_PCB156 = rtrunc(n, mean = 10.6,  sd = 6,    lo = 2.5,   hi = 40),
  M_PCB170 = rtrunc(n, mean = 37.3,  sd = 28,   lo = 2.5,   hi = 162),
  M_PCB180 = rtrunc(n, mean = 87.3,  sd = 65,   lo = 8.9,   hi = 400),
  M_PCB183 = rtrunc(n, mean = 8.71,  sd = 8,    lo = 2.5,   hi = 53.6),
  M_PCB187 = rtrunc(n, mean = 27.2,  sd = 25,   lo = 2.5,   hi = 198.6),

  # --- Cordonal metals -------------------------------------------------------
  C_Se  = rtrunc(n, mean = 51.3, sd = 12,   lo = 33,    hi = 94),
  C_Hg  = rtrunc(n, mean = 0.63, sd = 0.55, lo = 0.2,   hi = 3.1),
  C_As  = rtrunc(n, mean = 1.41, sd = 2.8,  lo = 0.5,   hi = 20.1),
  C_Zn  = rtrunc(n, mean = 876,  sd = 140,  lo = 560,   hi = 1188),
  C_Cu  = rtrunc(n, mean = 411,  sd = 220,  lo = 100,   hi = 2234),

  # --- Cordonal POPs ---------------------------------------------------------
  C_HCB  = rtrunc(n, mean = 24.4,  sd = 18,   lo = 5.0,   hi = 132.7),
  C_TNC  = rtrunc(n, mean = 2.59,  sd = 1.2,  lo = 2.5,   hi = 10.6),
  C_DDE  = rtrunc(n, mean = 107.0, sd = 140,  lo = 20.0,  hi = 1046),
  C_PCB74  = rtrunc(n, mean = 2.83,  sd = 2.0,  lo = 2.5,   hi = 13.4),
  C_PCB118 = rtrunc(n, mean = 4.45,  sd = 3.5,  lo = 2.5,   hi = 21.2),
  C_PCB138 = rtrunc(n, mean = 12.7,  sd = 10,   lo = 2.5,   hi = 66.3),
  C_PCB153 = rtrunc(n, mean = 21.5,  sd = 18,   lo = 2.5,   hi = 118.2),
  C_PCB156 = rtrunc(n, mean = 2.65,  sd = 1.2,  lo = 2.5,   hi = 11.8),
  C_PCB170 = rtrunc(n, mean = 5.57,  sd = 5.5,  lo = 2.5,   hi = 38.1),
  C_PCB180 = rtrunc(n, mean = 13.6,  sd = 13,   lo = 2.5,   hi = 87.4),
  C_PCB183 = rtrunc(n, mean = 2.59,  sd = 1.0,  lo = 2.5,   hi = 8.9),
  C_PCB187 = rtrunc(n, mean = 3.83,  sd = 3.5,  lo = 2.5,   hi = 24.4)
)


# =============================================================================
# 3. TOTAL LIPIDS DATASET
# =============================================================================

lipids_synthetic <- data.frame(
  shortCode = subject_ids,
  M_TL_1 = rtrunc(n, mean = 7.48, sd = 2.8, lo = 1.48, hi = 14.06),
  C_TL_1 = rtrunc(n, mean = 2.54, sd = 1.0, lo = 1.48, hi = 9.04)
)

# =============================================================================
# 4. CLINICAL DATA DATASET
# =============================================================================

clinical_synthetic <- data.frame(
  shortCode    = subject_ids,
  AGE_MUM      = sample(17:40, n, replace = TRUE),
  BMI          = round(rnorm(n, mean = 23.55, sd = 3), 2),
  GEST_FIN     = sample(247:291, n, replace = TRUE),
  Genere       = sample(c("maschile", "femminile"), n, replace = TRUE),
  Smoke        = sample(c("si", "no"), n, replace = TRUE),
  Parity       = sample(c("multipara", "primipara"), n, replace = TRUE),
  Titolo_Studio = sample(c("high", "medium", "low"), n, replace = TRUE),
  Parto        = sample(c("cesarean", "vaginal"), n, replace = TRUE)
)

# =============================================================================
# SAVE
# =============================================================================

dir.create("data/example", recursive = TRUE, showWarnings = FALSE)

writexl::write_xlsx(ds_synthetic,   "data/example/synthetic_miRNA.xlsx")
writexl::write_xlsx(poll_synthetic, "data/example/synthetic_poll.xlsx")
writexl::write_xlsx(lipids_synthetic, "data/example/synthetic_lipids.xlsx")
writexl::write_xlsx(clinical_synthetic, "data/example/synthetic_clinical.xlsx")

cat("Datasets saved to data/example/\n")
cat(sprintf("  synthetic_ds.xlsx     : %d subjects x %d miRNAs\n",
            nrow(ds_synthetic), ncol(ds_synthetic) - 1))
cat(sprintf("  synthetic_poll.xlsx   : %d subjects x %d pollutants\n",
            nrow(poll_synthetic), ncol(poll_synthetic) - 1))
cat(sprintf("  synthetic_lipids.xlsx : %d subjects x 2 total lipid variables\n",
            nrow(lipids_synthetic)))

cat(sprintf("  synthetic_clinical.xlsx : %d subjects x %d clinical variables\n",
            nrow(clinical_synthetic), ncol(clinical_synthetic) - 1))
