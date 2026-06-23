# =============================================================================
# Bayesian Weighted Quantile Sum (BWQS) Regression
# -----------------------------------------------------------------------------
# Description:
#   Runs BWQS models for three analyses:
#     1. TOT  — full mixture (all metals + POPs together)
#     2. EEs  — endocrine disruptors sub-mixture (metals: Se, Hg, As, Zn, Cu)
#     3. OCs  — organochlorines sub-mixture (HCB, TNC, DDE, PCBs)
#   Models are run for:
#     - miRNAs significant in linear regression (from script 2_0_linearRegression_analysis.R)
#     - All remaining miRNAs (non-significant)
#   For each miRNA, cordonal (_CD) models include "Parto" (which stand for delivery) as confounder;
#   maternal models do not.
#
# Input:
#   - input_dataset.xlsx        (output of 1_miRNA_normalization.R)
#   - regression_MM.xlsx        (output of 2_0_linearRegression_analysis.R)
#   - regression_CC.xlsx        (output of 2_0_linearRegression_analysis.R)
#
# Output (written to output_bWQS/):
#   - FINAL_MODELS_significant_TOT.xlsx
#   - FINAL_MODELS_significant_divided.xlsx
#   - FINAL_MODELS_nonsignificant.xlsx
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(writexl)
library(BWQS)
library(data.table)
library(parallel)
library(foreach)
library(doParallel)
library(progressr)


# --- Load data ---------------------------------------------------------------
DS_TOT <- as.data.frame(readxl::read_excel("./input_dataset.xlsx"))
row.names(DS_TOT) <- DS_TOT$shortCode
DS_TOT$shortCode  <- NULL

df_MM <- readxl::read_excel("./regression_MM.xlsx")
df_CC <- readxl::read_excel("./regression_CC.xlsx")


# --- Define variable groups --------------------------------------------------
clinical_vars <- c("AGE_MUM", "BMI", "GEST_FIN", "Genere", "Smoke",
                   "Parity", "Titolo_Studio", "Parto")

met <- colnames(DS_TOT)[grepl("^M_|^C_|^SUM_", colnames(DS_TOT))]
miR <- colnames(DS_TOT)[!colnames(DS_TOT) %in% c(met, clinical_vars)]

# Maternal and cordonal pollutants (excluding derived variables)
met_M <- met[grepl("^M_", met) & !grepl("^M_WQS", met)]
met_C <- met[grepl("^C_", met) & !grepl("^C_WQS", met)]

# Sub-mixture definitions
# EEs = endocrine disruptors (metals): Se, Zn, Cu
# OCs = organochlorines: HCB, TNC, DDE, PCBs
EEs_M <- met_M[met_M %in% c("M_Se", "M_Zn", "M_Cu")]
OCs_M <- met_M[!met_M %in% EEs_M]
EEs_C <- met_C[met_C %in% c("C_Se", "C_Zn", "C_Cu")]
OCs_C <- met_C[!met_C %in% EEs_C]

# miRNAs significant in linear regression (FDR < 0.05)
miR_sig <- unique(c(
  df_CC$miRNA[which(df_CC$p.adjust.FDR < 0.05)],
  df_MM$miRNA[which(df_MM$p.adjust.FDR < 0.05)]
))
miR_sig    <- unique(gsub("_CD$", "", miR_sig))
miR_sig    <- unique(c(miR_sig, paste0(miR_sig, "_CD")))
miR_sig    <- intersect(miR_sig, miR)

# Non-significant miRNAs
miR_nonsig <- miR[!miR %in% miR_sig]


# --- Helper: build formula string --------------------------------------------
make_formula <- function(mirna) {
  conf <- if (grepl("_CD$", mirna)) {
    "AGE_MUM + BMI + GEST_FIN + Genere + Smoke + Parity + Titolo_Studio + Parto"
  } else {
    "AGE_MUM + BMI + GEST_FIN + Genere + Smoke + Parity + Titolo_Studio"
  }
  paste0("`", mirna, "` ~ ", conf)
}

# --- Helper: select mixture based on miRNA type and analysis -----------------
get_mix <- function(mirna, tp) {
  is_cd <- grepl("_CD$", mirna)
  if (tp == "TOT") return(if (is_cd) met_C else met_M)
  if (tp == "EEs") return(if (is_cd) EEs_C else EEs_M)
  if (tp == "OCs") return(if (is_cd) OCs_C else OCs_M)
}

# --- Helper: run one BWQS model ----------------------------------------------
run_bwqs <- function(mirna, tp, data) {
  mt_nm      <- get_mix(mirna, tp)
  bwqs_model <- bwqs(
    formula  = as.formula(make_formula(mirna)),
    mix_name = mt_nm,
    data     = data,
    q        = NULL,
    seed     = 2025,
    iter     = 1000,
    c_int    = c(0.025, 0.95)
  )
  out          <- as.data.frame(bwqs_model$summary_fit)
  out$mix_name <- gsub("W_", "", row.names(out))
  out$miRNA    <- mirna
  out$ANALYSIS <- tp
  out
}

# --- Parallel setup ----------------------------------------------------------
n_cores <- max(1, detectCores() - 5)
cl      <- makeCluster(n_cores)
registerDoParallel(cl)

dir.create("output_bWQS", recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# SECTION 1: SIGNIFICANT miRNAs — TOT mixture
# =============================================================================

cat("Running BWQS for significant miRNAs (TOT mixture)...\n")

clusterExport(cl, c("miR_sig", "met_C", "met_M", "DS_TOT",
                    "make_formula", "get_mix", "run_bwqs", "EEs_C", "EEs_M", "OCs_C", "OCs_M"))
clusterEvalQ(cl, { library(BWQS); library(data.table) })

RANK_sig_TOT <- foreach(
  cc        = miR_sig,
  .combine  = rbind,
  .packages = c("BWQS", "data.table")
) %dopar% {
  run_bwqs(cc, "TOT", DS_TOT)
}

writexl::write_xlsx(RANK_sig_TOT, "./output_bWQS/FINAL_MODELS_significant_TOT.xlsx")


# =============================================================================
# SECTION 2: SIGNIFICANT miRNAs — EEs and OCs sub-mixtures
# =============================================================================

cat("Running BWQS for significant miRNAs (EEs and OCs)...\n")

jobs_sig <- expand.grid(cc = miR_sig, tp = c("EEs", "OCs"), stringsAsFactors = FALSE)

RANK_sig_divided <- foreach(
  i         = seq_len(nrow(jobs_sig)),
  .combine  = rbind,
  .packages = c("BWQS", "data.table")
) %dopar% {
  run_bwqs(jobs_sig$cc[i], jobs_sig$tp[i], DS_TOT)
}

writexl::write_xlsx(RANK_sig_divided, "./output_bWQS/FINAL_MODELS_significant_divided.xlsx")


# =============================================================================
# SECTION 3: NON-SIGNIFICANT miRNAs — TOT, EEs, OCs
# =============================================================================

cat("Running BWQS for non-significant miRNAs...\n")

jobs_nonsig <- expand.grid(cc = miR_nonsig, tp = c("TOT", "EEs", "OCs"), stringsAsFactors = FALSE)

clusterExport(cl, c("miR_nonsig", "jobs_nonsig"))

RANK_nonsig <- foreach(
  i         = seq_len(nrow(jobs_nonsig)),
  .combine  = rbind,
  .packages = c("BWQS", "data.table", "writexl")
) %dopar% {
  cc       <- jobs_nonsig$cc[i]
  tp       <- jobs_nonsig$tp[i]
  out_file <- paste0("output_bWQS/nonsignificant/bWQS_", cc, "_", tp, ".xlsx")

  # Skip if already computed
  if (file.exists(out_file)) return(NULL)

  res <- run_bwqs(cc, tp, DS_TOT)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(res, out_file)
  res
}

RANK_nonsig <- RANK_nonsig[!is.null(RANK_nonsig), ]
writexl::write_xlsx(RANK_nonsig, "./output_bWQS/FINAL_MODELS_nonsignificant.xlsx")

stopCluster(cl)


# =============================================================================
# SECTION 4: COMBINE RESULTS AND FLAG SIGNIFICANT EFFECTS
# =============================================================================
# A BWQS beta1 effect is considered significant when the credible interval
# does not include zero (both bounds same sign).

flag_significant <- function(df) {
  df$mean_flagged <- ifelse(
    (df$`2.5%` < 0 & df$`95%` < 0) | (df$`2.5%` > 0 & df$`95%` > 0),
    df$mean, NA
  )
  df
}

# Extract beta1 rows only
extract_beta1 <- function(df) subset(df, mix_name == "beta1")

WQS_sig_TOT    <- flag_significant(extract_beta1(RANK_sig_TOT))
WQS_sig_div    <- flag_significant(extract_beta1(RANK_sig_divided))
WQS_nonsig     <- flag_significant(extract_beta1(RANK_nonsig))

# Add label columns for downstream plotting
add_labels <- function(df) {
  df$MET_LIP   <- ifelse(grepl("_CD$", df$miRNA),
                         paste0("C_WQS_MIX_", df$ANALYSIS),
                         paste0("M_WQS_MIX_", df$ANALYSIS))
  df$miRNA_fin <- gsub("_CD$", "", df$miRNA)
  df
}

WQS_sig_TOT <- add_labels(WQS_sig_TOT)
WQS_sig_div <- add_labels(WQS_sig_div)
WQS_nonsig  <- add_labels(WQS_nonsig)

# Final combined dataset
WQS_all <- rbind(WQS_sig_TOT, WQS_sig_div, WQS_nonsig)

writexl::write_xlsx(WQS_all, "./output_bWQS/WQS_all_results.xlsx")

cat("Done. Results saved to output_bWQS/\n")
