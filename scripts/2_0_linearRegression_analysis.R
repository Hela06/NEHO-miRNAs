# =============================================================================
# Linear Regression Analysis: miRNA ~ Pollutants
# -----------------------------------------------------------------------------
# Description:
#   For each miRNA-pollutant pair, fits a linear regression model adjusted
#   for confounders. Results are split into three comparisons:
#     - MM: maternal miRNA ~ maternal pollutants
#     - CC: cordonal miRNA ~ cordonal pollutants
#     - MC: cordonal miRNA ~ maternal pollutants
#   FDR correction is applied per miRNA within each comparison.
#
# Input:
#   - input_dataset.xlsx   (output of 1.miRNA_normalization.R)
#
# Output:
#   - regression_MM.xlsx   (maternal miRNA ~ maternal pollutants)
#   - regression_CC.xlsx   (cordonal miRNA ~ cordonal pollutants)
#   - regression_MC.xlsx   (cordonal miRNA ~ maternal pollutants)
#
# Note: For maternal miRNA models, "Parto" is excluded from confounders
#       as it is a birth outcome not relevant to maternal exposure models.
#
# Author: Ilaria Cosentini
# Date:   June 2026 (last update)
# =============================================================================


# --- Libraries ---------------------------------------------------------------
library(readxl)
library(writexl)
library(dplyr)
library(broom)
library(data.table)


# --- Load data ---------------------------------------------------------------
DS_TOT <- as.data.frame(readxl::read_excel("./input_dataset.xlsx"))
row.names(DS_TOT) <- DS_TOT$shortCode
DS_TOT$shortCode  <- NULL


# --- Define variable groups --------------------------------------------------
# miRNA columns: all columns not matching pollutant prefixes or clinical names
clinical_vars <- c("AGE_MUM", "BMI", "GEST_FIN", "Genere", "Smoke",
                   "Parity", "Titolo_Studio", "Parto")

met <- colnames(DS_TOT)[grepl("^M_|^C_|^SUM_", colnames(DS_TOT))]
miR <- colnames(DS_TOT)[!colnames(DS_TOT) %in% c(met, clinical_vars)]

# Split pollutants into maternal and cordonal
met_M <- met[grepl("^M_", met)]
met_C <- met[grepl("^C_", met)]

confounders <- clinical_vars


# =============================================================================
# LINEAR REGRESSION: miRNA ~ pollutant + confounders
# =============================================================================

set.seed(2024)
df_fin <- list()
cat("Running linear regressions...\n")

for (cc in miR) {
  for (mm in met) {

    # Exclude "Parto" from confounders for maternal miRNA models
    conf <- if (!grepl("_CD$", cc)) confounders[confounders != "Parto"] else confounders

    cols     <- c(cc, mm, conf)
    reg_data <- na.omit(DS_TOT[, cols, drop = FALSE])

    model    <- lm(reg_data[, cc] ~ ., data = reg_data[, -1, drop = FALSE])
    adjusted <- broom::tidy(model, conf.int = TRUE)
    adjusted$miRNA   <- cc
    adjusted$pollutant <- mm

    df_fin[[length(df_fin) + 1]] <- adjusted
  }
}

df_fin <- do.call(rbind, df_fin)


# --- Keep only the pollutant term (remove intercept and confounders) ---------
coef_to_remove <- c("(Intercept)", confounders,
                    "GenereMaschile", "Smokesi", "PartoVAGINAL",
                    "Titolo_StudioMEDIUM", "Titolo_StudioLOW", "ParityPRIMIPARA")

df_fin_sub <- subset(df_fin, !term %in% coef_to_remove)


# =============================================================================
# SPLIT INTO MM / CC / MC AND APPLY FDR CORRECTION
# =============================================================================

fdr_adjust <- function(df, mirna_col = "miRNA", pval_col = "p.value") {
  df$p.adjust.FDR <- NA
  for (mir in unique(df[[mirna_col]])) {
    idx <- df[[mirna_col]] == mir
    df$p.adjust.FDR[idx] <- p.adjust(df[[pval_col]][idx], method = "fdr")
  }
  df
}

# Maternal miRNA ~ maternal pollutants
df_MM <- df_fin_sub[!grepl("_CD$", df_fin_sub$miRNA) & df_fin_sub$pollutant %in% met_M, ]
df_MM <- fdr_adjust(df_MM)

# Cordonal miRNA ~ cordonal pollutants
df_CC <- df_fin_sub[grepl("_CD$", df_fin_sub$miRNA) & df_fin_sub$pollutant %in% met_C, ]
df_CC <- fdr_adjust(df_CC)

# Cordonal miRNA ~ maternal pollutants
df_MC <- df_fin_sub[grepl("_CD$", df_fin_sub$miRNA) & df_fin_sub$pollutant %in% met_M, ]
df_MC <- fdr_adjust(df_MC)


# --- Save results ------------------------------------------------------------
writexl::write_xlsx(df_MM, "regression_MM.xlsx")
writexl::write_xlsx(df_CC, "regression_CC.xlsx")
writexl::write_xlsx(df_MC, "regression_MC.xlsx")

cat("Done.\n")
cat(sprintf("  MM pairs: %d\n", nrow(df_MM)))
cat(sprintf("  CC pairs: %d\n", nrow(df_CC)))
cat(sprintf("  MC pairs: %d\n", nrow(df_MC)))
