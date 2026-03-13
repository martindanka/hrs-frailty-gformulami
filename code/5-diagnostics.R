
# Packages --------------------------------------------------------------------------------------------------------

# Core packages
library(tidyverse)

# Imputation and modelling
library(gFormulaMI)
library(mice)

# Plot layout
library(grid)
library(gridExtra)

# Shared functions used across scripts ---------------------------------------------------------------------------
## Before sourcing, ensure your working directory is set to the project root (i.e. the HRS frailty folder)
## Example: setwd("/Users/yourname/Documents/HRS Frailty")
source("code/helpers.R")

# Reconstruct the g-formula predictor matrix ----------------------------------------------------------------------
# Replicates the setup from 3-gformulaMI.R to obtain the predictor matrix used in the g-formula.
# This is DIFFERENT from the MICE imputation predictor matrix in 2-imputations.R.

data_wide_clean <- readRDS("data/data_wide_clean.rds")

data_wide_main <- data_wide_clean |>
  dplyr::select(-c(pid), -starts_with(c("shlt", "psyche", "conde")))

# Use exposure A as representative (frailst predictor structure is identical across exposures)
data_wide_arta <- data_wide_main |>
  dplyr::select(-starts_with(c("music_activity", "craft_hobbies")))

# Empty MICE run + helper gFormulaImpute to extract the default predictor matrix
imp_helper_a <- mice(data_wide_arta, maxit = 0, m = 1)

gimp_single_arta <- gFormulaImpute(
  data = imp_helper_a,
  M = 1,
  trtVars = c("cultural_events7", "cultural_events8", "cultural_events9",
              "cultural_events10", "cultural_events11", "cultural_events12", "cultural_events13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1))
)

# Extract predictor matrix and method vector
predmat_a <- gimp_single_arta$predictorMatrix
meth_a <- gimp_single_arta$method

# Remove auxiliary variables (same approach as the auxiliary variable removal in 3-gformulaMI.R)
auxvar_5_6 <- tibble(
  base = c("mari", "employ", "hhown", "hincome", "nhsafe", "frail", "cog", "moneya", "hh2nd"),
  wave5 = paste0(base, 5),
  wave6 = paste0(base, 6)
) |>
  dplyr::select(-base)

auxvar_7_13 <- tibble(
  base = c("cog", "moneya", "hh2nd"),
  wave7 = paste0(base, 7),
  wave8 = paste0(base, 8),
  wave9 = paste0(base, 9),
  wave10 = paste0(base, 10),
  wave11 = paste0(base, 11),
  wave12 = paste0(base, 12),
  wave13 = paste0(base, 13)
) |>
  dplyr::select(-base)

aux_vars <- unique(c(unlist(auxvar_5_6), unlist(auxvar_7_13)))

pred_meth_a <- remove_vars(predmat_a, meth = meth_a, vec = aux_vars)
predmat_a <- pred_meth_a$predmat
meth_a <- pred_meth_a$meth

# Extract frailst predictors from the g-formula predictor matrix
frailst_preds_base <- names(which(predmat_a["frailst", ] == 1))
frailst_preds_base <- setdiff(frailst_preds_base, "regime") # regime is added by gFormulaMI, not in the imputed data

cat("Frailst predictors in the g-formula outcome model:\n")
print(frailst_preds_base)

# Exposure configuration ------------------------------------------------------------------------------------------
# The outcome models share the same covariates; only the exposure variable differs.
# We extract predictors from exposure A and swap the exposure variable for B and C.

exp_a_vars <- paste0("cultural_events", 7:13)

exposures <- list(
  list(flag = "a", label = "Cultural events",
       supp_label = "Concerts/movies/lectures",
       exp_new = NULL),
  list(flag = "b", label = "Music activity",
       supp_label = "Sing/play music",
       exp_new = paste0("music_activity", 7:13)),
  list(flag = "c", label = "Arts and crafts",
       supp_label = "Arts and crafts",
       exp_new = paste0("craft_hobbies", 7:13))
)

# Helper: complete one imputation, fit the outcome model, return diagnostics
extract_diagnostics <- function(imp, m, preds) {
  dat <- mice::complete(imp, m)
  fit <- lm(frailst ~ ., data = dat[, c("frailst", preds)])
  tibble(
    imputation = m,
    observed = dat$frailst,
    fitted = fitted(fit),
    residuals = residuals(fit)
  )
}

# Shared plot elements --------------------------------------------------------------------------------------------

theme_diag <- theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7)
  )

theme_supp <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 11),
    strip.text = element_text(size = 11)
  )

# Plot functions --------------------------------------------------------------------------------------------------

make_resid_plot <- function(df, ncol = 5) {
  ggplot(df, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.2, size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue", linewidth = 0.5) +
    facet_wrap(~ paste0("m = ", imputation), ncol = ncol) +
    labs(x = "Fitted", y = "Residuals", title = "Residuals vs Fitted") +
    theme_diag
}

make_calib_plot <- function(df, ncol = 5) {
  ggplot(df, aes(x = fitted, y = observed)) +
    geom_point(alpha = 0.2, size = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue", linewidth = 0.5) +
    facet_wrap(~ paste0("m = ", imputation), ncol = ncol) +
    labs(x = "Predicted", y = "Observed", title = "Observed vs Predicted") +
    theme_diag
}

make_decile_plot <- function(df, ncol = 5, point_size = 1) {
  decile_df <- df |>
    group_by(imputation) |>
    mutate(decile = ntile(fitted, 10)) |>
    group_by(imputation, decile) |>
    summarise(
      mean_predicted = mean(fitted),
      mean_observed = mean(observed),
      .groups = "drop"
    )

  ggplot(decile_df, aes(x = mean_predicted, y = mean_observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_point(size = point_size, colour = "steelblue") +
    facet_wrap(~ paste0("m = ", imputation), ncol = ncol) +
    labs(x = "Mean predicted (decile)", y = "Mean observed (decile)",
         title = "Decile Calibration") +
    theme_diag
}

# Create output directories ---------------------------------------------------------------------------------------

dir.create("output/diagnostics/exploratory", recursive = TRUE, showWarnings = FALSE)
dir.create("output/diagnostics/supplement", recursive = TRUE, showWarnings = FALSE)

# Selected imputations for the supplement figure
supplement_m <- c(1, 2, 3, 4)

# Toggle: include the main header in the compact supplement figure?
show_header <- TRUE

# Loop over exposures ---------------------------------------------------------------------------------------------

walk(exposures, \(exp) {
  flag <- exp$flag
  label <- exp$label
  exp_new <- exp$exp_new

  cat("\n--- Exposure", toupper(flag), "(", label, ") ---\n")

  # 1. Build predictor set for this exposure
  preds <- frailst_preds_base
  if (!is.null(exp_new)) {
    # Replace exposure A variable names with this exposure's names
    idx <- match(exp_a_vars, preds)
    idx_present <- !is.na(idx)
    preds[idx[idx_present]] <- exp_new[idx_present]
  }

  # 2. Load imputed data
  imp <- readRDS(paste0("data/imp_", flag, "_main.rds"))
  n_imp <- imp$m

  # Keep only predictors present in the completed data
  preds <- intersect(preds, colnames(mice::complete(imp, 1)))
  cat("Number of predictors:", length(preds), "\n")

  # 3. Extract diagnostics across all imputations
  all_diag <- map_dfr(seq_len(n_imp), \(m) extract_diagnostics(imp, m, preds))

  # 4. Exploratory output: all imputations (PNG) -------------------------------------------------------------------
  # PNG is used because the faceted scatter plots are too dense for vector PDF
  # (thousands of points × 50 panels overwhelms the PDF device on macOS).

  ggsave(paste0("output/diagnostics/exploratory/residuals_", flag, ".png"),
         make_resid_plot(all_diag), width = 12, height = 24, dpi = 150)

  ggsave(paste0("output/diagnostics/exploratory/calibration_", flag, ".png"),
         make_calib_plot(all_diag), width = 12, height = 24, dpi = 150)

  ggsave(paste0("output/diagnostics/exploratory/decile_", flag, ".png"),
         make_decile_plot(all_diag), width = 12, height = 24, dpi = 150)


  cat("Saved: exploratory plots for exposure", toupper(flag), "\n")

  # 5. Supplement figures: selected imputations --------------------------------------------------------------------
  supp_diag <- all_diag |>
    filter(imputation %in% supplement_m) |>
    mutate(imputation = factor(
      paste0("Imputation ", imputation),
      levels = paste0("Imputation ", supplement_m)
    ))

  # Row 1: Residuals vs Fitted
  p_resid_supp <- ggplot(supp_diag, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue", linewidth = 0.7) +
    facet_wrap(~ imputation, ncol = 4) +
    labs(x = "Fitted values", y = "Residuals") +
    theme_supp

  # Row 2: Observed vs Predicted (scatter)
  p_calib_supp <- ggplot(supp_diag, aes(x = fitted, y = observed)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_smooth(method = "loess", se = FALSE, colour = "steelblue", linewidth = 0.7) +
    facet_wrap(~ imputation, ncol = 4) +
    labs(x = "Predicted values", y = "Observed values") +
    theme_supp

  # Row 3: Decile calibration
  supp_decile <- supp_diag |>
    group_by(imputation) |>
    mutate(decile = ntile(fitted, 10)) |>
    group_by(imputation, decile) |>
    summarise(
      mean_predicted = mean(fitted),
      mean_observed = mean(observed),
      .groups = "drop"
    )

  p_decile_supp <- ggplot(supp_decile, aes(x = mean_predicted, y = mean_observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_point(size = 1.5, colour = "steelblue") +
    facet_wrap(~ imputation, ncol = 4) +
    labs(x = "Mean predicted (decile)", y = "Mean observed (decile)") +
    theme_supp

  supp_title <- paste0("Outcome model diagnostics \u2014 ", label, " (frailst)")

  # 3-row supplement (exploratory folder, PDF only)
  p_supplement_full <- arrangeGrob(
    p_resid_supp, p_calib_supp, p_decile_supp,
    nrow = 3,
    top = textGrob(supp_title, gp = gpar(fontsize = 12, fontface = "bold"))
  )

  ggsave(
    filename = paste0("output/diagnostics/exploratory/supplement_", flag, ".pdf"),
    plot = p_supplement_full,
    width = 12,
    height = 9
  )

  cat("Saved: 3-row supplement (exploratory) for exposure", toupper(flag), "\n")

  # Compact 2-row supplement (supplement folder, PDF + EPS + PNG)
  supp_title_compact <- paste0("Frailty outcome model diagnostics \u2013 ", exp$supp_label)

  p_resid_compact <- p_resid_supp +
    ggtitle("(a) Residuals vs fitted values")

  p_decile_compact <- p_decile_supp +
    ggtitle("(b) Decile-based calibration plot")

  p_supplement_compact <- arrangeGrob(
    p_resid_compact, p_decile_compact,
    nrow = 2,
    top = if (show_header) textGrob(supp_title_compact, gp = gpar(fontsize = 14, fontface = "bold")) else NULL
  )

  walk(c("pdf", "eps", "png"), \(ext) {
    ggsave(
      filename = paste0("output/diagnostics/supplement/supplement_compact_", flag, ".", ext),
      plot = p_supplement_compact,
      width = 12,
      height = 6,
      dpi = if (ext == "png") 600 else 72
    )
  })

  cat("Saved: compact supplement (PDF/EPS/PNG) for exposure", toupper(flag), "\n")

  # Free memory
  rm(imp, all_diag, supp_diag, supp_decile)
  gc()
})

cat("\nAll diagnostics complete.\n")
