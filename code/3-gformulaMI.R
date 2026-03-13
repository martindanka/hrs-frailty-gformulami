
# Packages --------------------------------------------------------------------------------------------------------

# Core packages
library(tidyverse)

# Imputation and modelling
library(gFormulaMI)
library(mice)

# Reporting and table formatting
library(flextable)
library(officer)

# Utilities
library(tictoc)
library(gridExtra)

# Shared functions used across scripts ---------------------------------------------------------------------------
## Before sourcing, ensure your working directory is set to the project root (i.e. the HRS frailty folder)
## Example: setwd("/Users/yourname/Documents/HRS Frailty")
source("code/helpers.R")

# Prepare data ----------------------------------------------------------------------------------------------------

data_wide_clean <- readRDS("data/data_wide_clean.rds")

data_wide_main <- data_wide_clean |>
  dplyr::select(-c(pid), -starts_with(c("shlt", "psyche", "conde")))

data_wide_arta <- data_wide_main |>
  dplyr::select(-starts_with(c("music_activity", "craft_hobbies")))

data_wide_artb <- data_wide_main |>
  dplyr::select(-starts_with(c("cultural_events", "craft_hobbies")))

data_wide_artc <- data_wide_main |>
  dplyr::select(-starts_with(c("cultural_events", "music_activity")))

# Import imputed datasets
imp_a <- readRDS("data/imp_a_main.rds")
imp_b <- readRDS("data/imp_b_main.rds")
imp_c <- readRDS("data/imp_c_main.rds")

# Model specification (gFormulaMI) --------------------------------------------------------------------------------

# 1) Run gFormulaImpute (empty run) to get the default predictor matrix
#    - This will only be used to specify the correct models.

# Create M = 1 imputed dataset with maxit = 0
# This will only be used to extract the predictor matrix
imp_helper_a <- mice(data_wide_arta, maxit = 0, m = 1)
imp_helper_b <- mice(data_wide_artb, maxit = 0, m = 1)
imp_helper_c <- mice(data_wide_artc, maxit = 0, m = 1)

# Run gFormulaImpute as a helper to extract the predictor matrix
gimp_single_arta <- gFormulaImpute(
  data = imp_helper_a,
  M = 1,
  trtVars = c("cultural_events7","cultural_events8","cultural_events9", "cultural_events10","cultural_events11",
              "cultural_events12", "cultural_events13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1))
)

gimp_single_artb <- gFormulaImpute(
  data = imp_helper_b,
  M = 1,
  trtVars = c("music_activity7","music_activity8","music_activity9", "music_activity10","music_activity11",
              "music_activity12", "music_activity13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1))
)

gimp_single_artc <- gFormulaImpute(
  data = imp_helper_c,
  M = 1,
  trtVars = c("craft_hobbies7","craft_hobbies8","craft_hobbies9", "craft_hobbies10","craft_hobbies11","craft_hobbies12", 
              "craft_hobbies13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1))
)

# Extract the matrix
predmat_a <- gimp_single_arta$predictorMatrix
predmat_b <- gimp_single_artb$predictorMatrix
predmat_c <- gimp_single_artc$predictorMatrix

# Specify the models

# Auxiliary variables should not be used
## 1 ────────────────────────────────────────────────────────────────
## Collect the auxiliary-variable names once, keep only uniques

## From sweeps 5-6
auxvar_5_6 <- tibble(
  base = c("mari", "employ", "hhown", "hincome", "nhsafe", "frail", "cog", "moneya", "hh2nd"),
  wave5 = paste0(base, 5),
  wave6 = paste0(base, 6)
) |>
  dplyr::select(-base)

## From sweeps 7-13
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

aux_vars <- unique(c(unlist(auxvar_5_6),
                     unlist(auxvar_7_13)))

# Disable the auxiliary variables 
meth_a <- gimp_single_arta$method
meth_b <- gimp_single_artb$method
meth_c <- gimp_single_artc$method

pred_meth_a <- remove_vars(predmat_a, meth = meth_a, vec = aux_vars)
pred_meth_b <- remove_vars(predmat_b, meth = meth_b, vec = aux_vars)
pred_meth_c <- remove_vars(predmat_c, meth = meth_c, vec = aux_vars)

# Store new predictor matrices and methods vectors separately
predmat_a <- pred_meth_a$predmat
predmat_b <- pred_meth_b$predmat
predmat_c <- pred_meth_c$predmat
meth_a <- pred_meth_a$meth
meth_b <- pred_meth_b$meth
meth_c <- pred_meth_c$meth

# Estimation (gFormulaMI) -----------------------------------------------------------------------------------------

# Seed for g-formula estimation
set.seed(46826)
gform_mi_a <- gFormulaImpute(
  data = imp_a, 
  M = 50,
  trtVars = c("cultural_events7","cultural_events8","cultural_events9", "cultural_events10","cultural_events11",
              "cultural_events12", "cultural_events13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1)),
  predictorMatrix = predmat_a,
  method = meth_a[names(meth_a) != "regime"]
)

gform_mi_b <- gFormulaImpute(
  data = imp_b, 
  M = 50,
  trtVars = c("music_activity7","music_activity8","music_activity9", "music_activity10","music_activity11",
              "music_activity12", "music_activity13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1)),
  predictorMatrix = predmat_b,
  method = meth_b[names(meth_b) != "regime"]
)

gform_mi_c <- gFormulaImpute(
  data = imp_c, 
  M = 50,
  trtVars = c("craft_hobbies7","craft_hobbies8","craft_hobbies9", "craft_hobbies10","craft_hobbies11","craft_hobbies12", 
              "craft_hobbies13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1)),
  predictorMatrix = predmat_c,
  method = meth_c[names(meth_c) != "regime"]
)

# Save files
saveRDS(gform_mi_a, "data/gform_mi_a.rds")
saveRDS(gform_mi_b, "data/gform_mi_b.rds")
saveRDS(gform_mi_c, "data/gform_mi_c.rds")

# On re-runs, comment out estimation above and load from saved results instead
gform_mi_a <- readRDS("data/gform_mi_a.rds")
gform_mi_b <- readRDS("data/gform_mi_b.rds")
gform_mi_c <- readRDS("data/gform_mi_c.rds")

results_a <- fit_pool_gform(gform_mi_a)
results_b <- fit_pool_gform(gform_mi_b)
results_c <- fit_pool_gform(gform_mi_c)

results_a
results_b
results_c

# Export results

export_results_doc(
  results_list = list(
    "Results for attending concerts, movies, lectures, or visiting museums" = results_a,
    "Results for singing or playing instruments" = results_b,
    "Results for arts and crafts" = results_c
  ),
  filename = "output/results_main.docx"
)

# Plot results ----------------------------------------------------------------------------------------------------

results_list <- list(results_a, results_b, results_c)
p <- plot_contrasts(results_list)
p

# Export the plot
ggsave(
  filename = "output/plot_main.pdf",
  plot = p,
  width = 6,
  height = 4
)

ggsave(
  filename = "output/plot_main.eps",
  plot = p,
  width = 6,
  height = 4
)

ggsave(
  filename = "output/plot_main.png",
  plot = p,
  width = 6,
  height = 4,
  dpi = 600
)

# Sensitivity analysis: Excluding frailty -------------------------------------------------------------------------
frailty_vec <- paste0("frail", 7:13)
a_nofrail <- remove_vars(predmat_a, meth_a, vec = frailty_vec)
b_nofrail <- remove_vars(predmat_b, meth_b, vec = frailty_vec)
c_nofrail <- remove_vars(predmat_c, meth_c, vec = frailty_vec)

set.seed(46826)
gform_mi_a_nofrail <- gFormulaImpute(
  data = imp_a,
  M = 50,
  trtVars = c("cultural_events7","cultural_events8","cultural_events9", "cultural_events10","cultural_events11",
              "cultural_events12", "cultural_events13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1)),
  predictorMatrix = a_nofrail$predmat,
  method = a_nofrail$meth[names(a_nofrail$meth) != "regime"]
)

gform_mi_b_nofrail <- gFormulaImpute(
  data = imp_b,
  M = 50,
  trtVars = c("music_activity7","music_activity8","music_activity9", "music_activity10","music_activity11",
              "music_activity12", "music_activity13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1)),
  predictorMatrix = b_nofrail$predmat,
  method = b_nofrail$meth[names(b_nofrail$meth) != "regime"]
)

gform_mi_c_nofrail <- gFormulaImpute(
  data = imp_c,
  M = 50,
  trtVars = c("craft_hobbies7","craft_hobbies8","craft_hobbies9", "craft_hobbies10","craft_hobbies11","craft_hobbies12",
              "craft_hobbies13"),
  trtRegimes = list(c(0,0,0,0,0,0,0), c(1,1,1,0,0,0,0), c(0,0,0,0,1,1,1), c(1,1,1,1,1,1,1)),
  predictorMatrix = c_nofrail$predmat,
  method = c_nofrail$meth[names(c_nofrail$meth) != "regime"]
)

saveRDS(gform_mi_a_nofrail, "data/gform_mi_a_nofrail.rds")
saveRDS(gform_mi_b_nofrail, "data/gform_mi_b_nofrail.rds")
saveRDS(gform_mi_c_nofrail, "data/gform_mi_c_nofrail.rds")

# On re-runs, comment out estimation above and load from saved results instead
gform_mi_a_nofrail <- readRDS("data/gform_mi_a_nofrail.rds")
gform_mi_b_nofrail <- readRDS("data/gform_mi_b_nofrail.rds")
gform_mi_c_nofrail <- readRDS("data/gform_mi_c_nofrail.rds")

results_a_nofrail <- fit_pool_gform(gform_mi_a_nofrail)
results_b_nofrail <- fit_pool_gform(gform_mi_b_nofrail)
results_c_nofrail <- fit_pool_gform(gform_mi_c_nofrail)

results_a_nofrail
results_b_nofrail
results_c_nofrail

# Export results
export_results_doc(
  results_list = list(
    "Results for attending concerts, movies, lectures, or visiting museums" = results_a_nofrail,
    "Results for singing or playing instruments" = results_b_nofrail,
    "Results for arts and crafts" = results_c_nofrail
  ),
  filename = "output/results_nofrail.docx"
)

# Plot results and export 
nofrail_list <- list(results_a_nofrail, results_b_nofrail, results_c_nofrail)
p <- plot_contrasts(nofrail_list)
p

ggsave(
  filename = "output/plot_nofrail.pdf",
  plot = p,
  width = 6,
  height = 4
)

ggsave(
  filename = "output/plot_nofrail.eps",
  plot = p,
  width = 6,
  height = 4
)

ggsave(
  filename = "output/plot_nofrail.png",
  plot = p,
  width = 6,
  height = 4, 
  dpi = 600
)

# Combined plot: Model I (a) + Model II (b) --------------------------------------------------
p_nofrail <- plot_contrasts(nofrail_list) +
  ggtitle("(a)") +
  theme(plot.title = element_text(size = 12, face = "plain", hjust = 0.5))

p_main <- plot_contrasts(results_list) +
  ggtitle("(b)") +
  theme(plot.title = element_text(size = 12, face = "plain", hjust = 0.5))

p_main_combined <- arrangeGrob(p_nofrail, p_main, nrow = 2)

ggsave("output/plot_main_combined.pdf", p_main_combined, width = 6, height = 8)
ggsave("output/plot_main_combined.eps", p_main_combined, width = 6, height = 8)
ggsave("output/plot_main_combined.png", p_main_combined, width = 6, height = 8, dpi = 600)
ggsave("output/plot_main_combined.svg", p_main_combined, width = 6, height = 8)