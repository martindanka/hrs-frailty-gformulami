
# Packages --------------------------------------------------------------------------------------------------------

# Core packages
library(tidyverse)

# Imputation and modelling
library(gFormulaMI)
library(mice)

# Parallel processing
library(future)
library(furrr)
library(parallel)

# Utilities
library(tictoc)

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


# Imputation set up -----------------------------------------------------------------------------------------------
count_missing(data_wide_clean)

### SET UP ###
imp_empty_a <- mice(data_wide_arta, m = 1, maxit = 0, printFlag = FALSE, seed = 55446)
imp_empty_b <- mice(data_wide_artb, m = 1, maxit = 0, printFlag = FALSE, seed = 55446)
imp_empty_c <- mice(data_wide_artc, m = 1, maxit = 0, printFlag = FALSE, seed = 55446)

## Check if all parameters are set correctly
imp_empty_a$predictorMatrix
imp_empty_a$method

### PREDICTION MATRIX ###

# Extract prediction matrices
predmat_a <- imp_empty_a$predictorMatrix
predmat_b <- imp_empty_b$predictorMatrix
predmat_c <- imp_empty_c$predictorMatrix


# Auxiliary variables ---------------------------------------------------------------------------------------------

# Auxiliary variables at t will be used for t+1 only.

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

## Later auxiliary variables
later_aux <- tibble(
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

# Assemble the auxiliary-variable vectors
aux_5_6 <- unlist(auxvar_5_6, use.names = FALSE)
aux_7_13 <- unlist(auxvar_7_13, use.names = FALSE)

aux_vec <- union(aux_5_6, aux_7_13)
aux_vec <- aux_vec[aux_vec %in% colnames(predmat_a)]
later_aux <- unlist(later_aux, use.names = FALSE)


# Update the prediction matrix ------------------------------------------------------------------------------------

# Build predictor-matrix + method vector
update_pred_meth <- function(imp_empty,
                             aux_vec,
                             later_aux,
                             frailty_vec,
                             mari_vec) {

  pm <- imp_empty$predictorMatrix
  meth <- imp_empty$method
  
  # 1 ─ auxiliaries: auxiliary predictor → next wave only  -----------------------------
  pm <- restrict_lag_predictor(pm, aux_vec)
  
  # 2 ─ late-wave auxiliaries: do not use --------------------
  drop <- remove_vars(pm, meth, later_aux)
  pm <- drop$predmat
  meth <- drop$meth
  
  # 3 ─ frailty series  ------------------------------------------------------
  meth[frailty_vec] <- "norm"
  meth[c("frailst", "schlyrs")] <- "norm"
  
  #   a) remove pre-baseline frail5 / frail6
  rm_frail <- remove_vars(pm, meth, c("frail5", "frail6"))
  pm <- rm_frail$predmat
  meth <- rm_frail$meth
  
  #   b) frailty itself is imputed using its own-history but not future values (earlier → later)
  pm <- impose_self_past(pm, frailty_vec)
  
  # 4 ─ marriage series  -----------------------------------------------------
  pm <- impose_self_past(pm, mari_vec)
  pm[c("mari5","mari6"), mari_vec] <- 0L
  
  # 5 ─ household ownership series  -----------------------------------------
  hhown_vec <- paste0("hhown", 7:13)
  pm <- impose_self_past(pm, hhown_vec)
  
  # 6 ─ custom predictors for frailst  --------------------------------------
  # Ensures frailst uses frailty and marital history — already set by default, kept as safeguard
  pm["frailst", c(frailty_vec, mari_vec)] <- 1L
  
  # 7 ─ return ---------------------------------------------------------------
  list(pred = pm, meth = meth)
}

# Common vectors needed for the function above
frailty_vec <- paste0("frail", 7:13)
mari_vec <- paste0("mari", 7:13)


# Common set up ---------------------------------------------------------------------------------------------------

# The original set up on a Macbook with 6 performance and 2 efficiency cores (32 GB RAM)
options(parallelly.fork.enable = "TRUE")
plan(multicore)  # use multicore; always swap to 'multisession' on Windows (see below)
ncore <- 5       # adjust if needed

# General set up for other machines (run this instead)
# plan(multisession)
# ncore <- detectCores() - 1 # Can use -2 instead if you want to continue other work on the machine.

# Number of imputations and iterations
m <- 50
maxit <- 40


# Imputation block for dataset a ----------------------------------------------------------------------------------

# 0. Choose dataset
flag <- "a"
dat_name <- paste0("data_wide_art", flag)
dat <- get(dat_name)

# 1. Empty mids
imp0 <- mice(dat, m = 1, maxit = 0, printFlag = FALSE, seed = 55446)

# 2. Get predictor matrix + method
mods <- update_pred_meth(imp0,
                         aux_vec = aux_vec,
                         later_aux = later_aux,
                         frailty_vec = frailty_vec,
                         mari_vec = mari_vec)

# Computation heavy — comment out this block and load saved RDS on re-runs
# 3. Run the actual imputation
tic()
imp <- progressr::with_progress(
  futuremice(dat,
             maxit = maxit,
             m = m,
             n.core = ncore,
             predictorMatrix = mods$pred,
             method = mods$meth,
             parallelseed = 55446)
)
timing_a <- toc(quiet = TRUE)


outfile <- paste0("data/imp_", flag, "_main.rds")
saveRDS(imp, file = outfile)
rm(imp); gc()

# Imputation block for dataset b ----------------------------------------------------------------------------------

# 0. Choose dataset
flag <- "b"
dat_name <- paste0("data_wide_art", flag)
dat <- get(dat_name)

# 1. Empty mids
imp0 <- mice(dat, m = 1, maxit = 0, printFlag = FALSE, seed = 55446)

# 2. Get predictor matrix + method
mods <- update_pred_meth(imp0,
                         aux_vec = aux_vec,
                         later_aux = later_aux,
                         frailty_vec = frailty_vec,
                         mari_vec = mari_vec)

# Computation heavy — comment out this block and load saved RDS on re-runs
# 3. Run the actual imputation
tic()
imp <- progressr::with_progress(
  futuremice(dat,
             maxit = maxit,
             m = m,
             n.core = ncore,
             predictorMatrix = mods$pred,
             method = mods$meth,
             parallelseed = 55446)
)
timing_b <- toc(quiet = TRUE)


outfile <- paste0("data/imp_", flag, "_main.rds")
saveRDS(imp, file = outfile)
rm(imp); gc()

# Imputation block for dataset c ----------------------------------------------------------------------------------

# 0. Choose dataset
flag <- "c"
dat_name <- paste0("data_wide_art", flag)
dat <- get(dat_name)

# 1. Empty mids
imp0 <- mice(dat, m = 1, maxit = 0, printFlag = FALSE, seed = 55446)

# 2. Get predictor matrix + method
mods <- update_pred_meth(imp0,
                         aux_vec = aux_vec,
                         later_aux = later_aux,
                         frailty_vec = frailty_vec,
                         mari_vec = mari_vec)

# Computation heavy — comment out this block and load saved RDS on re-runs
# 3. Run the actual imputation
tic()
imp <- progressr::with_progress(
  futuremice(dat,
             maxit = maxit,
             m = m,
             n.core = ncore,
             predictorMatrix = mods$pred,
             method = mods$meth,
             parallelseed = 55446)
)
timing_c <- toc(quiet = TRUE)


outfile <- paste0("data/imp_", flag, "_main.rds")
saveRDS(imp, file = outfile)
rm(imp); gc()

# Diagnostic ------------------------------------------------------------------------------------------------------

imp_a <- readRDS("data/imp_a_main.rds")
imp_b <- readRDS("data/imp_b_main.rds")
imp_c <- readRDS("data/imp_c_main.rds")

get_psrf(imp_a)
plot(imp_a)

get_psrf(imp_b)
plot(imp_b)

get_psrf(imp_c)
plot(imp_c)

# The trace plots and PSRFs are acceptable.
# Further improvements could be achieved by restricting the predictor matrix,
# but this risks bias from the imputation models being uncongenial to the analysis.
# We could also use a higher M, but this would increase the computational burden 
# for gains that are likely to be negligible.
