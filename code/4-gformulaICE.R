
# Packages --------------------------------------------------------------------------------------------------------

# Core packages
library(tidyverse)
library(furrr)

# Imputation and modelling
library(gFormulaMI)
library(ltmle)
library(mice)

# Reporting and table formatting
library(flextable)
library(officer)

# Utilities
library(tictoc)

# Shared functions used across scripts ---------------------------------------------------------------------------
source("code/helpers.R")

# Import datasets -------------------------------------------------------------------------------------------------
imp_a <- readRDS("data/imp_a_main.rds") # treatment: cultural_events
imp_b <- readRDS("data/imp_b_main.rds") # treatment: music_activity
imp_c <- readRDS("data/imp_c_main.rds") # treatment: craft_hobbies

# Define auxiliaries ----------------------------------------------------------------------------------------------
aux_vars <- {
  aux_5_6 <- tibble(
    base = c("mari","employ","hhown","hincome","nhsafe","frail",
             "cog","moneya","hh2nd"),
    `5` = str_c(base, 5),
    `6` = str_c(base, 6)
  ) %>% select(-base) %>% unlist(use.names = FALSE)
  
  aux_7_13 <- tibble(
    base = c("cog","moneya","hh2nd"),
    `7` = str_c(base, 7),
    `8` = str_c(base, 8),
    `9` = str_c(base, 9),
    `10` = str_c(base, 10),
    `11` = str_c(base, 11),
    `12` = str_c(base, 12),
    `13` = str_c(base, 13)
  ) %>% select(-base) %>% unlist(use.names = FALSE)
  
  union(aux_5_6, aux_7_13)
}

waves <- 7:13
l_vars <- c("mari","employ","hhown","hincome","nhsafe","frail")
y_node <- "frailst"



# Helper functions ------------------------------------------------------------------------------------------------
## Clean one mice object (run once per imputation) ----------------------
clean_imp <- function(imp, treat_prefix, aux_vars = aux_vars) {
  imp_list <- map(seq_len(imp$m), ~complete(imp, .x))
  imp_list |>
    map(~ .x |>
          select(-all_of(aux_vars)) |>
          mutate(across(starts_with(treat_prefix),
                        ~ if_else(as.numeric(.x) == 2, 1L, 0L))))
  }

## Single-run gformulaICE (everything explicit) ------------------------------
fit_single_ltmle <- function(dat,
                             Anodes,
                             Lnodes,
                             Ynodes,
                             abar_treat,
                             abar_control,
                             quiet = TRUE,
                             return_numeric = TRUE) {
  
  runner <- if (quiet) suppressMessages else identity
  
  fit <- runner(ltmle(
    dat,
    Anodes = Anodes,
    Lnodes = Lnodes,
    Ynodes = Ynodes,
    abar = list(treatment = abar_treat, control = abar_control),
    gcomp = TRUE, # Essential toggle — use gformulaICE instead of LTMLE
    estimate.time = FALSE
  ))
  
  if (return_numeric) {
    return(as.numeric(summary(fit)$effect.measures$ATE["estimate"]))
  } else {
    return(fit)
  }
}

## Estimate all strategies for a cleaned list --------------------------

estimate_strategies <- function(df_list,
                                treat_prefix,
                                waves = waves,
                                l_vars = l_vars,
                                y_node = y_node,
                                quiet = TRUE) {
  
  Anodes <- str_c(treat_prefix, waves)
  Lnodes <- unlist(map(waves, \(w) str_c(l_vars, w)))
  control <- rep(0L, length(Anodes))
  
  strategies <- list(
    always = rep(1L, length(Anodes)),
    early = c(rep(1L, 3), rep(0L, length(Anodes) - 3)),
    late = c(rep(0L, length(Anodes) - 3), rep(1L, 3))
  )
  
  imap_dfr(df_list, \(dat, imp_id) {
    imap_dfr(strategies, \(abar_treat, strat) {
      tibble(
        imp = imp_id,
        strategy = strat,
        ATE = fit_single_ltmle(
          dat = dat,
          Anodes = Anodes,
          Lnodes = Lnodes,
          Ynodes = y_node,
          abar_treat = abar_treat,
          abar_control = control,
          quiet = quiet)
      )
    })
  })
}

## Simple mean-pool across imputations ----------------------------------
pool_imputations <- function(res_list) {
  map_dfr(names(res_list), \(nm) {
    res_list[[nm]] %>%
      group_by(strategy) %>%
      summarise(dataset = nm,
                ATE_bar = mean(ATE),
                .groups = "drop")
  })
}


# Optional validation: single run test (can be skipped) -----------------------------------------------------------

## 1. first completed & cleaned data set
dat_example <- clean_imp(imp_a, "cultural_events", aux_vars)[[1]]

## 2. node vectors for the cultural_events exposure
Anodes_ex <- str_c("cultural_events", waves)
Lnodes_ex <- unlist(map(waves, \(w) str_c(l_vars, w)))

## 3. treatment strategy = always treat
abar_treat_ex <- rep(1L, length(Anodes_ex))
abar_ctrl_ex <- rep(0L, length(Anodes_ex))

## 4. run once
ATE_example <- fit_single_ltmle(
  dat = dat_example,
  Anodes = Anodes_ex,
  Lnodes = Lnodes_ex,
  Ynodes = y_node,
  abar_treat = abar_treat_ex,
  abar_control = abar_ctrl_ex,
  quiet = TRUE,
  return_numeric = FALSE # Toggle to see the full output. Default is TRUE, which only returns the numeric point estimate.
)
summary(ATE_example)

# Note: warning can be ignored, this refers to the variance estimation procedure (variance will not be used for this
# sensitivity analysis, CIs would need to be bootstrapped.).

# Run gformulaICE -------------------------------------------------------------------------------------------------
# Computation heavy — on re-runs, skip to the readRDS() call below
prefix_map <- c(imp_a = "cultural_events",
                imp_b = "music_activity",
                imp_c = "craft_hobbies")

# Run for all imputed datasets
clean_imps <- list(
  imp_a = clean_imp(imp_a, prefix_map[["imp_a"]], aux_vars),
  imp_b = clean_imp(imp_b, prefix_map[["imp_b"]], aux_vars),
  imp_c = clean_imp(imp_c, prefix_map[["imp_c"]], aux_vars)
)

# Run gformulaICE
# The original set up
options(parallelly.fork.enable = "TRUE")
plan(multicore, workers = 3) # use multicore; swap to multisession on Windows

# General set up for other machines (run this instead)
# plan(multisession, workers = 3)

# Parallelise and fit to all imputed datasets
ice_results <- future_imap(clean_imps, \(lst, nm) {
  estimate_strategies(
    df_list = lst,
    treat_prefix = prefix_map[[nm]],
    waves = waves,
    l_vars = l_vars,
    y_node = y_node,
    quiet = TRUE
  )
},
.progress = TRUE)
# As before, warnings can be ignored (variance estimation procedure, analytical CIs will not be used).

# Pool across imputations (ICE point estimates)
ice_results <- pool_imputations(ice_results)
ice_results

saveRDS(ice_results, "data/gformula_ice_point.rds")

# Export results --------------------------------------------------------------------------------------------------
ice_results <- readRDS("data/gformula_ice_point.rds")

gform_mi_a <- readRDS("data/gform_mi_a.rds")
gform_mi_b <- readRDS("data/gform_mi_b.rds")
gform_mi_c <- readRDS("data/gform_mi_c.rds")

# Uses fit_pool_gform() from helpers.R
# Labels: "Never treat", "Early exposure", "Late exposure", "Sustained exposure"

results_a <- fit_pool_gform(gform_mi_a)
results_b <- fit_pool_gform(gform_mi_b)
results_c <- fit_pool_gform(gform_mi_c)

ice_results

results_a
results_b
results_c

## 1. put the three MI result objects into a named list
res_list <- list(imp_a = results_a,
                 imp_b = results_b,
                 imp_c = results_c)

## 2. tidy them up, retaining just the three treatment rows
estimates <- imap_dfr(
  res_list,
  ~ .x %>%
    filter(!grepl("Never", strategy)) %>% # drop the intercept row
    mutate(
      strategy = case_when( # match spelling in ice_results
        grepl("Sustained", strategy) ~ "always",
        grepl("Early", strategy) ~ "early",
        grepl("Late", strategy) ~ "late"
      ),
      dataset = .y # .y is the list-name (imp_a …)
    ) %>%
    select(strategy, dataset, Estimate)
)

## 3. append (join) the new column
all_results <- ice_results %>%
  left_join(estimates, by = c("strategy", "dataset")) %>%
  rename(gformulaMI = Estimate,
         gformulaICE = ATE_bar) %>%
  relocate(gformulaMI, .before = gformulaICE) %>%
  relocate(dataset, .before = strategy) |>
  mutate(diff = gformulaMI - gformulaICE)

all_results

# Export as a word table

# ── 1. Relabel variables ───────────────────────────────────────────────────
# a. strategy labels → title-case phrasing
strategy_lookup <- c(
  "never" = "Never treat",
  "early" = "Treat early",
  "late" = "Treat late",
  "always" = "Always treat"
)
# b. dataset → Exposure names
exposure_lookup <- c(
  "imp_a" = "Cultural events",
  "imp_b" = "Singing or instruments",
  "imp_c" = "Arts and crafts"
)

table_df <- all_results %>% 
  mutate(
    Strategy = recode(strategy, !!!strategy_lookup),
    Exposure = recode(dataset, !!!exposure_lookup)
  ) %>% 
  select(Exposure, Strategy, gformulaMI, gformulaICE, diff) %>% # order cols
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# ── 2. Build a Word table with officer + flextable ─────────────────────────
ft <- flextable(table_df) |>
  set_header_labels(
    Exposure = "Exposure",
    Strategy = "Strategy",
    gformulaMI = "gformulaMI estimate",
    gformulaICE = "gformula ICE estimate",
    diff = "Difference"
  ) |>
  autofit()

# ── 3. Write the .docx file ────────────────────────────────────────────────
doc <- read_docx() |> body_add_flextable(ft)

print(doc, target = "output/results_ice.docx")