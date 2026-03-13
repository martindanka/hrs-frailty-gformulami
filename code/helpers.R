
# Pseudorandom seed -----------------------------------------------------------------------------------------------
pseudorandom_seed <- function() {
  initial_seed <- Sys.time() 
  initial_seed <- as.integer(initial_seed)
  the_seed <- initial_seed %% 100000
  print(the_seed)
}

# Descriptives ----------------------------------------------------------------------------------------------------
# Check proportions of missing data in each variable
count_missing <- function(data) {
  count.missing <- as.matrix(sapply(data, function(data) sum(is.na(data))))
  percent.missing <- as.matrix(round(count.missing / as.numeric(count(data)) * 100, 1))
  missing_percent_combined <- as.matrix(cbind(count.missing, percent.missing))

  count.complete <- as.matrix(sapply(data, function(data) sum(!is.na(data))))
  percent.complete <- as.matrix((count.complete / as.numeric(count(data))) * 100)
  complete_percent_combined <- as.matrix(cbind(count.complete, percent.complete))
  missing_var <- cbind.data.frame(missing_percent_combined, complete_percent_combined)
  missing_var$Names <- names(data)
  colnames(missing_var) <- c("N missing", "% missing", "N complete", "% complete", "Names")
  (missing_var <- missing_var[c("Names", setdiff(names(missing_var), "Names"))])
  missing_var[order(missing_var$`N missing`, decreasing = TRUE), ]
}

# Imputation helpers ----------------------------------------------------------------------------------------------

# Disable variables from being used as part of imputation (or gformulaMI)
remove_vars <- function(predmat, meth, vec = NULL){
  # Disable the variables in the predictor matrix
  predmat[, vec] <- 0
  predmat[vec, ] <- 0
  # Disable the variables in the methods vector
  meth[vec] <- ""
  # Return a new list with the predictor matrix and the methods vector
  return(list(predmat = predmat, meth = meth))
}

# Utility - get wave
get_wave <- function(x) as.numeric(stringr::str_extract(x, "\\d+$"))

### RESTRICT A PREDICTOR TO USING ITS LAG 1 ONLY (e.g. PRED T2 -> ALL VARS T3) ###
restrict_lag_predictor <- function(pmat, aux_vars) {
  row_waves <- get_wave(rownames(pmat))          # wave of each *target* variable
  for (aux in aux_vars) {
    # skip silently if the auxiliary is not in this particular data set
    if (!aux %in% colnames(pmat)) next
    aux_wave <- get_wave(aux)
    # start by switching the whole column off …
    pmat[ , aux] <- 0L
    # … then switch it back on only for targets measured at wave (n + 1)
    if (!is.na(aux_wave)) {
      pmat[row_waves == (aux_wave + 1L), aux] <- 1L
    }
  }
  pmat
}

### RESTRICT ALL TIME-VARYING PREDICTORS TO USING LAG FOR IMPUTING TARGET VARIABLE(S) (e.g. ALL PREDS T2 -> VAR X T3)
restrict_lags_target <- function(
    target_vars,
    pred,
    lag = 1L,
    other_tv = character()) {
  
  ## ── helper: extract trailing digits ───────────────────────────────────
  wave_number <- function(x) {
    as.integer(sub(".*?(\\d+)$", "\\1", x, perl = TRUE))
  }
  
  ## ── checks ────────────────────────────────────────────────────────────
  all_tv <- unique(c(target_vars, other_tv))
  stopifnot(is.matrix(pred),
            all(target_vars %in% rownames(pred)),
            all(all_tv      %in% colnames(pred)),
            lag >= 0L)
  
  tgt_wave <- wave_number(target_vars)
  oth_wave <- wave_number(other_tv)
  
  if (anyNA(tgt_wave) || anyNA(oth_wave))
    stop("Every time-varying variable must end in digits denoting its wave.")
  
  ## ── main loop: alter only target rows ─────────────────────────────────
  for (j in seq_along(target_vars)) {
    
    tgt            <- target_vars[j]
    current_wave   <- tgt_wave[j]
    
    ## 1.  Zero out *all* time-varying predictors for this row
    pred[tgt, all_tv] <- 0L
    
    ## 2.  Identify waves inside the allowed retrospective window
    if (lag > 0L) {
      allowed_waves <- (current_wave - lag):(current_wave - 1L)
      ## enable own-series history
      pred[tgt, target_vars[tgt_wave %in% allowed_waves]] <- 1L
      ## enable other time-varying series history
      pred[tgt, other_tv[oth_wave %in% allowed_waves]] <- 1L
    }
    ## (mice keeps the diagonal zero, so no need to touch it)
  }
  
  pred
}


### IMPOSE A TIME LAG FOR A GIVEN VARIABLE (e.g. VAR X AT T2 -> VAR X AT T3) ###
# Note: not used in the current pipeline; kept for reference / future use.
impose_self_lag <- function(pmat, vars, lag = 1) {
  ## basic checks
  stopifnot(is.matrix(pmat),
            all(vars %in% rownames(pmat)),
            all(vars %in% colnames(pmat)),
            lag >= 0L)
  
  ## loop over each time-point
  for (j in seq_along(vars)) {
    target <- vars[j]
    
    ## zero-out all predictors drawn from the vars series
    pmat[target, vars] <- 0L
    
    ## admissible retrospective window
    if (lag > 0L && j > 1L) {
      start_idx   <- max(1L, j - lag)
      allowed     <- vars[start_idx:(j - 1L)]
      pmat[target, allowed] <- 1L
    }
    ## (by convention mice keeps the diagonal at zero, so no change needed)
  }
  pmat
}

### DISABLE FUTURE VALUES FROM IMPUTING PAST VALUES OF A VARIABLE 
# (e.g. VAR X T1 + VAR X T2 + ... + VAR X Tn-> VAR X T(n+1) )
# For a given variable, disable its future values from being used to impute its past values.
# Supply a vector with this variable and the time indices (e.g., "frail7", "frail8", etc.)
impose_self_past <- function(pmat, vars) {
  # Extract numeric suffixes and sort vars by them
  nums <- as.integer(sub(".*?(\\d+)$", "\\1", vars))
  vars_sorted <- vars[order(nums)]
  
  # Get indices of those variables in the predictor‐matrix
  idx <- match(vars_sorted, rownames(pmat))
  
  # For each “kth” variable, zero out predictors from itself and any future ones
  for (k in seq_along(idx)) {
    pmat[idx[k], idx[k:length(idx)]] <- 0
  }
  
  pmat
}


### GET PSRF FOR IMPUTED DATASETS ### 
get_psrf <- function(imp_dat, filter = TRUE, threshold = 1.1){
  niter <- imp_dat$iteration
  psrf_mean <- convergence(
    imp_dat,
    diagnostic = "psrf",
    parameter = "mean", 
    it = 1:niter   # Use all iterations
  )
  psrf_sd <- convergence(
    imp_dat,
    diagnostic = "psrf",
    parameter = "sd", 
    it = 1:niter   # Use all iterations
  )
  if (filter) {
    psrf_mean <- psrf_mean |>
      filter(psrf >= threshold & .it == niter)
    psrf_sd <- psrf_sd |>
      filter(psrf >= threshold & .it == niter)
  }
  cat("\nAll PSRF ≥", threshold, "for the mean (at final iteration):\n")
  print(psrf_mean)
  
  cat("\nAll PSRF ≥", threshold, "for the SD (at final iteration):\n")
  print(psrf_sd)
}

# COMPARE TWO PREDICTOR MATRICES ###
compare_pms <- function(a, b, n = 20) {
  idx <- which(a != b, arr.ind = TRUE)
  if (!length(idx)) return(invisible(TRUE))
  diff <- data.frame(
    row = rownames(a)[idx[,1]],
    col = colnames(a)[idx[,2]],
    old = a[idx],
    new = b[idx]
  )
  print(head(diff, n))
  invisible(FALSE)
}


# gformulaMI helpers ----------------------------------------------------------------------------------------------

### EXTRACT GFORMULAMI RESULTS ### 
fit_pool_gform <- function(gdat) {
  po_fit <- with(gdat, lm(frailst ~ factor(regime)))
  results <- as.data.frame(syntheticPool(po_fit)) |>
    mutate(strategy = c("Never treat", "Early exposure", "Late exposure", "Sustained exposure")) |>
    relocate(strategy, .before = Estimate)
  results
}


# Tables and plots ------------------------------------------------------------------------------------------------

### EXPORT RESULTS INTO A WORD TABLE ###
# Helper: convert a single results‐data.frame into a formatted flextable
make_ft <- function(df) {
  # Ensure 'p' exists
  if (!"p" %in% names(df)) {
    stop("Each data frame must contain a column named 'p'.")
  }
  # 1. Format p‐values: "<0.001" if below 0.001, else round to three decimals
  p_orig <- df$p
  df$p <- ifelse(
    p_orig < 0.001,
    "<0.001",
    sprintf("%.3f", round(p_orig, 3))
  )
  
  # 2. Round all other numeric columns to 3 d.p.
  is_num <- sapply(df, is.numeric)
  is_num["p"] <- FALSE
  df[is_num] <- lapply(df[is_num], function(x) round(x, 3))
  
  # 3. Build flextable
  ft <- flextable(df)
  
  # 4. Relabel headers (only for columns that exist)
  header_list <- list(
    strategy = "Strategy",
    Term = "Term",
    Estimate = "Estimate",
    Within = "Var(W)",
    Between = "Var(B)",
    Total = "Var(T)",
    df = "df",
    `95% CI L` = "95% CI Lower",
    `95% CI U` = "95% CI Upper",
    p = "p-value"
  )
  header_list <- header_list[names(header_list) %in% colnames(df)]
  ft <- do.call(set_header_labels, c(list(ft), header_list))
  
  # 5. Force any remaining numeric columns to display three decimals
  num_cols_now <- which(colnames(df) %in% names(df)[sapply(df, is.numeric)])
  if (length(num_cols_now) > 0) {
    ft <- colformat_num(
      ft,
      j = num_cols_now,
      digits = 3,
      na_str = ""
    )
  }
  
  ft
}

# Generic function: write an arbitrary number of results‐tables into one .docx
export_results_doc <- function(results_list, filename) {
  #' results_list: a named list, where each element is a data.frame (or matrix)
  #'               and the name is the desired heading for that table.
  #' filename:     path to the .docx to create (including directory if needed).
  
  # Check that 'results_list' is a list with names
  if (!is.list(results_list) || is.null(names(results_list))) {
    stop("`results_list` must be a named list: names = table headings, values = data.frames.")
  }
  
  # Start a single Word document
  doc <- read_docx()
  heading_fp  <- fp_text(font.size = 12, bold = TRUE)
  heading_par <- fp_par(text.align = "center")
  
  # Iterate over each named element
  for (title_text in names(results_list)) {
    df_raw <- results_list[[title_text]]
    df     <- as.data.frame(df_raw)  # ensure it’s a data.frame
    
    # Create flextable for this table
    ft <- make_ft(df)
    
    # Add a blank line before each heading (except the very first)
    if (length(doc) > 0) {
      doc <- body_add_par(doc, "", style = "Normal")
    }
    
    # Add centered, 12pt‐bold heading
    heading_fpar <- fpar(
      ftext(title_text, prop = heading_fp),
      fp_p = heading_par
    )
    doc <- body_add_fpar(doc, heading_fpar)
    doc <- body_add_par(doc, "", style = "Normal")  # blank line
    
    # Insert the flextable
    doc <- body_add_flextable(doc, value = ft)
  }
  
  # Save the combined .docx
  print(doc, target = filename)
}

### PLOT RESULTS ###
# General plotting function for any number of results objects.
# If exposure_labels isn’t supplied, it defaults to:
#   "Cultural events", "Singing or instruments", "Arts and crafts"
plot_contrasts <- function(
    results_list,
    exposure_labels = c("Cultural events", "Singing or instruments", "Arts and crafts"),
    common_scale = TRUE,
    ylims = c(-1.3, 1.3)
) {
  # Ensure that results_list and exposure_labels have the same length
  if (length(results_list) != length(exposure_labels)) {
    stop("`results_list` and `exposure_labels` must have the same length.")
  }
  
  # Define the common factor levels for “strategy”
  strategy_levels <- c("Never treat", "Early exposure", "Late exposure", "Sustained exposure")
  
  # 1. Process each results data.frame: rename CIs, add exposure label, re-factor strategy
  processed <- lapply(seq_along(results_list), function(i) {
    df <- as.data.frame(results_list[[i]])
    
    df <- df %>%
      rename(
        CI_L = `95% CI L`,
        CI_U = `95% CI U`
      ) %>%
      mutate(
        exposure = exposure_labels[i],
        strategy = factor(strategy, levels = strategy_levels)
      )
    
    df
  })
  
  # 2. Combine all and drop the “Never treat” reference group
  all_results <- bind_rows(processed) %>%
    filter(strategy != "Never treat") %>%
    mutate(exposure = factor(exposure, levels = exposure_labels))
  
  # 3. Build the ggplot
  p <- ggplot(all_results, aes(x = exposure, y = Estimate, colour = strategy)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_errorbar(
      aes(ymin = CI_L, ymax = CI_U),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    labs(
      x      = "Exposure",
      y      = "Mean difference in frailty (95% CI)",
      colour = "Exposure strategy"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.margin = margin(t = 5.5, r = 20, b = 5.5, l = 5.5)
    )
  if (common_scale) {
    p <- p + ylim(ylims[1], ylims[2])
  }
  
  return(p)
}
