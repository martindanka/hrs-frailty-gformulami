
# Packages --------------------------------------------------------------------------------------------------------

library(tidyverse)
library(haven)

# Shared functions used across scripts ---------------------------------------------------------------------------
## Before sourcing, ensure your working directory is set to the project root (i.e. the HRS frailty folder)
## Example: setwd("/Users/yourname/Documents/HRS Frailty")
source("code/helpers.R")

# Functions -------------------------------------------------------------------------------------------------------
zap_stata <- function(data) {
  data |>
    zap_labels() |>
    zap_label() |>
    zap_formats()
}

# Import data -----------------------------------------------------------------------------------------------------

covars_long_raw <- read_dta("data/cams_covars00-18.dta")
frailty_long_raw <- read_dta("data/frailty_full.dta")
cams_wide_raw <- read_dta("data/cams05_17.dta")
tracker <- read_dta("data/trk2022tr_r.dta")
  
# Deaths ----------------------------------------------------------------------------------------------------------

# Derive death indicator and flag participants who died during follow-up
# (reproduces logic from the pre-processing pipeline)
covars_long_alive <- covars_long_raw %>%
  ## 1.  Sort by person and wave, add a mid-month “day” and build full dates
  arrange(wave) %>% 
  mutate(day = 15,
         iwdate = make_date(iwyear, iwmonth, day),
         ddate = make_date(knowndeceasedyr, knowndeceasedmo, day)) %>%
  
  ## 2.  Initial death flag: interview on/after death  →  1
  mutate(
    death = case_when(
      # 1) interview date exists AND interview ≥ death → dead
      !is.na(iwdate) & iwdate >= ddate ~ 1L,
      # 2) interview date exists (but interview < death, or no death date) → alive
      !is.na(iwdate) ~ 0L,
      # 3) otherwise → leave missing
      .default = NA
    ),
    
    ## 3.  Use alive codes where death is still missing
    death = case_when(
      is.na(death) & alive %in% 5:6 ~ 1, 
      is.na(death) & alive %in% 1:2 ~ 0,
      .default = death
    )
  ) %>% 
  
  ## 4.  Carry a “dead” flag forward within each person
  group_by(pid) %>% 
  arrange(wave, .by_group = TRUE) %>% 
  mutate(
    death = case_when(
      # 1) if still missing but last wave was dead, stay dead
      is.na(death) & cummax(replace_na(death, 0)) == 1 ~ 1L,
      # 2) if still missing and no death date, assume alive
      is.na(death) & is.na(ddate) ~ 0L,
      # 3) otherwise leave whatever death already is (0 or 1)
      TRUE ~ death
    )
  ) %>% 
  ungroup() %>% 
  
  ## 5.  Keep only the FIRST wave where death==1 for any person
  group_by(pid, death) %>% 
  mutate(n = row_number()) %>% 
  ungroup() %>% 
  filter(!(death == 1 & n > 1)) %>% 
  select(-n) %>% 
  
  ## 6.  Drop every participant who ever died during follow-up
  group_by(pid) %>% 
  mutate(mk = max(death, na.rm = TRUE)) %>% 
  ungroup() %>%

  ## 7.  Tidy up
  select(-day, -iwdate, -ddate)

# Consistency check against pre-processing script (Feifei Bu's Stata code)
# t <- covars_long_alive |>
#   filter(mk == 0)
# table(covars_long_alive$death)
# Gives the same counts.


# Additional variables --------------------------------------------------------------------------------------------

# Frailty
frailty_long_sel <- frailty_long_raw |>
  dplyr::select(pid, wave, frail) 

data_long_full <- covars_long_alive |>
  left_join(frailty_long_sel, by = c("pid", "wave"))

# Age
## Define main wave years for deterministic age calculation
wave_years <- c(`5` = 2001, `6` = 2003, `7` = 2005, `8` = 2007,
                `9` = 2009, `10` = 2011, `11` = 2013, `12` = 2015,
                `13` = 2017, `14` = 2019)

## Derive age
data_long_full <- data_long_full %>%
  mutate(
    age_d = case_when(
      # — if age is observed (i.e. not the 999 “missing” code) 
      age != 999 ~ age,
      # — if age==999 *and* interview year is known 
      age == 999 & !is.na(iwyear) ~ iwyear - birthyr,
      # — otherwise (age==999 & iwyear is NA) 
      age == 999 & is.na(iwyear) ~ wave_years[as.character(wave)] - birthyr
    )
  )

# Data cleaning ---------------------------------------------------------------------------------------------------

# Select variables & recode factors
data_long_clean_step1 <- data_long_full |>
  # Select variables
  dplyr::select(pid, wave, schlyrs, ethnic, male, usborn, mari, employ,
                hhown, nhsafe, hincome, frail, age_d, shlt, psyche, conde, cog27, moneya, 
                hh2nd, mk, birthyr) |>
  mutate(
    # Recode unordered factors
    across(
      c(male, ethnic, starts_with("art"), usborn, mari, employ, hhown, moneya), 
      ~as_factor(.x)
    ),
    # Recode ordinal variables
    across(
      c(nhsafe, hincome, shlt),
      ~factor(.x, ordered = TRUE)
    )
  ) |>
  # Zap Stata labels
  zap_stata() |>
  # Rename arts variables and cognition
  rename(
    cog = cog27
  )

# One-off fixes
data_long_clean_step2 <- data_long_clean_step1 |>
  mutate(
    psyche = as_factor(psyche),
    moneya = case_when(
      moneya == "0.no" ~ "no",
      moneya == "1.yes" ~ "yes",
      .default = NA
    ),
    moneya = as_factor(moneya),
    age_d = as.integer(age_d), 
    hh2nd = case_when(
      hh2nd == 0 ~ "no",
      hh2nd == 1 ~ "yes",
      .default = NA
    ),
    hh2nd = as_factor(hh2nd)
  )


# Sort data based on time ordering
## This will be needed for gFormulaMI
data_long_clean_step3 <- data_long_clean_step2  |>
  relocate(pid, wave, mk, birthyr, schlyrs, male, ethnic, usborn, age_d, mari, employ, hhown, hincome, nhsafe, conde, psyche, shlt,
           frail, cog, moneya, hh2nd)

data_long_clean <- data_long_clean_step3 %>%
  group_by(pid) %>%
  mutate(
    ## grab the value once per person
    frailst = frail[wave == 14][1]
    ) %>%
  ungroup() 

rm(data_long_clean_step1, data_long_clean_step2, data_long_clean_step3)

# Convert into wide format
data_wide <- data_long_clean |>
  filter(wave != 14) |>
  pivot_wider(
    id_cols = c(pid, mk, birthyr, schlyrs, male, ethnic, usborn, frailst),
    names_from = wave,
    names_sep = "",
    values_from = c(age_d, mari, employ, hhown, hincome, nhsafe, conde, psyche, shlt, frail, cog, moneya, hh2nd)
  )

# Merge with exposure ---------------------------------------------------------------------------------------------

## ------------------------------------------------------------------
## 1.  Wave-number lookup:  "05" → "7",  "07" → "8", … "17" → "13"
## ------------------------------------------------------------------
wave_lookup <- setNames(
  as.character(7:13),
  sprintf("%02d", seq(5, 17, 2))
)

## ------------------------------------------------------------------
## 2.  Rename art- and qtype-variables
## ------------------------------------------------------------------
cams_wide <- cams_wide_raw %>% 
  
  ## ---- a)  art1_/art2_/art3_  -------------------------------------
rename_with(
  .cols = matches("^art[123]_"),
  .fn = function(v) {
    v %>%
      str_replace("^art1_", "cultural_events") %>%
      str_replace("^art2_", "music_activity") %>%
      str_replace("^art3_", "craft_hobbies") %>%
      str_replace_all(wave_lookup)
  }
) %>%

  ## ---- b)  qtypeXX  ----------------------------------------------
rename_with(
  .cols = matches("^qtype\\d{2}$"),
  .fn = ~ str_replace_all(.x, wave_lookup)
)

cams_wide_sel <- cams_wide |>
  select(pid, starts_with(c("cultural_events", "music_activity", "craft_hobbies"))) |>
  mutate(across(starts_with(c("cultural_events", "music_activity", "craft_hobbies")), 
                ~ as_factor(.x))) |>
  zap_stata()

data_wide <- data_wide |>
  left_join(cams_wide_sel, by = "pid")


# Exclusions ------------------------------------------------------------------------------------------------------

# Determine who participated in CAMS W1, W2, or W3 
tracker <- tracker |> 
  # Generate unique id
  mutate(pid = paste0(HHID, PN)) |>
  relocate(pid, .after = PN)

table(tracker$CAMS01)

tracker <- tracker |>
  mutate(
    # Participated in CAMS W1, W2, or W3?
    eligible = case_when(
      CAMS01 == 1 | CAMS03 == 1 | CAMS05 %in% c(1, 2) ~ 1,
      .default = 0
    ),
    # Eligible to participate in CAMS W1, W2, or W3?
    eligible_broad = case_when(
      CAMS01 %in% c(1, 5) | CAMS03 %in% c(1, 5) | CAMS05 %in% c(1, 2, 5, 6) ~ 1,
      .default = 0
    )
  ) |>
  relocate(eligible, eligible_broad, .after = pid)

table(tracker$eligible)
table(tracker$eligible_broad)

id <- tracker |>
  filter(eligible == 1) |>
  select(pid)

id_broad <- tracker |>
  filter(eligible_broad == 1) |>
  select(pid)

# Check if all eligible are already included
all(id_broad$pid %in% data_wide$pid)
all(id$pid %in% data_wide$pid)
miss_id <- setdiff(id$pid, data_wide$pid)
miss_id

# Four participants not included.

# Prepare those eligible
## Defined by having participated in either CAMS 2001, CAMS 2003 or CAMS 2005.
data_wide_elig <- data_wide |>
  filter(pid %in% id$pid)
length(unique(data_wide_elig$pid))
## 6,880

# Exclude those who died before Wave 14
data_wide_alive <- data_wide_elig |>
  filter(mk != 1) |>
  select(-mk)

length(unique(data_wide_alive$pid))
## 4,042

# Exclude those below the age of 50 at baseline (CAMS 2005)
data_wide_sample <- data_wide_alive |>
  filter(birthyr < 1954) |>
  select(-birthyr)
length(unique(data_wide_sample$pid))
## 3,775

# Derive info on missingness
## Missing outcome
data_wide_compl_outcome <- data_wide_sample |>
  filter(!is.na(frailst))

## Missing time-fixed variables
data_wide_compl_base <- data_wide_compl_outcome |>
  filter(!is.na(schlyrs) & !is.na(male) & !is.na(ethnic) & !is.na(usborn))

## Missing time-varying variables (used in substantive models)
tvc_bases <- c("mari", "employ", "hhown", "hincome", "nhsafe", "frail", "cultural_events", "music_activity", "craft_hobbies")
time_idx <- 7:13
tvc_names <- as.vector(outer(tvc_bases, time_idx, paste0))
data_wide_compl_timevar <- data_wide_compl_base %>%
  filter(if_all(all_of(tvc_names), ~ !is.na(.)))

## Get numbers for flowchart
length(unique(data_wide_compl_outcome$pid)) # 2753
length(unique(data_wide_compl_base$pid)) # 2750
length(unique(data_wide_compl_timevar$pid)) # 1672

# Derive analytical dataset ---------------------------------------------------------------------------------------

# Reordering columns ----------------------------------------------------------------------------------------------
## gFormulaMI uses the column order to determine the correct timing.
## Therefore, the columns must be sorted. 

# 1) Identify the time-fixed variables.
fixed_vars <- c("pid", "schlyrs", "male", "ethnic", "usborn", "frailst")

# 2) Identify time-varying columns 
time_var_cols <- setdiff(colnames(data_wide_sample), fixed_vars)

# 3) Parse each time‐varying column:
#    - 'base' = variable name (e.g. "mari", "nhsafe", "cultural_events")
#    - 'time' = numeric index from the column name (e.g. "5", "6", "7", ...)
#    - 'rank_treat' is 1 for exposure variables, which we want to appear last within each time block.
parsed <- tibble(col = time_var_cols) %>%
  mutate(
    # Extract everything up to (but not including) the digits for 'base'
    base = str_extract(col, "^.*?(?=\\d)"),
    # Extract the trailing digits as the 'time'
    time = str_extract(col, "\\d+$"),
    time = as.integer(time),
    
    # Put exposure variables last within each time block
    rank_treat = case_when(
      base %in% c("cultural_events", "music_activity", "craft_hobbies") ~ 1,
      .default = 0
    )
  )

# 4) Sort the time‐varying columns: 
#    a) by the time index ascending
#    b) within each time index, exposure variables last
#    c) otherwise, keep the original order (default behaviour)
parsed_ordered <- parsed %>%
  arrange(time, rank_treat)

# 5) Build the new column order 
#    - all baseline columns first,
#    - then the newly ordered time‐varying columns.
new_col_order <- c(fixed_vars, parsed_ordered$col)

# 6) Rearrange 
data_wide_clean <- data_wide_sample %>%
  dplyr::select(all_of(new_col_order)) |>
  relocate(c(frailst), .after = last_col()) |>
  relocate(male, ethnic, usborn, .before = schlyrs) 

# Check the new column order
str(data_wide_clean)

# Check missing data
count_missing(data_wide_clean)

# Remove unused age
data_wide_clean <- data_wide_clean |>
  select(-c(age_d5, age_d6, age_d8, age_d9, age_d10, age_d11, age_d12, age_d13))

# Export
saveRDS(data_wide_clean, "data/data_wide_clean.rds")
