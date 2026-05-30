# ==============================================================================
# State Capacity x Fragility Distribution Table
# V-Dem and WGI classifications, both including rule of law
# ==============================================================================
# Author: Cedric Antunes (FGV-CEPESP) ------------------------------------------
# Date: May, 2026 --------------------------------------------------------------

# Assumes standardized_data_final is loaded 

# Required packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(knitr)
})

# ------------------------------------------------------------------------------
# Required variables -----------------------------------------------------------
# ------------------------------------------------------------------------------
required_vars <- c(
  "COUNTRY_NAME",
  "ISO_CODE_3",
  "YEAR",
  "STATE_CAPACITY",
  "VDEM_FRAGILE_IDEAL",
  "VDEM_NON_FRAGILE_IDEAL",
  "WGI_FRAGILE_IDEAL",
  "WGI_LOADING_FACTOR_1_NORMALIZED"
)

# ------------------------------------------------------------------------------
# Helpers ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Re-orienting fragility index: the higher, the more fragile
orient_higher_quality <- function(x, fragile_binary) {
  cc <- complete.cases(x, fragile_binary)
  
  if (sum(cc) < 10) return(x)
  
  r <- suppressWarnings(
    cor(x[cc], fragile_binary[cc], use = "complete.obs")
  )
  
  if (is.na(r)) return(x)
  
  # If the WGI index is positively correlated with fragility,
  # reverse it so higher = better state quality.
  if (r > 0) -x else x
}

# ------------------------------------------------------------------------------
# Preparing the data -----------------------------------------------------------
# ------------------------------------------------------------------------------
df_base <- standardized_data_final |>
  ungroup() |>
  transmute(
    country = COUNTRY_NAME,
    iso3 = as.character(ISO_CODE_3),
    year = as.integer(YEAR),
    
    state_capacity = STATE_CAPACITY,
    
    # V-Dem classification, preferred measure including rule of law
    vdem_fragile = VDEM_FRAGILE_IDEAL,
    vdem_non_fragile = VDEM_NON_FRAGILE_IDEAL,
    
    # WGI classification, preferred measure including rule of law
    wgi_fragile = WGI_FRAGILE_IDEAL,
    wgi_raw = WGI_LOADING_FACTOR_1_NORMALIZED
  ) |>
  mutate(
    wgi_quality = orient_higher_quality(wgi_raw, wgi_fragile)
  )

# ------------------------------------------------------------------------------
# Constructing V-Dem and WGI status variables ----------------------------------
# ------------------------------------------------------------------------------
df_status <- df_base |>
  group_by(year) |>
  mutate(
    # WGI quality percentile is only meaningful where WGI is observed.
    wgi_quality_pctl = if_else(
      !is.na(wgi_quality),
      percent_rank(wgi_quality),
      NA_real_
    )
  ) |>
  ungroup() |>
  mutate(
    vdem_status = case_when(
      vdem_fragile == 1 ~ "Fragile",
      vdem_non_fragile == 1 ~ "Non-fragile",
      !is.na(vdem_fragile) | !is.na(vdem_non_fragile) ~ "Transitioning / other",
      TRUE ~ NA_character_
    ),
    
    # WGI non-fragile is defined as above the 35th percentile
    # of the yearly WGI quality distribution when WGI is observed.
    wgi_status = case_when(
      is.na(wgi_quality) & is.na(wgi_fragile) ~ NA_character_,
      wgi_fragile == 1 ~ "Fragile",
      wgi_quality_pctl >= 0.35 ~ "Non-fragile",
      !is.na(wgi_quality_pctl) ~ "Transitioning / other",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# Classifying capacity-fragility combinations ----------------------------------
# ------------------------------------------------------------------------------
make_capacity_fragility_table <- function(data, status_var, source_label) {
  
  data %>%
    filter(
      !is.na(state_capacity),
      !is.na(.data[[status_var]])
    ) %>%
    group_by(year) %>%
    mutate(
      capacity_pctl = percent_rank(state_capacity)
    ) %>%
    ungroup() %>%
    mutate(
      category = case_when(
        capacity_pctl >= 0.70 & .data[[status_var]] == "Fragile" ~
          "High capacity but fragile",
        
        capacity_pctl <= 0.30 & .data[[status_var]] == "Non-fragile" ~
          "Low capacity but non-fragile",
        
        capacity_pctl >= 0.70 & .data[[status_var]] == "Non-fragile" ~
          "High capacity and non-fragile",
        
        capacity_pctl <= 0.30 & .data[[status_var]] == "Fragile" ~
          "Low capacity and fragile",
        
        TRUE ~ "Middle / mixed"
      )
    ) %>%
    count(category, name = "country_years") %>%
    mutate(
      source = source_label,
      share = country_years / sum(country_years)
    )
}

# ------------------------------------------------------------------------------
# Building table for V-Dem and WGI ---------------------------------------------
# ------------------------------------------------------------------------------
table_vdem <- make_capacity_fragility_table(
  data = df_status,
  status_var = "vdem_status",
  source_label = "VDEM"
)

table_wgi <- make_capacity_fragility_table(
  data = df_status,
  status_var = "wgi_status",
  source_label = "WGI"
)

capacity_fragility_distribution <- bind_rows(
  table_vdem,
  table_wgi
) |>
  mutate(
    category = factor(
      category,
      levels = c(
        "High capacity and non-fragile",
        "Low capacity and fragile",
        "Low capacity but non-fragile",
        "High capacity but fragile",
        "Middle / mixed"
      )
    ),
    share_print = percent(share, accuracy = 0.1)
  ) |>
  arrange(category, source)

print(capacity_fragility_distribution)

# ------------------------------------------------------------------------------
# Wide version for paper appendix ----------------------------------------------
# ------------------------------------------------------------------------------
capacity_fragility_distribution_wide <- capacity_fragility_distribution |>
  select(source, category, country_years, share_print) |>
  pivot_wider(
    names_from = source,
    values_from = c(country_years, share_print),
    names_glue = "{source}_{.value}"
  ) |>
  arrange(category)

print(capacity_fragility_distribution_wide)

# ------------------------------------------------------------------------------
# LaTeX-ready table ------------------------------------------------------------
# ------------------------------------------------------------------------------
latex_table <- capacity_fragility_distribution_wide |>
  mutate(
    category = as.character(category)
  ) |>
  select(
    Category = category,
    `V-Dem country-years` = VDEM_country_years,
    `V-Dem share` = VDEM_share_print,
    `WGI country-years` = WGI_country_years,
    `WGI share` = WGI_share_print
  )

print(latex_table)

knitr::kable(
  latex_table,
  format = "latex",
  booktabs = TRUE,
  caption = "State Capacity and Fragility: High/Low Country-Year Distribution",
  label = "tab:capacity_fragility_quadrants"
)
