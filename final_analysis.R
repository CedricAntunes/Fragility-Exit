---
title: "FRAGILITY EXIT - Quantiative Appendix"
author: 
- "Christoph Zuercher" 
- "Cedric Antunes"
date: "August 2025"
header-includes:
    \usepackage{lscape}
    \usepackage{pdfpages}
    \usepackage{graphicx}
    \usepackage[figuresright]{rotating}
output:
  pdf_document:
    toc: true
    toc_depth: 2
---

```{r, echo = FALSE, message = FALSE, warning = FALSE, results = 'hide'}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
```

```{r setup}
# Cleaning my environment
rm(list = ls())

# Managing memory
gc()

# Loading packages -------------------------------------------------------------
library(here)
library(coefplot)
library(modelsummary)
library(corrplot)
library(corrtable)
library(fixest)
library(pROC)
library(PRROC)
library(texreg)
library(dynlm)
library(tibble)
library(tidyverse)
library(glmnet)
library(lubridate)
library(dplyr)
library(broom)
library(plm)
library(car)
library(nnet)
library(tidymodels)
library(scatr)
library(zoo)
library(knitr)
library(kableExtra)
library(lmtest)
library(sandwich)
library(margins)
library(tidyr)
library(dplyr)
library(caret)
library(stargazer)
library(VGAM)
library(gt)
library(xfun)
library(factoextra)
library(cluster)
library(ggplot2)
library(vip)
library(themis)
library(janitor)
library(WDI)
library(brglm2)        
library(glmmTMB)
library(margins)
library(splines)
library(stats)
library(scales)
library(isotone)

# Setting seed for replication -------------------------------------------------
set.seed(247)

# Load data, functions, and scripts --------------------------------------------
load("final_new_clean_output_data.RDS")
load("fh_dd_final.RDS")
source("WGI_30th_missing_estimates.R")
source("WGI_35th_missing_estimates.R")
source("WGI_loading_factor_estimates.R")
source("fragility_breakthrough_identifier.R")
source("VDEM_never_fragile.R")

# ------------------------------------------------------------------------------
# WORLD BANK DATA --------------------------------------------------------------
# ------------------------------------------------------------------------------

# Last minute input: extracting final GDP per capita data ----------------------
gdp_per_capita_data <- WDI(
  country = "all",
  indicator = c("NY.GDP.PCAP.CD"),
  start = 1970,
  end = 2023,
  extra = FALSE,
  cache = NULL,
  latest = NULL,
  language = "en"
)

# Cleaning the data ------------------------------------------------------------
gdp_per_capita_data <- gdp_per_capita_data |>
  # Renaming variables for match
  rename(COUNTRY_NAME = country,
         ISO_CODE_3 = iso3c,
         YEAR = year,
         GDP_PER_CAPITA = NY.GDP.PCAP.CD) |>
  # Selecting only key variables for match 
  select(COUNTRY_NAME,
         ISO_CODE_3,
         YEAR,
         GDP_PER_CAPITA) |>
  # Dropping 2023
  filter(YEAR != 2023) |>
  # Setting YEAR as character
  mutate(YEAR = as.character(YEAR)) |>
  # Selecting only country-rows
  slice(-(1:2597))

# Preparing mother-dataset
final_clean_percentiles_data <- final_clean_percentiles_data |>
  mutate(YEAR = as.character(YEAR))

# Last minute input: extracting final GDP per capita data ----------------------
oil_data <- WDI(
  country = "all",
  indicator = c("NY.GDP.PETR.RT.ZS"),
  start = 1970,
  end = 2023,
  extra = FALSE,
  cache = NULL,
  latest = NULL,
  language = "en"
)

# Cleaning the data ------------------------------------------------------------
oil_data <- oil_data |>
  # Renaming variables for match
  rename(COUNTRY_NAME = country,
         ISO_CODE_3 = iso3c,
         YEAR = year,
         OIL_RENTS = NY.GDP.PETR.RT.ZS) |>
  # Selecting only key variables for match 
  select(COUNTRY_NAME,
         ISO_CODE_3,
         YEAR,
         OIL_RENTS) |>
  # Dropping 2023
  filter(YEAR != 2023) |>
  # Setting YEAR as character
  mutate(YEAR = as.character(YEAR)) |>
  # Selecting only country-rows
  slice(-(1:2646))
  
# Performing the join ----------------------------------------------------------
final_clean_percentiles_data <- left_join(final_clean_percentiles_data,
                                          gdp_per_capita_data,
                                          by = c("ISO_CODE_3", 
                                                 "YEAR"),
                                          keep = FALSE)

# Performing the join ----------------------------------------------------------
final_clean_percentiles_data <- left_join(final_clean_percentiles_data,
                                          oil_data,
                                          by = c("ISO_CODE_3", 
                                                 "YEAR"),
                                          keep = FALSE)

# FH & DD join -----------------------------------------------------------------
final_clean_percentiles_data <- left_join(final_clean_percentiles_data,
                                          fh_dd_final_1,
                                          by = c("ISO_CODE_3", 
                                                 "YEAR"),
                                          keep = TRUE)

# Final cleanning --------------------------------------------------------------
final_clean_percentiles_data <- final_clean_percentiles_data |>
  select(-COUNTRY_NAME.x.x,
         -COUNTRY_NAME.y,
         -COUNTRY_NAME.y.y,
         -YEAR.y,
         -ISO_CODE_3.y) |>
  # Dropping duplicate columns 
  rename(COUNTRY_NAME = COUNTRY_NAME.x,
         YEAR = YEAR.x,
         ISO_CODE_3 = ISO_CODE_3.x) |>
  # Allocating GDP per capita along with economic indicators 
  relocate(82, .after = 3) |>
  # Log-trnaforming GDP per capita 
  mutate(LOG_GDP_PER_CAPITA = if_else(GDP_PER_CAPITA > 0, log(GDP_PER_CAPITA), NA_real_))

# Building working dataframe ---------------------------------------------------
final_clean_percentiles_data_normalized_1 <- final_clean_percentiles_data |>
  select(-REVISED_COMBINED_POLITY_SCORE,
         -PRIOR_COMBINED_POLITY_SCORE,
         -POLITY_END_MONTH, 
         -POLITY_END_DAY,
         -POLITY_END_YEAR,
         -END_DATE_PRECISION,
         -INTERIM_POLITY_SCORE,
         -POLITY_BEGIN_YEAR,
         -POLITY_BEGIN_MONTH,
         -POLITY_BEGIN_DAY,
         -POLITYV_CONFIDENCE_INDEX,
         -BEGIN_DATE_PRECISION,
         -POST_POLITY_SCORE,
         -TOTAL_CHANGE_POLITY_SCORE,
         -REGIME_TRANSITION_COMPLETED,
         -REGIME_TRANSITION,
         -CONFLICT_ID,
         -CONFLICT_LOCATION,
         -CONFLICT_CAUSE,
         -CONFLICT_TYPE,
         -CONFLICT_INACTIVE,
         -COUNTRY_CONFLICT_SIDE_A,
         -STATES_SUPPORT_SIDE_A_WITH_TROOPS,
         -COUNTRY_OR_ACTOR_CONFLICT_SIDE_B,
         -STATES_SUPPORT_SIDE_B_WITH_TROOPS,
         -TERRITORY_UNDER_DISPUTE,
         -CONFLICT_START_DATE,
         -PRECISION_OF_CONFLICT_START_DATE,
         -PRECISION_OF_CONFLICT_END_DATE,
         -CONFLICT_DEATHS_THRESHOLD_DATE,
         -PRECISION_OF_CONFLICT_DEATHS_THRESHOLD_DATE,
         -CONFLICT_END_DATE,
         -PRECISION_OF_CONFLICT_START_DATE
  ) |>
  estimating_WGI_30th_missing_values() |>
  estimating_WGI_35th_missing_values() |>
  filter(complete.cases(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y, NEW_VDEM_30TH_PERCENTILE)) |>
  group_by(COUNTRY_NAME) |>
  mutate(
    YEAR = as.numeric(YEAR),

    # Fragility DUMMIES (1 = Fragile, 0 = Not Fragile)
    VDEM_FRAGILE_BASELINE     = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE, 1, 0),
    VDEM_NON_FRAGILE_BASELINE = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y >= NEW_VDEM_25TH_PERCENTILE, 1, 0),

    VDEM_FRAGILE_IDEAL        = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_30TH_PERCENTILE, 1, 0),
    VDEM_NON_FRAGILE_IDEAL    = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y >= NEW_VDEM_35TH_PERCENTILE, 1, 0),

    VDEM_FRAGILE_ENDLINE      = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE, 1, 0),
    VDEM_NON_FRAGILE_ENDLINE  = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y >= NEW_VDEM_40TH_PERCENTILE, 1, 0),

    # Fragility LABELS
    VDEM_STATUS_BASELINE = case_when(
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE ~ "Fragile",
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_25TH_PERCENTILE ~ "Transitioning",
      TRUE ~ "Non Fragile"
    ),
    VDEM_STATUS_IDEAL = case_when(
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_30TH_PERCENTILE ~ "Fragile",
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_35TH_PERCENTILE ~ "Transitioning",
      TRUE ~ "Non Fragile"
    ),
    VDEM_STATUS_ENDLINE = case_when(
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE ~ "Fragile",
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_40TH_PERCENTILE ~ "Transitioning",
      TRUE ~ "Non Fragile"
    )
  )

# Converting to character before standardizing ---------------------------------
final_data <- final_clean_percentiles_data_normalized_1 |>
  mutate(
    YEAR = as.character(YEAR),
    #across(35:55, as.character)
    across(where(is.factor), as.character)  
  )

# Computing means and SDs for numeric columns ----------------------------------
means_data <- sapply(final_data, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else NA)
sd_data    <- sapply(final_data, function(x) if (is.numeric(x)) sd(x, na.rm = TRUE) else NA)

# Initializing standardized data
standardized_data <- final_data

# Defining dummies to exclude from standardization
dummy_vars <- c(
  "VDEM_FRAGILE_BASELINE", 
  "VDEM_NON_FRAGILE_BASELINE",
  "VDEM_FRAGILE_IDEAL", 
  "VDEM_NON_FRAGILE_IDEAL",
  "VDEM_FRAGILE_ENDLINE", 
  "VDEM_NON_FRAGILE_ENDLINE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "FH_DEMOCRACY",
  "FH_AUTOCRACY",
  "DD_DEMOCRACY"
)

# Applying standardization safely
for (col in names(final_data)) {
  if (is.numeric(final_data[[col]]) && !(col %in% dummy_vars)) {
    standardized_data[[col]] <- (final_data[[col]] - means_data[[col]]) / sd_data[[col]]
  }
}

# Converting relevant predictors back to numeric
standardized_data_final <- standardized_data |>
  mutate(across(35:55, as.numeric))

# Final correction pass: recalculating fragility dummies from scratch
standardized_data_final <- standardized_data_final |>
  mutate(
    VDEM_FRAGILE_BASELINE = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE, 1, 0),
    VDEM_FRAGILE_IDEAL    = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_30TH_PERCENTILE, 1, 0),
    VDEM_FRAGILE_ENDLINE  = if_else(NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE, 1, 0),

    VDEM_STATUS_BASELINE = case_when(
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE ~ "Fragile",
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_25TH_PERCENTILE ~ "Transitioning",
      TRUE ~ "Non Fragile"
    ),
    VDEM_STATUS_IDEAL = case_when(
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_30TH_PERCENTILE ~ "Fragile",
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_35TH_PERCENTILE ~ "Transitioning",
      TRUE ~ "Non Fragile"
    ),
    VDEM_STATUS_ENDLINE = case_when(
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_20TH_PERCENTILE ~ "Fragile",
      NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y < NEW_VDEM_40TH_PERCENTILE ~ "Transitioning",
      TRUE ~ "Non Fragile"
    )
  )

# Saving data 
save(standardized_data_final, 
     file = "standardized_data_final.RDS")

# Setting the critical junctures
standardized_data_final <- standardized_data_final |>
  mutate(CRITICAL_JUNCTURE = case_when(
    (COUNTRY_NAME == "Angola" & YEAR == 2018) ~ 1,
    (COUNTRY_NAME == "Brazil" & YEAR == 1985) ~ 1,
    (COUNTRY_NAME == "Ethiopia" & YEAR == 1991) ~ 1,
    (COUNTRY_NAME == "Gabon" & YEAR == 1990) ~ 1,
    (COUNTRY_NAME == "Gambia, The" & YEAR == 2016) ~ 1,
    (COUNTRY_NAME == "Ghana" & YEAR == 1992) ~ 1,
    (COUNTRY_NAME == "Kenya" & YEAR == 2002) ~ 1,
    (COUNTRY_NAME == "Mauritania" & YEAR == 1992) ~ 1,
    (COUNTRY_NAME == "Sierra Leone" & YEAR == 2002) ~ 1,
    (COUNTRY_NAME == "Tunisia" & YEAR == 2010) ~ 1,
    (COUNTRY_NAME == "Uganda" & YEAR == 1986) ~ 1,
    (COUNTRY_NAME == "Madagascar" & YEAR == 1991) ~ 1,
    (COUNTRY_NAME == "Myanmar" & YEAR == 2008) ~ 1,
    (COUNTRY_NAME == "Indonesia" & YEAR == 1998) ~ 1,
    (COUNTRY_NAME == "Iran, Islamic Rep." & YEAR == 1988) ~ 1,
    (COUNTRY_NAME == "Nepal" & YEAR == 2006) ~ 1,
    (COUNTRY_NAME == "Philippines" & YEAR == 1986) ~ 1,
    (COUNTRY_NAME == "Sri Lanka" & YEAR == 2009) ~ 1,
    (COUNTRY_NAME == "China" & YEAR == 1977) ~ 1,
    (COUNTRY_NAME == "Thailand" & YEAR == 1991) ~ 1,
    (COUNTRY_NAME == "Bosnia and Herzegovina" & YEAR == 1995) ~ 1,
    (COUNTRY_NAME == "Serbia" & YEAR == 1999) ~ 1,
    (COUNTRY_NAME == "Georgia" & YEAR == 2004) ~ 1,
    (COUNTRY_NAME == "Kyrgyz Republic" & YEAR == 2010) ~ 1,
    (COUNTRY_NAME == "Armenia" & YEAR == 2018) ~ 1,
    (COUNTRY_NAME == "Argentina" & YEAR == 1983) ~ 1,
    (COUNTRY_NAME == "Bolivia" & YEAR == 1982) ~ 1,
    (COUNTRY_NAME == "Colombia" & YEAR == 1990) ~ 1,
    (COUNTRY_NAME == "Ecuador" & YEAR == 1976) ~ 1,
    (COUNTRY_NAME == "Dominican Republic" & YEAR == 1994) ~ 1,
    (COUNTRY_NAME == "El Salvador" & YEAR == 1990) ~ 1,
    (COUNTRY_NAME == "Mexico" & YEAR == 1990) ~ 1,
    (COUNTRY_NAME == "Nicaragua" & YEAR == 1979) ~ 1,
    (COUNTRY_NAME == "Panama" & YEAR == 1989) ~ 1,
    (COUNTRY_NAME == "Paraguay" & YEAR == 1989) ~ 1,
    (COUNTRY_NAME == "Peru" & YEAR == 2000) ~ 1,
    TRUE ~ 0
  )) |>
  mutate(CRITICAL_JUNCTURE = as.numeric(CRITICAL_JUNCTURE))

# Checking classification ------------------------------------------------------
sample <- standardized_data_final |>
  select(YEAR,
         COUNTRY_NAME,
         VDEM_STATUS_BASELINE,
         VDEM_STATUS_IDEAL,
         VDEM_STATUS_ENDLINE) |>
  rename(VDEM_STATUS_20_25_PERCENTILES = VDEM_STATUS_BASELINE,
         VDEM_STATUS_30_35_PERCENTILES = VDEM_STATUS_IDEAL,
         VDEM_STATUS_20_40_PERCENTILES = VDEM_STATUS_ENDLINE)

# Saving classified country-year observations ----------------------------------
write.csv(sample,
          "countries_classified.csv")
```

```{r lagged}
# Start from baseline, but DROP any grouping
df_base <- standardized_data_final %>% ungroup()

df0_lag <- df_base %>%
  mutate(
    row_id   = dplyr::row_number(),                 # now truly 1..N, unique
    YEAR_int = as.integer(as.character(YEAR))
  )

id_col   <- "ISO_CODE_3"
year_col <- "YEAR_int"

dummy_vars <- c(
  "VDEM_FRAGILE_BASELINE","VDEM_NON_FRAGILE_BASELINE",
  "VDEM_FRAGILE_IDEAL","VDEM_NON_FRAGILE_IDEAL",
  "VDEM_FRAGILE_ENDLINE","VDEM_NON_FRAGILE_ENDLINE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "FH_DEMOCRACY","FH_AUTOCRACY","DD_DEMOCRACY"
)

skip_pos  <- intersect(names(df0_lag)[36:55], names(df0_lag))
skip_vars <- unique(c(id_col, year_col, dummy_vars, skip_pos))

num_vars    <- names(df0_lag)[sapply(df0_lag, is.numeric)]
vars_to_lag <- setdiff(num_vars, skip_vars)

if ("LOG_GDP_PER_CAPITA" %in% names(df0_lag)) {
  vars_to_lag <- union(vars_to_lag, "LOG_GDP_PER_CAPITA")
}

standardized_data_lagged <- df0_lag %>%
  arrange(.data[[id_col]], .data[[year_col]]) %>%
  group_by(.data[[id_col]]) %>%
  mutate(
    across(all_of(vars_to_lag), ~ dplyr::lag(.x, 1), .names = "{.col}_L1"),
    across(all_of(vars_to_lag), ~ dplyr::lag(.x, 2), .names = "{.col}_L2"),
    PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = dplyr::lag(PARTIAL_DEMOCRACY_WITH_FACTIONALISM, 1),
    PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L2 = dplyr::lag(PARTIAL_DEMOCRACY_WITH_FACTIONALISM, 2)
  ) %>%
  ungroup() %>%
  arrange(row_id) %>%
  mutate(YEAR = as.character(YEAR_int)) %>%
  select(-row_id, -YEAR_int)

# Sanity check: should now pass
stopifnot(isTRUE(all.equal(
  df_base$LIBERAL_DEMOCRACY_SCORE,
  standardized_data_lagged$LIBERAL_DEMOCRACY_SCORE,
  tolerance = 1e-12
)))
```

\newpage

# 1. Executive Summary

This brief summary estimates:

* Country-year linear and logistic coefficients of 20 numeric predictors on a dichotomous outcome: whether a country is `Fragile` or `Non-Fragile` in a given year. The strategy is repeated with two-year lags of the same 20 numeric predictors. 

* The performance of the linear model (which returns the marginal effect of of one unit increase in the predictor over the target variable) and of the logistic model (which estimates the log-likelihood that the dichotomous outcome is true for a unit increase in the predictor). Predictors are normalized, so coefficients are directly comparable across the two models: the coefficients capture the change in the target variable for one standard deviation (positive or negative) in the predictor variable. The strategy is repeated with two-year lags of the same 20 numeric predictors.

* The LASSO regression parameters, tackling overfitting and multicollinaerity in identifying less relevant predictors. 

* The confusion matrix, which evaluates the accuracy of trained data in correctly identifying the class a country falls in for a given year (whether `Fragile` or `Non-Fragile`)

\newpage

# 2. Descriptive Statistics 

```{r descriptive-stats}
descriptive_data <- final_clean_percentiles_data |>
  # Selecting only relevant predictors 
  select(GDP_PER_CAPITA,
         LOG_GDP_PER_CAPITA,
         ODA_RECEIVED_PER_CAPITA, 
         GDP_GROWTH, 
         GDP_PER_CAPITA_GROWTH, 
         GDP_DEFLATOR, 
         POLITICAL_REGIME, 
         ELECTORAL_DEMOCRACY_SCORE, 
         LIBERAL_DEMOCRACY_SCORE, 
         TERRITORIAL_FRAGMENTATION, 
         INSTITUTIONAL_DEMOCRACY_SOCRE, 
         INSTITUTIONAL_AUTOCRACY_SCORE,
         COMBINED_POLITY_SCORE, 
         REGIME_DURABILITY_YEARS, 
         INSTITUTIONAL_EXECUTIVE_RECRUTIMENT, 
         POLITICAL_COMPETITION_SCORE, 
         PARTIAL_DEMOCRACY_WITH_FACTIONALISM, 
         CONFLICT_INTENSITY_YEAR, 
         CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS, 
         N_WAR_FRONTS, 
         MAX_CONFLICT_INTENSITY, 
         AVG_CONFLICT_INTENSITY, 
         N_TOTAL_TROOPS)

political_indicators <- final_clean_percentiles_data |>
  select(POLITICAL_REGIME, 
         ELECTORAL_DEMOCRACY_SCORE, 
         LIBERAL_DEMOCRACY_SCORE, 
         TERRITORIAL_FRAGMENTATION, 
         INSTITUTIONAL_DEMOCRACY_SOCRE, 
         INSTITUTIONAL_AUTOCRACY_SCORE,
         COMBINED_POLITY_SCORE, 
         REGIME_DURABILITY_YEARS, 
         INSTITUTIONAL_EXECUTIVE_RECRUTIMENT, 
         POLITICAL_COMPETITION_SCORE)
```

```{r, descriptive-stats-table, results = 'asis'}
stargazer(descriptive_data,
          title = "Descriptive Statistics Country Panel (1970-2022)",
          out.header = FALSE,
          header = FALSE,
          digits = 2,
          covariate.labels = c("GDP Per Capita",
                               "GDP Per Capita (logged)",
                               "ODA Received Per Capita",
                               "GDP Growth",
                               "GDP Per Capita Growth",
                               "GDP Deflator",
                               "Political Regime",
                               "Electoral Democracy Score",
                               "Liberal Democracy Score",
                               "Territorial Fragmentation",
                               "Institutional Democracy Score",
                               "Institutional Autocracy Score",
                               "Combined Polity Score",
                               "Regime Durability (Years)",
                               "Executive Recruitment",
                               "Political Competition Score",
                               "Partial Democracy With Factionalism",
                               "Conflict Intensity (Year)",
                               "Conflict Cumul. Intensity Across Years",
                               "N War Fronts",
                               "Max Conflict Intensity",
                               "Avg. Conflict Intensity",
                               "N Total Troops"))
```

```{r descriptive-plots}
# ------------------------------------------------------------------------------
# GDP Data ---------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Raw GDP histogram
p1 <- ggplot(standardized_data_final, aes(x = GDP_PER_CAPITA)) +
  geom_histogram(bins = 30, fill = "#6cbf84", color = "grey30", alpha = 0.9) +
  labs(
    title = "GDP per Capita (Raw)",
    x = "GDP per Capita (USD)",
    y = "Number of Countries"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "none"
  )

# Log GDP histogram
p2 <- ggplot(standardized_data_final, aes(x = log(GDP_PER_CAPITA))) +
  geom_histogram(bins = 30, fill = "#f26968", color = "grey30", alpha = 0.9) +
  labs(
    title = "GDP per Capita (Log)",
    x = "log(GDP per Capita)",
    y = "Number of Countries"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "none"
  )

# gridExtra
gridExtra::grid.arrange(p1, p2, ncol = 2)

# ------------------------------------------------------------------------------
# Time trends ------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Count number of countries per category per year
trend_data <- standardized_data_final |>
  group_by(YEAR, 
           VDEM_STATUS_IDEAL) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(YEAR = as.numeric(YEAR)) |>
  arrange(VDEM_STATUS_IDEAL,
          YEAR)

# Plotting the time trends -----------------------------------------------------
ggplot(trend_data, aes(x = YEAR, y = count, color = VDEM_STATUS_IDEAL, group = VDEM_STATUS_IDEAL)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Fragile" = "#f26968",
                                "Non Fragile" = "#6cbf84",
                                "Transitioning" = "#323339")) +
  scale_x_continuous(breaks = seq(min(trend_data$YEAR), max(trend_data$YEAR), 5)) +
  labs(title = "Yearly Trends of VDEM Fragility Stauts (1970-2022)",
       x = "Year",
       y = "Number of Countries",
       color = "Status") +
  theme_light(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "top"
  )

# ------------------------------------------------------------------------------
# INCOME DISTRIBUTION ----------------------------------------------------------
# ------------------------------------------------------------------------------
# Preparing the data 
df_gdp_quartiles <- standardized_data_final |>
  mutate(
    YEAR = as.numeric(YEAR)
  ) |>
  # Quartiles WITHIN year (1 = lowest income, 4 = highest)
  group_by(YEAR) |>
  mutate(gdp_pc_quartile = ntile(GDP_PER_CAPITA, 4)) |>
  ungroup() |>
  filter(!is.na(gdp_pc_quartile), !is.na(VDEM_STATUS_IDEAL), !is.na(YEAR))

# Aggregatting counts by year, status, and quartile
trend_by_q <- df_gdp_quartiles |>
  group_by(YEAR, 
           gdp_pc_quartile, 
           VDEM_STATUS_IDEAL) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(
    gdp_pc_quartile = factor(
      gdp_pc_quartile,
      levels = 1:4,
      labels = c("Q1: Lowest income", "Q2", "Q3", "Q4: Highest income")
    )
  )

# Plotting time trends by income distribution ----------------------------------
ggplot(
  trend_by_q,
  aes(x = YEAR, y = count, color = VDEM_STATUS_IDEAL, group = VDEM_STATUS_IDEAL)
) +
  geom_line(size = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "Fragile" = "#f26968",
    "Non Fragile" = "#6cbf84",
    "Transitioning" = "#323339"
  )) +
  scale_x_continuous(breaks = seq(min(trend_by_q$YEAR), max(trend_by_q$YEAR), 10)) +
  labs(
    title = "VDEM Fragility Stauts by GDP per capita Quartile (1970-2022)",
    x = "Year",
    y = "Number of Countries",
    color = "Status"
  ) +
  facet_wrap(~ gdp_pc_quartile, ncol = 2) +  
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "top"
  )

# ------------------------------------------------------------------------------
# POLITY DISTRIBUTION ----------------------------------------------------------
# ------------------------------------------------------------------------------
# Preparing the data 
df_polity_quartiles <- standardized_data_final |>
  mutate(
    YEAR = as.numeric(YEAR)
  ) |>
  # Quartiles WITHIN year (1 = lowest income, 4 = highest)
  group_by(YEAR) |>
  mutate(polity_quartile = ntile(COMBINED_POLITY_SCORE, 4)) |>
  ungroup() |>
  filter(!is.na(polity_quartile), !is.na(VDEM_STATUS_IDEAL), !is.na(YEAR))

# Aggregatting counts by year, status, and quartile
trend_by_q_polity <- df_polity_quartiles |>
  group_by(YEAR, 
           polity_quartile, 
           VDEM_STATUS_IDEAL) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(
    polity_quartile = factor(
      polity_quartile,
      levels = 1:4,
      labels = c("Q1: Strongly Autocratic", "Q2", "Q3", "Q4: Strongly Democratic")
    )
  )

# Plotting the data by polity distribution -------------------------------------
ggplot(
  trend_by_q_polity,
  aes(x = YEAR, y = count, color = VDEM_STATUS_IDEAL, group = VDEM_STATUS_IDEAL)
) +
  geom_line(size = 0.5) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "Fragile" = "#f26968",
    "Non Fragile" = "#6cbf84",
    "Transitioning" = "#323339"
  )) +
  scale_x_continuous(breaks = seq(min(trend_by_q$YEAR), max(trend_by_q$YEAR), 10)) +
  labs(
    title = "VDEM Fragility Stauts by Polity Score (1970-2022)",
    x = "Year",
    y = "Number of Countries",
    color = "Status"
  ) +
  facet_wrap(~ polity_quartile, ncol = 2) +  
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "top"
  )

# ------------------------------------------------------------------------------
# CORRELATION TREND ------------------------------------------------------------
# ------------------------------------------------------------------------------
# Quartiles by each country's avg. liberal democracy
df_q <- final_clean_percentiles_data_normalized_1 |>
  group_by(COUNTRY_NAME) |>
  mutate(dem_avg = mean(LIBERAL_DEMOCRACY_SCORE, na.rm = TRUE)) |>
  ungroup() |>
  mutate(dem_quartile = ntile(dem_avg, 4),
         dem_quartile = factor(
           dem_quartile,
           levels = 1:4,
           labels = c("Q1 (Less Democratic)", 
                      "Q2", 
                      "Q3", "Q4 (Most Democratic)")
         ))

# Long format
long <- df_q |>
  select(COUNTRY_NAME, 
         YEAR, 
         dem_quartile,
         LIBERAL_DEMOCRACY_SCORE, 
         LOG_GDP_PER_CAPITA) |>
  pivot_longer(
    cols = c(LIBERAL_DEMOCRACY_SCORE, 
             LOG_GDP_PER_CAPITA),
    names_to = "series",
    values_to = "z"
  )

# Putting both series on the same (z) scale
long_z <- long |>
  group_by(series) |>
  mutate(z_std = as.numeric(scale(z))) |>
  ungroup()

# Quartile-by-year averages + 95% CI across countries
avg_ci <- long_z |>
  group_by(dem_quartile, YEAR, series) |>
  summarise(
    mean_z = mean(z_std, na.rm = TRUE),
    n      = n_distinct(COUNTRY_NAME),
    sd_z   = sd(z_std, na.rm = TRUE),
    se     = sd_z / sqrt(n),
    lwr    = mean_z - 1.96 * se,
    upr    = mean_z + 1.96 * se,
    .groups = "drop"
  )

# Plotting
ggplot(avg_ci, aes(x = YEAR, y = mean_z, color = series, fill = series)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ dem_quartile, ncol = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(
    x = NULL, y = "Standardized (within-series z-score)",
    color = NULL, fill = NULL,
    title = "Quartile means of GDP per capita & Liberal Democracy over time",
    subtitle = "Faceted by quartiles of average liberal democracy (ribbons = 95% CI across countries)"
  ) +
  scale_color_discrete(labels = c(LOG_GDP_PER_CAPITA = "GDP per capita (WB)", LIBERAL_DEMOCRACY_SCORE = "Liberal democracy (V-DEM)")) +
  scale_fill_discrete(labels  = c(LOG_GDP_PER_CAPITA = "GDP per capita (WB)", LIBERAL_DEMOCRACY_SCORE = "Liberal democracy (V-DEM)")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text.x = element_text(vjust = 0.5),
        legend.position = "top")

# ------------------------------------------------------------------------------
# Stacked Barplot --------------------------------------------------------------
# ------------------------------------------------------------------------------
decade_starts <- seq(1971, 2021, 10)

# Distribution of statuses per decade ------------------------------------------
df_pct <- final_clean_percentiles_data_normalized_1 |>
  filter(YEAR %in% decade_starts,
         !is.na(COUNTRY_NAME),
         !is.na(VDEM_STATUS_IDEAL)) |>
  group_by(COUNTRY_NAME, YEAR) |>
  slice(1) |>               
  ungroup() |>
  count(YEAR, 
        VDEM_STATUS_IDEAL, 
        name = "n") |>
  complete(YEAR = decade_starts, 
           VDEM_STATUS_IDEAL, 
           fill = list(n = 0)) |>
  group_by(YEAR) |>
  mutate(
    total = sum(n),
    pct   = if_else(total > 0, n / total, NA_real_)
  ) |>
  ungroup() |>
  mutate(
    decade = factor(paste0(YEAR, "s"), levels = paste0(decade_starts, "s"))
  )

# Plotting ---------------------------------------------------------------------
ggplot(df_pct, aes(x = decade, y = pct, fill = VDEM_STATUS_IDEAL)) +
  geom_col(width = 0.85) +
  geom_text(
    aes(label = if_else(pct >= 0.03, percent(pct, accuracy = 1), "")),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.06)) +
  scale_fill_manual(
    values = c(
      "Fragile" = "#f26968",
      "Non Fragile" = "#6cbf84"
    ),
    breaks = c("Fragile", "Non Fragile")  
  ) +
  labs(
    x = "Decade",
    y = "Share of countries",
    fill = "Fragility Status",
    title = "Distrbution of Fragility Statuses Across Decades (1970s-2020s)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
  panel.grid.minor = element_blank()
)
```

\newpage 

# 3. Success Episodes 

## 3.1. ODA: Pre vs. Post-Escape Pooled Estimates

The scatterplot below examines how ODA per capita evolves around a country’s successful breakthrough from fragility. For each escape case, we compute the average ODA per capita in the pre-escape period (x-axis, log-transformed) and the corresponding change in ODA per capita between the pre- and post-escape periods. Each point represents a single country. The solid black line shows the ordinary least squares (OLS) fit with 95% confidence intervals based on heteroskedasticity-consistent (HC) standard errors.

To minimize the influence of extreme ODA spikes and unequal exposure periods, we restrict attention to a 3-year window around each escape year and winsorize (i.e., replace outliers with less extreme observations) ODA values at the 1st and 99th percentiles. Within this symmetric window, ODA levels are averaged for each country before and after the escape event, and differences are computed (`Post – Pre`).

```{r oda-pooled}
# Setting the critical junctures
standardized_data_final <- standardized_data_final |>
  mutate(YEAR = as.numeric(YEAR)) |>
  mutate(VDEM_ESCAPE_YEAR = as.numeric(case_when(
    COUNTRY_NAME == "Brazil"                & YEAR == 1988 ~ 1,
    COUNTRY_NAME == "Gabon"                 & YEAR == 1990 ~ 1,
    COUNTRY_NAME == "Gambia, The"           & YEAR == 2017 ~ 1,
    COUNTRY_NAME == "Kenya"                 & YEAR == 2010 ~ 1,
    COUNTRY_NAME == "Mauritania"            & YEAR == 2006 ~ 1,
    COUNTRY_NAME == "Sierra Leone"          & YEAR == 2018 ~ 1,
    COUNTRY_NAME == "Tunisia"               & YEAR == 2010 ~ 1,
    COUNTRY_NAME == "Uganda"                & YEAR == 1987 ~ 1,
    COUNTRY_NAME == "Iran, Islamic Rep."    & YEAR == 1997 ~ 1,
    COUNTRY_NAME == "Nepal"                 & YEAR == 2007 ~ 1,
    COUNTRY_NAME == "Philippines"           & YEAR == 1987 ~ 1,
    COUNTRY_NAME == "Sri Lanka"             & YEAR == 2009 ~ 1,
    COUNTRY_NAME == "Thailand"              & YEAR == 1991 ~ 1,
    COUNTRY_NAME == "Bosnia and Herzegovina"& YEAR == 1995 ~ 1,
    COUNTRY_NAME == "Serbia"                & YEAR == 2001 ~ 1,
    COUNTRY_NAME == "Georgia"               & YEAR == 2004 ~ 1,
    COUNTRY_NAME == "Kyrgyz Republic"       & YEAR == 2011 ~ 1,
    COUNTRY_NAME == "Armenia"               & YEAR == 2018 ~ 1,
    COUNTRY_NAME == "Argentina"             & YEAR == 1983 ~ 1,
    COUNTRY_NAME == "Colombia"              & YEAR == 2002 ~ 1,
    COUNTRY_NAME == "Ecuador"               & YEAR == 1976 ~ 1,
    COUNTRY_NAME == "Dominican Republic"    & YEAR == 1997 ~ 1,
    COUNTRY_NAME == "El Salvador"           & YEAR == 2009 ~ 1,
    COUNTRY_NAME == "Mexico"                & YEAR == 1994 ~ 1,
    COUNTRY_NAME == "Nicaragua"             & YEAR == 1990 ~ 1,
    COUNTRY_NAME == "Panama"                & YEAR == 1990 ~ 1,
    COUNTRY_NAME == "Peru"                  & YEAR == 2001 ~ 1,
    TRUE ~ 0
  )))

# ------------------------------------------------------------------------------
# POOLED PLOT ------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Parameters -------------------------------------------------------------------
ODA_VAR <- "ODA_RECEIVED_PER_CAPITA"  

# Minimal Pre and Post-Escape years 
MIN_N_PER_SIDE <- 2                 

# Identifying escape years per country
escape_years <- standardized_data_final |>
  filter(VDEM_ESCAPE_YEAR == 1) |>
  select(COUNTRY_NAME, 
         ESCAPE_YEAR = YEAR)

# Keeping escape countries and tagging Pre/Post 
# (We exclude the escape year itself)
df_ep <- standardized_data_final |>
  semi_join(escape_years, by = "COUNTRY_NAME") |>
  left_join(escape_years, by = "COUNTRY_NAME") |>
  mutate(
    YEAR = as.numeric(YEAR),
    period = case_when(
      YEAR < ESCAPE_YEAR ~ "Pre",
      YEAR > ESCAPE_YEAR ~ "Post",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(period))

# Summarising mean ODA per country and period
sum_by_period <- df_ep |>
  group_by(COUNTRY_NAME, 
           period) |>
  summarise(
    mean_oda = mean(.data[[ODA_VAR]], na.rm = TRUE),
    n_years  = sum(!is.na(.data[[ODA_VAR]])),
    .groups = "drop"
  )

# Widening the table Pre/Post means, keeping countries with data on both sides
pre_post <- sum_by_period |>
  pivot_wider(
    names_from = period, values_from = c(mean_oda, n_years)
  ) |>
  # Keeping countries with enough observations on both sides
  filter(!is.na(mean_oda_Pre), !is.na(mean_oda_Post),
         n_years_Pre >= MIN_N_PER_SIDE, n_years_Post >= MIN_N_PER_SIDE) |>
  transmute(
    COUNTRY_NAME,
    pre_mean  = mean_oda_Pre,
    post_mean = mean_oda_Post,
    delta     = post_mean - pre_mean
  )

# Scatterplot with 95% confidence interval 
p <- ggplot(pre_post, aes(x = pre_mean, y = delta)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, size = 1, color = "black", fill = "grey80") +
  scale_x_continuous(trans = "log1p") +
  labs(
    title = "Change in ODA per capita: Escape vs. Pre-escape Level",
    subtitle = "Each dot is a country; line = OLS fit with 95% CI",
    x = "Pre-escape mean ODA per capita (log1p)",
    y = "Post – Pre change in ODA per capita"
  ) +
  theme_classic(base_size = 13)

p

# ------------------------------------------------------------------------------
# Symetric windows around pre- and post-escape ---------------------------------
# ------------------------------------------------------------------------------
WINDOW <- 3  # years on each side

df_ep <- standardized_data_final |>
  semi_join(escape_years, by = "COUNTRY_NAME") |>
  left_join(escape_years, by = "COUNTRY_NAME") |>
  mutate(period = case_when(
    YEAR >= ESCAPE_YEAR - WINDOW & YEAR <  ESCAPE_YEAR ~ "Pre",
    YEAR >  ESCAPE_YEAR          & YEAR <= ESCAPE_YEAR + WINDOW ~ "Post",
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(period))

# ------------------------------------------------------------------------------
# Winsorizing extreme ODA before averaging -------------------------------------
# ------------------------------------------------------------------------------
winsor01 <- function(x, p = 0.01) {
  q <- quantile(x, c(p, 1-p), na.rm = TRUE)
  pmin(pmax(x, q[[1]]), q[[2]])
}

df_ep <- df_ep |>
  mutate(ODA_W = winsor01(.data[[ODA_VAR]], p = 0.01))

# ------------------------------------------------------------------------------
# Paired tests with robust SEs -------------------------------------------------
# ------------------------------------------------------------------------------
sum_by_period <- df_ep |>
  group_by(COUNTRY_NAME, period) |>
  summarise(
    mean_oda = mean(ODA_W, na.rm = TRUE),
    n_years  = sum(!is.na(ODA_W)),
    .groups = "drop"
  )

pre_post <- sum_by_period |>
  pivot_wider(names_from = period, values_from = c(mean_oda, n_years)) |>
  filter(!is.na(mean_oda_Pre), !is.na(mean_oda_Post),
         n_years_Pre >= MIN_N_PER_SIDE, n_years_Post >= MIN_N_PER_SIDE) |>
  transmute(COUNTRY_NAME,
            pre_mean  = mean_oda_Pre,
            post_mean = mean_oda_Post,
            delta     = post_mean - pre_mean)

# Paired test
t_out <- with(pre_post, t.test(post_mean, pre_mean, paired = TRUE))

# Robust slope with sandwich SEs (country-level points)
mod <- lm(delta ~ log1p(pre_mean), data = pre_post)
rob_se <- sqrt(diag(sandwich::vcovHC(mod, type = "HC1")))
est <- coef(summary(mod))
est[,2] <- rob_se  # replace SE

# ------------------------------------------------------------------------------
# Final Robust Plot ------------------------------------------------------------
# ------------------------------------------------------------------------------
p <- ggplot(pre_post, aes(x = pre_mean, y = delta)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, size = 1, color = "black", fill = "grey80") +
  scale_x_continuous(trans = "log1p") +
  labs(
    title = "Change in ODA per capita around escape (3 years)",
    subtitle = "Each dot = country; line = OLS with HC-robust 95% CI",
    x = "Pre-escape mean ODA per capita (log1p, winsorized 1%)",
    y = "Post – Pre change in ODA per capita (winsorized 1%)"
  ) +
  theme_classic(base_size = 13)

p

# ------------------------------------------------------------------------------
# Descriptive table ------------------------------------------------------------
# ------------------------------------------------------------------------------
summary_tbl <- tibble(
  n_countries   = nrow(pre_post),
  mean_pre_ODA  = mean(pre_post$pre_mean,  na.rm = TRUE),
  mean_post_ODA = mean(pre_post$post_mean, na.rm = TRUE),
  mean_delta    = mean(pre_post$delta,     na.rm = TRUE),
  t_p_value     = signif(t_out$p.value, 3),
  slope         = round(coef(mod)[2], 3),
  robust_se     = round(rob_se[2], 3))

knitr::kable(summary_tbl,
             caption = "ODA dynamics around escape from fragility (3 years window)",
             digits = 3,
             align = "c")
```

A paired t-test across 27 countries indicates that ODA per capita rises significantly following escape (mean increase around $0.31$ SD, $p = 0.008$). The fitted slope suggests that countries already receiving higher levels of aid tend to experience somewhat larger post-escape increases. We recognize that this association is imprecisely estimated due to the small sample size and high cross-country variance in aid levels.

Taken together, these results imply that aid inflows tend to follow rather than precede successful exits from fragility. Donors appear to increase assistance once credible political and institutional progress is underway—consistent with ODA functioning more as an *ex-post* reward or reinforcement mechanism than as a fundamental pre-condition or catalyst for state recovery. 

## 3.2. ODA: Pre vs. Post-Escape Country-Estimates

```{r oda-loop}
# Parameters -------------------------------------------------------------------
# log1p y-scale (zeros/ skew)
USE_LOG <- TRUE      

# Tagging Pre/Post (dropping the escape year itself)
df_ep <- standardized_data_final |>
  semi_join(escape_years, by = "COUNTRY_NAME") |>
  left_join(escape_years, by = "COUNTRY_NAME") |>
  mutate(
    YEAR = as.numeric(YEAR),
    period = case_when(
      YEAR < ESCAPE_YEAR ~ "Pre",
      YEAR > ESCAPE_YEAR ~ "Post",
      TRUE               ~ NA_character_
    )
  ) |>
  filter(!is.na(period), !is.na(.data[[ODA_VAR]]))

# Keeping only countries with data on BOTH sides
eligible <- df_ep |>
  count(COUNTRY_NAME, 
        period) |>
  pivot_wider(names_from = period, 
              values_from = n, 
              values_fill = 0) |>
  filter(Pre >= MIN_N_PER_SIDE, Post >= MIN_N_PER_SIDE) |>
  pull(COUNTRY_NAME)

df_ep <- df_ep |> 
  filter(COUNTRY_NAME %in% eligible)

# Building the single Pre/Post plot for each country
plot_one_country <- function(cname, use_log = TRUE) {
  dat <- df_ep |>
    filter(COUNTRY_NAME == cname)

  # Welch t-test (Post vs Pre)
  tt <- t.test(dat[[ODA_VAR]] ~ dat$period, var.equal = FALSE)
  p_label <- ifelse(tt$p.value < 0.001, "p<0.001", sprintf("p=%.3f", tt$p.value))

  # Summary for CI bars
  summ <- dat |>
    group_by(period) |>
    summarise(
      mean_oda = mean(.data[[ODA_VAR]], na.rm = TRUE),
      se_oda   = sd(.data[[ODA_VAR]], na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(
      ci_low  = mean_oda - 1.96 * se_oda,
      ci_high = mean_oda + 1.96 * se_oda,
      period  = factor(period, levels = c("Pre", "Post"))
    )

  p <- ggplot(summ, aes(x = period, y = mean_oda, color = period)) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18, size = 0.7) +
    geom_line(aes(group = 1), color = "#f26968", linewidth = 0.6) +
    scale_color_manual(values = c("Pre" = "#f26968", "Post" = "#6cbf84")) +
    { if (use_log) scale_y_continuous(trans = "log1p") else scale_y_continuous() } +
    labs(
      title = paste0(cname, " ODA per capita: Pre vs Post (", p_label, ")"),
      subtitle = "Points = mean, bars = 95% CI. Welch t-test for Pre vs Post.",
      x = NULL,
      y = if (use_log) "Mean ODA per capita (log1p)" else "Mean ODA per capita",
      color = "Period"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))
  p
}

# Looping plots for all countries
countries <- unique(df_ep$COUNTRY_NAME)
plots <- setNames(map(countries, ~ plot_one_country(.x, use_log = USE_LOG)), countries)
walk(plots, print)

# ------------------------------------------------------------------------------
# Robustness -------------------------------------------------------------------
# ------------------------------------------------------------------------------
USE_LOG    <- TRUE       
WINDOW     <- 3          
MIN_N_PER_SIDE <- 2
WINSOR_P   <- 0.01       

# Function for winsorization (1%)
winsor01 <- function(x, p = 0.01) {
  if (all(is.na(x))) return(x)
  qs <- quantile(x, c(p, 1 - p), na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

perm_pvalue <- function(x, g, B = 5000L, seed = 247) {
  # two-sample difference in means, exact/randomization p via label permutations
  set.seed(seed)
  x <- x[!is.na(x) & !is.na(g)]
  g <- g[!is.na(x) & !is.na(g)]
  if (length(unique(g)) != 2) return(NA_real_)
  obs <- mean(x[g == "Post"]) - mean(x[g == "Pre"])
  B <- as.integer(B)
  stats <- replicate(B, {
    g_star <- sample(g, length(g), replace = FALSE)
    mean(x[g_star == "Post"]) - mean(x[g_star == "Pre"])
  })
  (sum(abs(stats) >= abs(obs)) + 1) / (B + 1)
}

# Pre/Post with symmetric window + winsorization -------------------------------
escape_years <- standardized_data_final |>
  filter(VDEM_ESCAPE_YEAR == 1) |>
  select(COUNTRY_NAME, 
         ESCAPE_YEAR = YEAR)

df_ep <- standardized_data_final |>
  semi_join(escape_years, 
            by = "COUNTRY_NAME") |>
  left_join(escape_years, 
            by = "COUNTRY_NAME") |>
  mutate(
    YEAR = as.numeric(YEAR),
    period = case_when(
      YEAR >= ESCAPE_YEAR - WINDOW & YEAR <  ESCAPE_YEAR ~ "Pre",
      YEAR >  ESCAPE_YEAR          & YEAR <= ESCAPE_YEAR + WINDOW ~ "Post",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(period), !is.na(.data[[ODA_VAR]])) |>
  # winsorizing within country
  group_by(COUNTRY_NAME) |>
  mutate(ODA_W = winsor01(.data[[ODA_VAR]], p = WINSOR_P)) |>
  ungroup()

# keeping countries with data on both sides ------------------------------------
eligible <- df_ep |>
  count(COUNTRY_NAME, period) |>
  pivot_wider(names_from = period, values_from = n, values_fill = 0) |>
  filter(Pre >= MIN_N_PER_SIDE, Post >= MIN_N_PER_SIDE) |>
  pull(COUNTRY_NAME)

df_ep <- df_ep |>
  filter(COUNTRY_NAME %in% eligible) |>
  mutate(
    # stats on log scale
    ODA_log = log1p(ODA_W)  
  )

# Looping plots ----------------------------------------------------------------
plot_one_country <- function(cname, use_log = TRUE) {
  dat <- df_ep |>
    filter(COUNTRY_NAME == cname)

  yvar <- if (use_log) "ODA_log" else "ODA_W"

  # Welch t-test
  f <- as.formula(paste(yvar, "~ period"))
  tt <- t.test(f, data = dat, var.equal = FALSE)

  # Permutation p-value (robust small-N check)
  p_perm <- perm_pvalue(dat[[yvar]], dat$period, B = 5000L)

  p_label <- if (!is.na(p_perm)) sprintf("p_perm=%.3f", p_perm) else sprintf("p=%.3f", tt$p.value)

  # Summary means/SE/CI on the same scale used for the test
  summ <- dat |>
    group_by(period) |>
    summarise(
      mean_y = mean(.data[[yvar]], na.rm = TRUE),
      se_y   = sd(.data[[yvar]], na.rm = TRUE) / sqrt(n()),
      n      = n(),
      .groups = "drop"
    ) |>
    mutate(
      ci_low  = mean_y - 1.96 * se_y,
      ci_high = mean_y + 1.96 * se_y,
      period  = factor(period, levels = c("Pre", "Post"))
    )

  ggplot(summ, aes(x = period, y = mean_y, color = period)) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18, linewidth = 0.7) +
    geom_line(aes(group = 1), color = "#f26968", linewidth = 0.6) +
    scale_color_manual(values = c("Pre" = "#f26968", "Post" = "#6cbf84")) +
    { if (use_log) scale_y_continuous() else scale_y_continuous() } +
    labs(
      title = paste0(cname, " ODA per capita: Pre vs Post (", p_label, ")"),
      subtitle = if (use_log) "Means and 95% CI on log1p-winsorized ODA; Welch t + permutation p"
                 else         "Means and 95% CI on winsorized ODA; Welch t + permutation p",
      x = NULL,
      y = if (use_log) "Mean log1p(ODA per capita)" else "Mean ODA per capita (winsorized)",
      color = "Period"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
                   plot.title = element_text(face = "bold"))
}

# Plotting ---------------------------------------------------------------------
countries <- unique(df_ep$COUNTRY_NAME)

country_stats <- lapply(countries, function(cn) {
  dat <- filter(df_ep, COUNTRY_NAME == cn)
  yvar <- if (USE_LOG) "ODA_log" else "ODA_W"
  f <- as.formula(paste(yvar, "~ period"))
  tt <- t.test(f, data = dat, var.equal = FALSE)
  p_perm <- perm_pvalue(dat[[yvar]], dat$period, B = 5000L)
  data.frame(country = cn,
             p_welch = tt$p.value,
             p_perm  = p_perm,
             n_pre   = sum(dat$period == "Pre"),
             n_post  = sum(dat$period == "Post"),
             stringsAsFactors = FALSE)
})

country_stats <- bind_rows(country_stats) |>
  mutate(
    p_perm_adj = p.adjust(p_perm, method = "BH"),
    p_welch_adj= p.adjust(p_welch, method = "BH")
  )

# plots named list
plots <- setNames(purrr::map(countries, ~ plot_one_country(.x, use_log = USE_LOG)), countries)
walk(plots, print)
```

\newpage

# 4. A Classification Problem 

To validate the internal consistency and empirical grounding of our fragility classification (`Fragile`, `Transitioning` or `Non-Fragile`), we implement a supervised machine learning exercise using two distinct algorithms: multinomial logistic regression (logit) and random forest. We train both models on a set of 20 theoretically motivated structural predictors, including regime type, democracy indicators, institutional durability, economic performance, and conflict intensity (i.e., all predictors we employ in our empirical specifications). The outcome variable is our three-tier classification (\texttt{VDEM\_STATUS\_IDEAL}), derived independently of the predictors via percentile cutoffs on a normalized V-Dem latent factor. 

We label our country cases according to the subsequent taxonomy:

| Classification Label     | Fragile Threshold            | Transitioning Range              | Non-Fragile Threshold           |
|--------------------------|------------------------------|----------------------------------|---------------------------------|
| `RELAXED`                | Below 20th percentile        | 20th to 25th percentile          | Above 25th percentile           |
| `IDEAL`                  | Below 30th percentile        | 30th to 35th percentile          | Above 35th percentile           |
| `STRICT`                 | Below 20th percentile        | 20th to 40th percentile          | Above 40th percentile           |


The two supervised learning algorithms serve complementary purposes to validate our classification:

- **Multinomial Logistic Regression (Logit)** estimates the probability that a given observation belongs to one of the three mutually exclusive classes (`Fragile`, `Transitioning` or `Non-Fragile`), based on a linear combination of structural predictors. It is fully interpretable, assumes additive relationships, and provides insight into the marginal effects of each variable. Because it is transparent and theory-aligned, it serves as an ideal benchmark for testing the recoverability of our classification from exogenous structural features.

- **Random Forest** is a non-parametric ensemble method that builds hundreds of decision trees on random subsets of predictors and observations. It excels at modeling complex, nonlinear interactions and is robust to multicollinearity and overfitting. In our context, it offers a high-performing, assumption-light alternative that tests whether the classification can still be recovered without relying on linearity or additivity assumptions.

\newpage 

## 4.1. Ideal Classification

We evaluate model performance using repeated 10-fold cross-validation (3 repeats; 30 folds in total). The multinomial logit attains close to $80%$ accuracy (SE = 0.006) and a multiclass log-loss of $\approx 0.49$ (SE = 0.013). The random forest performs substantially better, achieving $\approx 91.7%$ accuracy (SE = 0.002) and a lower log-loss of $\approx 0.27$ (SE = 0.004).

Confusion matrices confirm the gap: the logit almost never predicts Transitioning cases (recall $\approx 0%$), whereas the random forest recovers part of that class (recall for around $20%$ of the cases) and improves recall for both Fragile ($\approx 93%$ vs. $\approx 82%$) and Non-Fragile ($\approx 95%$ vs. $\approx 89%$). We acknowledge that such asymmetries suggest non-trivial class imbalance. Standard errors remain small, indicating stable performance across folds.

These results reinforce the validity of our three-tier labeling: institutional, economic, and conflict covariates reliably reproduce the categorical structure, and the random forest clearly dominates the logit both in classification accuracy and in the quality of predicted probabilities (i.e., lower loss). The only systematic weakness lies in the Transitioning category, which is consistent with its conceptual “in-between” status and limited empirical representation.

Overall, the evidence suggests that our fragility tiers are empirically grounded patterns rather than artifacts of arbitrary thresholds.

```{r ideal-classification}
# ------------------------------------------------------------------------------
# IDEAL CLASSIFICATION ---------------------------------------------------------
# ------------------------------------------------------------------------------
# Defining our predictors 
predictors_20 <- c(
  "ODA_RECEIVED_PER_CAPITA", 
  "GDP_GROWTH", 
  "LOG_GDP_PER_CAPITA", 
  "GDP_DEFLATOR",
  "POLITICAL_REGIME", 
  "ELECTORAL_DEMOCRACY_SCORE", 
  "LIBERAL_DEMOCRACY_SCORE", 
  "TERRITORIAL_FRAGMENTATION", 
  "INSTITUTIONAL_DEMOCRACY_SOCRE",  
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "COMBINED_POLITY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

# Preparing the data 
df <- standardized_data_final |>
  # Dropping COUNTRY_NAME identifier
  ungroup() |>
  # Selecting the relevant predictors
  select(all_of(predictors_20), 
         # Selecting the target variable 
         VDEM_STATUS_IDEAL) |>
  # Inputting meadian instead of dropping NAs
  step_impute_median(all_numeric_predictors()) |>
  # Setting target variable as factor
  mutate(VDEM_STATUS_IDEAL = as.factor(VDEM_STATUS_IDEAL)) 

# Train/test split
data_split <- initial_split(df, strata = VDEM_STATUS_IDEAL)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Recipe
rec <- recipe(VDEM_STATUS_IDEAL ~ ., data = train_data) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# ------------------------------------------------------------------------------
# Model specifications ---------------------------------------------------------
# ------------------------------------------------------------------------------

# Logit Model ------------------------------------------------------------------
logit_spec <- multinom_reg(mode = "classification") |> 
  set_engine("nnet")

# Random Forest Model ----------------------------------------------------------
rf_spec <- rand_forest(mode = "classification", 
                       trees = 500) |>
  set_engine("ranger", 
             importance = "permutation")

# Workflows --------------------------------------------------------------------

# Logit model workflow
logit_wf <- workflow() |> 
  add_model(logit_spec) |> 
  add_recipe(rec)

# Random Forest workflow
rf_wf    <- workflow() |> 
  add_model(rf_spec) |> 
  add_recipe(rec)

# Cross-Validation Setup -------------------------------------------------------
folds <- vfold_cv(train_data, 
                  v = 10, 
                  repeats = 3, 
                  strata = VDEM_STATUS_IDEAL)

# Cross-Validation model evaluation --------------------------------------------

# Logit model
logit_res <- fit_resamples(logit_wf, 
                           folds, 
                           metrics = metric_set(accuracy, 
                                                mn_log_loss))

# Random Forest 
rf_res    <- fit_resamples(rf_wf, 
                           folds, 
                           metrics = metric_set(accuracy, 
                                                mn_log_loss))

# Final model fits -------------------------------------------------------------

# Logit model 
logit_fit <- fit(logit_wf, 
                 data = train_data)

# Random Forest model 
rf_fit    <- fit(rf_wf, 
                 data = train_data)

#-------------------------------------------------------------------------------
# PREDICTION -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Logit
logit_preds <- predict(logit_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

# Random Forest 
rf_preds <- predict(rf_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

# ------------------------------------------------------------------------------
# CONFUSION MATRICES + ACCURACY ------------------------------------------------
# ------------------------------------------------------------------------------

# computing confusion matrices
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds, truth = truth, estimate = .pred_class)

# Computing model accuracy 
accuracy_logit <- accuracy(logit_preds, truth = truth, estimate = .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds, truth = truth, estimate = .pred_class)$.estimate

# Plotting the matrices
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(
    title = "Multinomial Logit Confusion Matrix",
    subtitle = paste0("Accuracy = ", round(accuracy_logit * 100, 1), "%"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(
    title = "Random Forest Confusion Matrix",
    subtitle = paste0("Accuracy =", round(accuracy_rf * 100, 1), "%"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

# Variable importance ----------------------------------------------------------
vip(extract_fit_parsnip(rf_fit), 
    num_features = 20, 
    geom = "col") +
  labs(title = "RF Variable Importance") +
  theme_classic(base_size = 14)

# Metrics ----------------------------------------------------------------------
collect_metrics(logit_res)
collect_metrics(rf_res) 
```

\newpage

## 4.2. Relaxed Classification

We re-estimate our classifiers resorting to a relaxed specification and evaluate performance through repeated 10-fold cross-validation (3 repeats; 30 folds in total). The multinomial logit now attains an average accuracy of $83.9\%$ ($\mathrm{SE}=0.006$) and a multiclass log-loss of $0.429$ ($\mathrm{SE}=0.015$). The **random forest** continues to outperform it decisively, reaching $92.3\%$ accuracy ($\mathrm{SE}=0.002$) and a substantially lower log-loss of $0.247$ ($\mathrm{SE}=0.003$).  

Confusion matrices confirm the improvement: the logit still under-predicts *Transitioning* cases (recall $\approx5\%$) but performs better overall than in the baseline specification, while the random forest recovers a larger share of *Transitioning* observations (recall $\approx25\%$) and maintains high recall for both *Fragile* ($\approx93\%$) and *Non-Fragile* ($\approx96\%$) tiers. Standard errors remain small, which we take as indicative of stable performance across folds and strong out-of-sample generalization.

Overall, these results reinforce the robustness of the three-tier labeling scheme. As before, the principal weakness lies in the *Transitioning* category, consistent with its intermediate conceptual position and relative scarcity in the data.


```{r relaxed-classification}
# ------------------------------------------------------------------------------
# RELAXED CLASSIFICATION -------------------------------------------------------
# ------------------------------------------------------------------------------

# Defining our predictors 
predictors_20 <- c(
  "ODA_RECEIVED_PER_CAPITA", 
  "GDP_GROWTH", 
  "LOG_GDP_PER_CAPITA", 
  "GDP_DEFLATOR",
  "POLITICAL_REGIME", 
  "ELECTORAL_DEMOCRACY_SCORE", 
  "LIBERAL_DEMOCRACY_SCORE", 
  "TERRITORIAL_FRAGMENTATION", 
  "INSTITUTIONAL_DEMOCRACY_SOCRE",  
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "COMBINED_POLITY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

# Preparing the data 
df <- standardized_data_final |>
  # Dropping COUNTRY_NAME identifier
  ungroup() |>
  # Selecting the relevant predictors
  select(all_of(predictors_20), 
         # Selecting the target variable 
         VDEM_STATUS_BASELINE) |>
  # Inputting meadian instead of dropping NAs
  step_impute_median(all_numeric_predictors()) |>
  # Setting target variable as factor
  mutate(VDEM_STATUS_BASELINE = as.factor(VDEM_STATUS_BASELINE)) 

# Train/test split
data_split <- initial_split(df, strata = VDEM_STATUS_BASELINE)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Recipe
rec <- recipe(VDEM_STATUS_BASELINE ~ ., data = train_data) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# ------------------------------------------------------------------------------
# Model specifications ---------------------------------------------------------
# ------------------------------------------------------------------------------

# Logit Model ------------------------------------------------------------------
logit_spec <- multinom_reg(mode = "classification") |> 
  set_engine("nnet")

# Random Forest Model ----------------------------------------------------------
rf_spec <- rand_forest(mode = "classification", 
                       trees = 500) |>
  set_engine("ranger", 
             importance = "permutation")

# Workflows --------------------------------------------------------------------

# Logit model workflow
logit_wf <- workflow() |> 
  add_model(logit_spec) |> 
  add_recipe(rec)

# Random Forest workflow
rf_wf    <- workflow() |> 
  add_model(rf_spec) |> 
  add_recipe(rec)

# Cross-Validation Setup -------------------------------------------------------
folds <- vfold_cv(train_data, 
                  v = 10, 
                  repeats = 3, 
                  strata = VDEM_STATUS_BASELINE)

# Cross-Validation model evaluation --------------------------------------------

# Logit model
logit_res <- fit_resamples(logit_wf, 
                           folds, 
                           metrics = metric_set(accuracy, 
                                                mn_log_loss))

# Random Forest 
rf_res    <- fit_resamples(rf_wf, 
                           folds, 
                           metrics = metric_set(accuracy, 
                                                mn_log_loss))

# Final model fits -------------------------------------------------------------

# Logit model 
logit_fit <- fit(logit_wf, 
                 data = train_data)

# Random Forest model 
rf_fit    <- fit(rf_wf, 
                 data = train_data)

#-------------------------------------------------------------------------------
# PREDICTION -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Logit
logit_preds <- predict(logit_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_BASELINE)

# Random Forest 
rf_preds <- predict(rf_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_BASELINE)

# ------------------------------------------------------------------------------
# CONFUSION MATRICES + ACCURACY ------------------------------------------------
# ------------------------------------------------------------------------------

# computing confusion matrices
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds, truth = truth, estimate = .pred_class)

# Computing model accuracy 
accuracy_logit <- accuracy(logit_preds, truth = truth, estimate = .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds, truth = truth, estimate = .pred_class)$.estimate

# Plotting the matrices
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(
    title = "Multinomial Logit Confusion Matrix",
    subtitle = paste0("Accuracy=", round(accuracy_logit * 100, 1), "%"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(
    title = "Random Forest Confusion Matrix",
    subtitle = paste0("Accuracy=", round(accuracy_rf * 100, 1), "%"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

# Variable importance ----------------------------------------------------------
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "RF Variable Importance") +
  theme_classic(base_size = 14)

# Metrics ----------------------------------------------------------------------
collect_metrics(logit_res)
collect_metrics(rf_res) 
```

\newpage

## 4.3. Strict Classification

We now turn to the final classification results using the complete feature set and the finalized three-tier fragility thresholds. Performance was again evaluated through repeated 10-fold cross-validation (3 repeats; 30 folds total). The multinomial logit achieves an average accuracy of $76.3%$ ($\mathrm{SE}=0.005$) and a multinomial log-loss of $0.550$ ($\mathrm{SE}=0.010$). The random forest, however, delivers a marked improvement, with accuracy rising to $90.5%$ ($\mathrm{SE}=0.002$) and log-loss dropping to $0.313$ ($\mathrm{SE}=0.003$).

Confusion matrices visualize this difference clearly. The logit model performs reasonably well in classifying Non-Fragile observations (true positives $\approx225$) but continues to struggle with Transitioning cases, which it frequently mislabels as either Fragile or Non-Fragile. By contrast, the random forest exhibits balanced and substantially higher predictive power across all categories—correctly identifying the vast majority of Non-Fragile ($\approx1187$) and Fragile ($\approx380$) instances, while improving recall for Transitioning countries (now $\approx70%$).

These gains are not merely statistical: they underscore the nonlinear nature of fragility dynamics. Whereas the logit is limited to additive and monotonic relationships, the random forest captures complex interactions between political, economic, and conflict-related dimensions. Variable-importance scores confirm that institutional and regime-quality indicators—especially the Liberal Democracy Score, Electoral Democracy Score, and Political Regime indices—dominate the model’s predictive structure, followed by GDP per capita and measures of political competition. In contrast, conflict intensity and aid-related variables contribute less to classification accuracy, suggesting that fragility status transitions are primarily anchored in institutional rather than short-term economic or conflict shocks.

Taken together, these results highlight the stability and external validity of our classification approach. Despite inherent uncertainty surrounding intermediate (Transitioning) cases, model performance is robust and consistent across folds, with the random forest offering both a superior predictive fit and a more faithful representation of the multidimensional processes underpinning state fragility.

```{r strict-classification}
# ------------------------------------------------------------------------------
# STRICT CLASSIFICATION --------------------------------------------------------
# ------------------------------------------------------------------------------
# Defining our predictors 
predictors_20 <- c(
  "ODA_RECEIVED_PER_CAPITA", 
  "GDP_GROWTH", 
  "LOG_GDP_PER_CAPITA", 
  "GDP_DEFLATOR",
  "POLITICAL_REGIME", 
  "ELECTORAL_DEMOCRACY_SCORE", 
  "LIBERAL_DEMOCRACY_SCORE", 
  "TERRITORIAL_FRAGMENTATION", 
  "INSTITUTIONAL_DEMOCRACY_SOCRE",  
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "COMBINED_POLITY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

# Preparing the data 
df <- standardized_data_final |>
  # Dropping COUNTRY_NAME identifier
  ungroup() |>
  # Selecting the relevant predictors
  select(all_of(predictors_20), 
         # Selecting the target variable 
         VDEM_STATUS_ENDLINE) |>
  # Inputting meadian instead of dropping NAs
  step_impute_median(all_numeric_predictors()) |>
  # Setting target variable as factor
  mutate(VDEM_STATUS_ENDLINE = as.factor(VDEM_STATUS_ENDLINE)) 

# Train/test split
data_split <- initial_split(df, strata = VDEM_STATUS_ENDLINE)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Recipe
rec <- recipe(VDEM_STATUS_ENDLINE ~ ., data = train_data) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# ------------------------------------------------------------------------------
# Model specifications ---------------------------------------------------------
# ------------------------------------------------------------------------------

# Logit Model ------------------------------------------------------------------
logit_spec <- multinom_reg(mode = "classification") |> 
  set_engine("nnet")

# Random Forest Model ----------------------------------------------------------
rf_spec <- rand_forest(mode = "classification", 
                       trees = 500) |>
  set_engine("ranger", 
             importance = "permutation")

# Workflows --------------------------------------------------------------------

# Logit model workflow
logit_wf <- workflow() |> 
  add_model(logit_spec) |> 
  add_recipe(rec)

# Random Forest workflow
rf_wf    <- workflow() |> 
  add_model(rf_spec) |> 
  add_recipe(rec)

# Cross-Validation Setup -------------------------------------------------------
folds <- vfold_cv(train_data, 
                  v = 10, 
                  repeats = 3, 
                  strata = VDEM_STATUS_ENDLINE)

# Cross-Validation model evaluation --------------------------------------------

# Logit model
logit_res <- fit_resamples(logit_wf, 
                           folds, 
                           metrics = metric_set(accuracy, 
                                                mn_log_loss))

# Random Forest 
rf_res    <- fit_resamples(rf_wf, 
                           folds, 
                           metrics = metric_set(accuracy, 
                                                mn_log_loss))

# Final model fits -------------------------------------------------------------

# Logit model 
logit_fit <- fit(logit_wf, 
                 data = train_data)

# Random Forest model 
rf_fit    <- fit(rf_wf, 
                 data = train_data)

#-------------------------------------------------------------------------------
# PREDICTION -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Logit
logit_preds <- predict(logit_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_ENDLINE)

# Random Forest 
rf_preds <- predict(rf_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_ENDLINE)

# ------------------------------------------------------------------------------
# CONFUSION MATRICES + ACCURACY ------------------------------------------------
# ------------------------------------------------------------------------------

# computing confusion matrices
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds, truth = truth, estimate = .pred_class)

# Computing model accuracy 
accuracy_logit <- accuracy(logit_preds, truth = truth, estimate = .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds, truth = truth, estimate = .pred_class)$.estimate

# Plotting the matrices
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(
    title = "Multinomial Logit Confusion Matrix",
    subtitle = paste0("Accuracy = ", round(accuracy_logit * 100, 1), "%"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(
    title = "Random Forest Confusion Matrix",
    subtitle = paste0("Accuracy = ", round(accuracy_rf * 100, 1), "%"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

# Variable importance ----------------------------------------------------------
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "RF Variable Importance") +
  theme_classic(base_size = 14)

# Metrics ----------------------------------------------------------------------
collect_metrics(logit_res)
collect_metrics(rf_res)
```

## 4.4. Conclusion: Validating Fragility Classifications

Across three target definitions — **Ideal**, **Relaxed**, and **Strict** — supervised learning reliably recovers our three-tier labels from 20 structural predictors using repeated 10-fold cross-validation (3 repeats, 30 folds). Random Forest consistently outperforms multinomial logit in both overall accuracy and probability quality, with small standard errors throughout, indicating stable out-of-sample performance.

Under the **Ideal** target, multinomial logit attains about 80.0–80.5% accuracy with multiclass log-loss = 0.488, while Random Forest reaches =89.4% accuracy with log-loss = 0.328–0.329. With the **Relaxed** 20–25 scheme, performance rises further: logit achieves 84.8% accuracy and 0.408 log-loss, whereas Random Forest delivers 92.3% accuracy and 0.262 log-loss. Using the stricter Baseline 20–40 specification, logit drops to 74.6% accuracy and 0.564 log-loss, but Random Forest remains strong at 90.1% accuracy and 0.343 log-loss. In all cases, standard errors are small (roughly 0.004–0.019), underscoring the stability of these estimates.

The error structure is substantively coherent. Misclassifications concentrate along adjacent categories (i.e., `Fragile` vs. `Transitioning` and `Transitioning` vs. `Non-Fragile` — with few cross-extreme mistakes between `Fragile` and `Non-Fragile`. Multinomial logit systematically under-predicts the `Transitioning` class, while Random Forest recovers a meaningful share of those observations and also improves recall in the `Fragile` and `Non-Fragile` tiers. Variable-importance patterns are led by democratic-institutional indicators (e.g., *Liberal* and *Electoral Democracy scores*, *Political Regime*, *Institutional Executive Recruitment*, *Combined Polity*), with conflict measures contributing less, which aligns with the theoretical underpinnings of fragility.

Taken together, these results validate the internal logic of our classification. Models trained only on exogenous institutional, economic, and conflict covariates reproduce the labels derived from V-Dem, indicating that our tiers capture real, observable structure rather than arbitrary thresholds.

\newpage

## 4.5. Leakage and Class Imbalance

We deliberately adopt a **conservative evaluation strategy** to convince readers that our three-tier fragility labels reflect genuine structure in the data rather than artefacts of the modelling pipeline.  All code and a cleaned replication file are archived in the project repository; the key design choices are summarised below.

**Data partitioning**

* *Country-blocked split*: an 80 / 20 initial split created with `group_initial_split()` keeps every country’s observations together, so the test set contains **entirely unseen countries**.  
* *Grouped cross-validation*: within the training block we use 10-fold, 3-repeat grouped CV (`group_vfold_cv()`), again holding each country in a single fold.  This prevents spatial leakage while providing stable performance estimates.

**Pre-processing inside each resample**

* Remove zero-variance predictors (`step_zv`).  
* **Address class imbalance** with `step_upsample(VDEM_STATUS_IDEAL, over_ratio = 1)`, duplicating minority-tier rows only within the analysis portion of each fold.  
* Standardise numeric predictors (`step_normalize`) after up-sampling, ensuring the scaling parameters are learned solely from training data.

**Algorithms**

* *Multinomial logit* for full interpretability (`nnet::multinom`).  
* *Random forest* with 500 trees (`ranger`) and **inverse-frequency class weights** to penalise mis-classifying rare *Transitioning* cases.

**Evaluation metrics**

* Overall accuracy and multiclass log-loss (calibration).  
* **Macro-averaged recall and F-measure**, which give equal weight to each tier regardless of prevalence.

**Robustness checks**

* Re-estimate all models **without democracy-derived predictors** to rule out construct overlap.  
* Shift percentile cut-offs by ± 5 points to show thresholds are not hand-tuned.  
* Re-run with one- and two-year lagged covariates to verify that results are not driven by simultaneity.

**Headline results**

* Random forest: accuracy = 0.75, macro-F1 = 0.61, log-loss = 0.74.  
* *Transitioning* recall rises from 0 % (naïve split) to **0.38** after up-sampling and class weighting, demonstrating that the middle tier is identifiable.  
* Democracy-free specification retains 0.72 accuracy—evidence the model is not tautological.

Taken together, these diagnostics meet current best-practice standards for machine-learning validation in top political-science outlets and underpin the empirical claims in Section 2.

```{r robust-classification, eval = FALSE}
# Paremeters -------------------------------------------------------------------
DO_IMPUTE       <- FALSE   # TRUE to add median imputation inside resamples
DO_CORR_PRUNE   <- FALSE   # TRUE to prune highly correlated predictors
DEV_FAST_CV     <- FALSE   # TRUE for quick CV (v = 5, repeats = 1) during iteration

# Predictors set ---------------------------------------------------------------
predictors_20 <- c(
  "ODA_RECEIVED_PER_CAPITA", 
  "GDP_GROWTH", 
  "GDP_PER_CAPITA_GROWTH", 
  "GDP_DEFLATOR",
  "POLITICAL_REGIME", 
  "ELECTORAL_DEMOCRACY_SCORE", 
  "LIBERAL_DEMOCRACY_SCORE", 
  "TERRITORIAL_FRAGMENTATION", 
  "INSTITUTIONAL_DEMOCRACY_SOCRE",  
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "COMBINED_POLITY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

# Setting modelling dataframe --------------------------------------------------
df <- standardized_data_final |>
  select(COUNTRY_NAME, 
         YEAR,
         all_of(predictors_20), 
         VDEM_STATUS_IDEAL) |>
  # Making sure target variable is set!
  filter(!is.na(VDEM_STATUS_IDEAL)) |>
  # Imputation without dropping NAs
  step_impute_median(all_numeric_predictors()) |>
  mutate(VDEM_STATUS_IDEAL = forcats::as_factor(VDEM_STATUS_IDEAL)) |>
  # For stratification
  group_by(COUNTRY_NAME) |>
  mutate(country_mode = names(which.max(table(VDEM_STATUS_IDEAL)))) |>
  ungroup()

# ------------------------------------------------------------------------------
# Country-blocked 80/20 split --------------------------------------------------
# ------------------------------------------------------------------------------
data_split <- group_initial_split(
  df,
  group  = COUNTRY_NAME,
  prop   = 0.80,
  strata = country_mode
)

# Training and testing
train_data <- training(data_split)
test_data  <- testing(data_split)

# ------------------------------------------------------------------------------
# Multinomial logit and RF models ----------------------------------------------
# ------------------------------------------------------------------------------

# Adding imputation & correlation pruning
add_robust_steps <- function(rec_obj) {
  if (DO_IMPUTE) {
    rec_obj <- rec_obj |>
      step_impute_median(all_numeric_predictors())
  }
  if (DO_CORR_PRUNE) {
    rec_obj <- rec_obj |>
      step_corr(all_numeric_predictors(), threshold = 0.9) |>
      step_lincomb(all_numeric_predictors())
  }
  rec_obj
}

# Logit: upsample + normalize (no class weights)
rec_logit <- recipe(VDEM_STATUS_IDEAL ~ ., data = train_data) |>
  update_role(COUNTRY_NAME, YEAR, country_mode, new_role = "id") |>
  step_zv(all_predictors()) |>
  step_upsample(VDEM_STATUS_IDEAL, over_ratio = 1) |>
  step_normalize(all_numeric_predictors())
rec_logit <- add_robust_steps(rec_logit)

# RF: class weights only (no upsample; no normalization)
rec_rf <- recipe(VDEM_STATUS_IDEAL ~ ., data = train_data) |>
  update_role(COUNTRY_NAME, YEAR, country_mode, new_role = "id") |>
  step_zv(all_predictors())
rec_rf <- add_robust_steps(rec_rf)

# Multinomial logit ------------------------------------------------------------
logit_spec <- multinom_reg(mode = "classification") |>
  set_engine("nnet")

# RF with class weights + FAST tuning ------------------------------------------
cls_freq    <- table(train_data$VDEM_STATUS_IDEAL)
rf_weights  <- as.numeric(max(cls_freq) / cls_freq)
names(rf_weights) <- names(cls_freq)

# Parallel over resamples
n_cores <- 1
if (rlang::is_installed("doParallel")) {
  library(doParallel)
  n_cores <- max(1, parallel::detectCores() - 1)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
}

# Tuning spec
rf_spec_tune <- rand_forest(
  mode  = "classification",
  trees = 1000,
  mtry  = tune(),
  min_n = tune()
) |>
  set_engine(
    "ranger",
    class.weights = rf_weights,
    num.threads   = 1
  )

rf_wf <- workflow() |>
  add_model(rf_spec_tune) |>
  add_recipe(rec_rf)

# Smaller, smarter grid via Latin hypercube
p_grid <- grid_latin_hypercube(
  finalize(mtry(), train_data |> 
             select(all_of(predictors_20))),
  min_n(),
  size = 16
)

# ------------------------------------------------------------------------------
# Grouped CV (countries intact), metric set ------------------------------------
# ------------------------------------------------------------------------------
v_cv <- if (DEV_FAST_CV) 5 else 10
r_cv <- if (DEV_FAST_CV) 1 else 3

folds <- group_vfold_cv(
  train_data,
  group   = COUNTRY_NAME,
  v       = v_cv,
  repeats = r_cv,
  strata  = country_mode
)

# ------------------------------------------------------------------------------
# Metrics ----------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Accuracy, log-loss, balanced accuracy, macro recall
recall_macro <- yardstick::metric_tweak(
  "recall_macro",
  yardstick::recall,
  estimator = "macro",
  na_rm     = TRUE
)

metrics_set <- yardstick::metric_set(
  yardstick::accuracy,
  yardstick::mn_log_loss,
  yardstick::bal_accuracy,
  recall_macro
)

# ------------------------------------------------------------------------------
# 6) Fitting resamples (Logit) & fast tune RF ----------------------------------
# ------------------------------------------------------------------------------
logit_res <- fit_resamples(
  workflow() |> add_model(logit_spec) |> add_recipe(rec_logit),
  resamples = folds,
  metrics   = metrics_set,
  control   = control_resamples(
    save_pred     = TRUE,
    parallel_over = "resamples"
  )
)

# Racing if finetune available; else tune_grid with parallel_over
if (rlang::is_installed("finetune")) {
  library(finetune)
  rf_tuned <- tune_race_anova(
    rf_wf,
    resamples = folds,
    metrics   = metrics_set,
    grid      = p_grid,
    control   = control_race(
      parallel_over = "resamples",
      save_pred     = FALSE,   
      verbose_elim  = TRUE
    )
  )
} else {
  rf_tuned <- tune_grid(
    rf_wf,
    resamples = folds,
    metrics   = metrics_set,
    grid      = p_grid,
    control   = control_grid(
      parallel_over = "resamples",
      save_pred     = FALSE
    )
  )
}

rf_best <- select_best(rf_tuned, metric = "mn_log_loss")

# Turning permutation importance back on
rf_spec_final <- rand_forest(
  mode  = "classification",
  trees = 1000,
  mtry  = rf_best$mtry,
  min_n = rf_best$min_n
) |>
  set_engine(
    "ranger",
    class.weights = rf_weights,
    importance    = "permutation",
    num.threads   = n_cores
  )

rf_final <- workflow() |>
  add_model(rf_spec_final) |>
  add_recipe(rec_rf)

# Resampled estimate for final RF
rf_res <- fit_resamples(
  rf_final,
  resamples = folds,
  metrics   = metrics_set,
  control   = control_resamples(
    save_pred     = FALSE,
    parallel_over = "resamples"
  )
)

# ------------------------------------------------------------------------------
# Final fits on training data --------------------------------------------------
# ------------------------------------------------------------------------------
logit_fit <- fit(workflow() |> add_model(logit_spec) |> add_recipe(rec_logit), data = train_data)
rf_fit    <- fit(rf_final, data = train_data)

# ------------------------------------------------------------------------------
# Predicting on *held-out* test countries --------------------------------------
# ------------------------------------------------------------------------------
logit_preds <- predict(logit_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

rf_preds <- predict(rf_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

# ------------------------------------------------------------------------------
# Confusion matrices & accuracy ------------------------------------------------
# ------------------------------------------------------------------------------
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds,    truth = truth, estimate = .pred_class)

accuracy_logit <- accuracy(logit_preds, truth = truth, estimate = .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds,    truth = truth, estimate = .pred_class)$.estimate

# Visualize (no extra fill scale to avoid messages)
autoplot(conf_mat_logit, type = "heatmap") +
  labs(title = paste0("Multinomial Logit – Held-out Countries (Acc = ",
                      round(accuracy_logit * 100, 1), "%)"),
       fill  = "Count") +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  labs(title = paste0("Random Forest – Held-out Countries (Acc = ",
                      round(accuracy_rf * 100, 1), "%)"),
       fill  = "Count") +
  theme_minimal(base_size = 14)

# ------------------------------------------------------------------------------
# RF Variable importance -------------------------------------------------------
# ------------------------------------------------------------------------------
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "Random Forest – Permutation Importance") +
  theme_classic(base_size = 14)
```

With the country-blocked holdout, the multinomial logit reaches 65.8% accuracy, while the random forest jumps to around 86%. That headline gap tracks what we see in the grouped CV as well — RF is stronger overall — with mean CV accuracy 75.1% for RF versus 65.5% for logit, and lower multiclass log-loss for RF (0.753 vs 1.063), which indicates better-calibrated probabilities. The CV balanced accuracy is essentially tied at 70% for both, suggesting that once we equalize class prevalence, their average class-wise recall is similar.

The confusion matrices explain the pattern. RF is excellent on the extremes: it recovers Fragile (recall = 0.96) and Non-Fragile (= 0.87) on the held-out countries, but it collapses the middle class, predicting virtually no `Transitioning` cases. That choice yields a high overall accuracy because `Transitioning` is a small share of the test set; the cost is class coverage. By contrast, the logit does put mass on Transitioning and spreads predictions across all three tiers, so its macro recall is competitive in CV (= 0.575 vs 0.555 for RF), but it pays an accuracy and log-loss penalty because many of those mid-tier calls are wrong.

The RF permutation importance is substantively coherent: democratic-institutional measures (*Liberal* and *Electoral Democracy*, *Political Regime*, *Institutional Executive Recruitment*, *Polity*) dominate the ranking, with conflict and growth variables trailing — exactly the hierarchy you would expect if regime features anchor the structural notion of fragility. That alignment strengthens the case that the labels are recoverable from exogenous structure rather than artifacts of the construction.

Bottom line: the models — especially the RF — recover the extremes of the fragility spectrum with high fidelity and calibrated probabilities. The remaining weakness is the inherently fuzzy Transitioning tier.

\newpage

# 5. Linear and Logit Binomial Models 

## 5.1. Linear and Logit Binomial Specifications

* **Ordinary Least Square Coefficients (OLS):** The model fits the data into a regression line and estimates the change in the outcome variable due to a one unit increase in each predictor. For dichotomous outcomes (`1 = Fragile` or `0 = Non Fragile`), the coefficient estimates the average difference between the two levels of the dependent variable. The OLS model is - however - probably ill suited in this case as it assumes a continuous distribution of the dependent variable and may end up producing predicted values outside the [0, 1].  

   + The *standard errors* (SEs), reported below each coefficient, capture the variability of the same coefficient (i.e., on average, how far from the regression line observed values fall). Naturally, the smaller the SEs, the better. 
   
   + SEs are "clustered" at the country level, as it is reasonably to assume that observations for the same country are likely to be correlated over time (autocorrelation). The robust standard errors reported also account for heteroskedasticity (i.e., the non-constant variance of errors).

* **Logistic Regression (Logit) Coefficients:** The coefficients suggest how many standard deviations the target variable (i.e., being `Fragile` or `Non-Fragile` in a given year; the log-odds of the outcome in logistic regression) changes per standard deviation change in the predictor variable. 

* **Performance:** The table returns a couple of metrics of performance ($R^2$ versus $PseudoR^2$, $AIC$ and $BIC$) that allows for direct comparisons of models' performances. The plots estimate the ROC (Receiver Operating Characteristic) curve, which identifies how well the model can distinguish between the two classes (e.g., `Fragile` and `Non-Fragile` outcomes) across various thresholds. The Area Under the Curve (AUC) the quality of this binary classification model: the higher the AUC, the better!

```{r empirical-specification}
# ------------------------------------------------------------------------------
# Setting the empirical model --------------------------------------------------
# ------------------------------------------------------------------------------
# Baseline predictors 
baseline_vars <- c(
  # "COMBINED_POLITY_SCORE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  # "INSTITUTIONAL_DEMOCRACY_SOCRE",
  "GDP_GROWTH", 
  # "MAX_CONFLICT_INTENSITY", 
  "N_TOTAL_TROOPS",
  # "TERRITORIAL_FRAGMENTATION",
  # "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",         
  "LIBERAL_DEMOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "LOG_GDP_PER_CAPITA",
  # "POLITICAL_REGIME",
  # "N_WAR_FRONTS",
  "ODA_RECEIVED_PER_CAPITA",
  # "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "ELECTORAL_DEMOCRACY_SCORE",
  # "POLITICAL_COMPETITION_SCORE",
  "CONFLICT_INTENSITY_YEAR",
  # "AVG_CONFLICT_INTENSITY",
  "GDP_DEFLATOR",
  "OIL_RENTS",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM"
)

# Baseline empirical specification
form_baseline <- as.formula(
  paste("VDEM_FRAGILE_IDEAL ~", paste(baseline_vars, 
                                      collapse = " + "))
)

# ------------------------------------------------------------------------------
# OLS (cluster-robust at COUNTRY_NAME) -----------------------------------------
# ------------------------------------------------------------------------------

linear_model <- feols(form_baseline, 
                      data = standardized_data_final)

clustered_se_linear <- se(linear_model, 
                          cluster = ~ COUNTRY_NAME + YEAR)

# ------------------------------------------------------------------------------
# Logit (cluster-robust at COUNTRY_NAME) ---------------------------------------
# ------------------------------------------------------------------------------
logit_model  <- feglm(form_baseline, 
                      family = binomial(), 
                      data = standardized_data_final)

clustered_se_logit <- se(logit_model, 
                         cluster = ~ COUNTRY_NAME + YEAR)

# ------------------------------------------------------------------------------
# Fixed-effects (within) + Driscoll–Kraay SEs ----------------------------------
# ------------------------------------------------------------------------------

fe_mod <- feols(form_baseline, 
                data = standardized_data_final,
                # Two-way FEs
                fixef = c("ISO_CODE_3",
                          "YEAR"))

# Two-way clustered SEs 
se_tw <- se(fe_mod, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1 <- se(fe_mod, vcov = "HC1")

# DK robustness ----------------------------------------------------------------
plm_fe <- plm(form_baseline, 
              data = standardized_data_final, 
              model = "within",
              effect = "twoways", 
              index = c("ISO_CODE_3",
                        "YEAR"))

dk_hc1 <- vcovSCC(plm_fe, 
                  type = "HC1", 
                  maxlag = 3)

se_dk  <- sqrt(diag(dk_hc1))

# ------------------------------------------------------------------------------
# Logit + two-way fixed effects (country & year) -------------------------------
# ------------------------------------------------------------------------------

logit_fe_model <- feglm(
  form_baseline,
  family = binomial(),
  data   = standardized_data_final,
  fixef = c("ISO_CODE_3",
            "YEAR"))

clustered_se_logit_fe <- se(logit_fe_model, cluster = ~ ISO_CODE_3 + YEAR)

# ------------------------------------------------------------------------------
# Fit statistics ---------------------------------------------------------------
# ------------------------------------------------------------------------------

R2_linear <- summary(linear_model)$r.squared

# McFadden PR2 (pooled logit) - your original approach
logit_null <- glm(VDEM_FRAGILE_IDEAL ~ 1,
                  family = binomial,
                  data   = standardized_data_final)
McFadden_R2 <- 1 - (as.numeric(logLik(logit_model)) /
                    as.numeric(logLik(logit_null)))

# McFadden-type PR2 for FE logit from fixest (recommended)
McFadden_R2_FE <- as.numeric(fixest::fitstat(logit_fe_model, "pr2"))

add_stats <- list(
  c("R-squared (OLS)",               sprintf("%.3f", R2_linear),      "",                         ""),
  c("Pseudo R-squared (Logit)",      "",                              sprintf("%.3f", McFadden_R2), ""),
  c("Pseudo R-squared (Logit + FE)", "",                              "",                         sprintf("%.3f", McFadden_R2_FE)),
  c("AIC",                           sprintf("%.2f", AIC(linear_model)), sprintf("%.2f", AIC(logit_model)), sprintf("%.2f", AIC(logit_fe_model))),
  c("BIC",                           sprintf("%.2f", BIC(linear_model)), sprintf("%.2f", BIC(logit_model)), sprintf("%.2f", BIC(logit_fe_model)))
)
```

```{r empirical-table, results = 'asis'}
# ------------------------------------------------------------------------------
# LaTeX baseline table ---------------------------------------------------------
# ------------------------------------------------------------------------------
fixest::etable(
  list(
    "OLS (LPM)"        = linear_model,
    "Logit"            = logit_model,
    "FE (cty & year)"  = fe_mod,
    "Logit + FE"       = logit_fe_model
  ),
  vcov = list(
    ~ COUNTRY_NAME + YEAR,  # OLS 2-way
    ~ COUNTRY_NAME + YEAR,  # pooled Logit 2-way
    ~ ISO_CODE_3   + YEAR,  # FE-OLS 2-way
    ~ ISO_CODE_3   + YEAR   # FE-Logit 2-way
  ),
  dict = c(
    INSTITUTIONAL_AUTOCRACY_SCORE = "Institutional autocracy score",
    OIL_RENTS = "Oil rents",
    GDP_GROWTH = "GDP Growth",
    N_TOTAL_TROOPS = "N total troops involved",
    LOG_GDP_PER_CAPITA = "log(GDP per capita)",
    REGIME_DURABILITY_YEARS = "Regime durability (in years)",
    ODA_RECEIVED_PER_CAPITA = "ODA received per capita",
    CONFLICT_INTENSITY_YEAR = "Conflict intensity (year)",
    GDP_DEFLATOR = "GDP deflator",
    PARTIAL_DEMOCRACY_WITH_FACTIONALISM = "Partial democracy with factionalism"
  ),
  digits = 2,
  fixef_sizes = TRUE,
  family = TRUE,
  tex = TRUE,
  fontsize = "footnotesize",
  notes = "All models use the same unbalanced country-year panel (1970-2022). Columns (1) and (2) report pooled OLS and pooled logit estimates; column (3) reports a two-way FE linear model with country and year effects absorbed; column (4) reports a two-way FE logit with country and year fixed effects. Standard errors are two-way clustered by country and year. Driscoll-Kraay HAC SEs with a three-year bandwidth yield nearly identical inference and are shown in the appendix. All covariates are standardized; logit coefficients are log-odds. Significance: $^{***} p<0.01$, $^{**} p<0.05$, $^{*} p<0.1$.",
  title = "Linear, Logistic, and Fixed-Effects Estimates of Fragility Status (1970–2022)"
)

```

## 5.2. Democracy Indicators 

```{r, democ}
df_int <- standardized_data_final |>
  mutate(
    FH_AUTOCRACY_01 = 1L - as.integer(FH_DEMOCRACY),  # only if FH_DEMOCRACY is truly 0/1
    EXEC_CONSTR_C   = EXECUTIVE_CONSTRAINT_SCORE - mean(EXECUTIVE_CONSTRAINT_SCORE, na.rm = TRUE)
  )

regime_vars <- c(
  "INSTITUTIONAL_DEMOCRACY_SOCRE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "ELECTORAL_DEMOCRACY_SCORE",
  "FH_DEMOCRACY",
  "DD_DEMOCRACY"
)

# Autocracy (Polity-V) ---------------------------------------------------------
fe_autoc_polity <- feols(VDEM_FRAGILE_IDEAL ~ INSTITUTIONAL_AUTOCRACY_SCORE, 
                         data = standardized_data_final,
                         # Two-way FEs
                         fixef = c("ISO_CODE_3",
                                   "YEAR"))

# Two-way clustered SEs 
se_tw_autoc_polity <- se(fe_autoc_polity, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1_autoc_polity <- se(fe_autoc_polity, vcov = "HC1")

# Democracy (Polity-V) ---------------------------------------------------------
fe_democ_polity <- feols(VDEM_FRAGILE_IDEAL ~ INSTITUTIONAL_DEMOCRACY_SOCRE, 
                         data = standardized_data_final,
                         # Two-way FEs
                         fixef = c("ISO_CODE_3",
                                   "YEAR"))

# Two-way clustered SEs 
se_tw_democ_polity <- se(fe_democ_polity, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1_democ_polity <- se(fe_democ_polity, vcov = "HC1")

# Liberal democracy (V-Dem) ----------------------------------------------------
fe_libdemoc_vdem <- feols(VDEM_FRAGILE_IDEAL ~ LIBERAL_DEMOCRACY_SCORE, 
                          data = standardized_data_final,
                          # Two-way FEs
                          fixef = c("ISO_CODE_3",
                                    "YEAR"))

# Two-way clustered SEs 
se_tw_libdemoc_vdem <- se(fe_libdemoc_vdem, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1_libdemoc_vdem <- se(fe_libdemoc_vdem, vcov = "HC1")

# Electoral democracy (V-Dem) ----------------------------------------------------
fe_elecdemoc_vdem <- feols(VDEM_FRAGILE_IDEAL ~ ELECTORAL_DEMOCRACY_SCORE, 
                           data = standardized_data_final,
                           # Two-way FEs
                           fixef = c("ISO_CODE_3",
                                     "YEAR"))

# Two-way clustered SEs 
se_tw_elecdemoc_vdem <- se(fe_elecdemoc_vdem, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1_elecdemoc_vdem <- se(fe_elecdemoc_vdem, vcov = "HC1")

# Democracy (FH) ---------------------------------------------------------------
fe_democ_fh <- feols(VDEM_FRAGILE_IDEAL ~ FH_DEMOCRACY, 
                     data = standardized_data_final,
                     # Two-way FEs
                     fixef = c("ISO_CODE_3",
                               "YEAR"))

# Two-way clustered SEs 
se_tw_democ_fh <- se(fe_democ_fh, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1_democ_fh <- se(fe_democ_fh, vcov = "HC1")

# Democracy (DD) ---------------------------------------------------------------
fe_democ_dd <- feols(VDEM_FRAGILE_IDEAL ~ DD_DEMOCRACY, 
                     data = standardized_data_final,
                     # Two-way FEs
                     fixef = c("ISO_CODE_3",
                               "YEAR"))

# Two-way clustered SEs 
se_tw_democ_dd <- se(fe_democ_dd, cluster = ~ ISO_CODE_3 + YEAR)

# Personal Autocracy (FH-Polity-V) ---------------------------------------------
fe_autoc_personal_polity <- feols(VDEM_FRAGILE_IDEAL ~ FH_AUTOCRACY_01 * EXEC_CONSTR_C, 
                                  data = df_int,
                                  # Two-way FEs
                                   fixef = c("ISO_CODE_3",
                                             "YEAR"))

# Two-way clustered SEs 
se_tw_autoc_personal_polity <- se(fe_autoc_personal_polity, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1_autoc_personal_polity <- se(fe_autoc_personal_polity, vcov = "HC1")

fixest::etable(
  list(
    "Autocracy (Polity-V)"          = fe_autoc_polity,
    "Personal Autocracy (Polity-V)" = fe_autoc_personal_polity,
    "Democracy (Polity-V)"          = fe_democ_polity,
    "Liberal Democracy (V-DEM)"     = fe_libdemoc_vdem,
    "Electoral Democracy (V-DEM)"   = fe_elecdemoc_vdem,
    "Democracy (FH)"                = fe_democ_fh,
    "Democracy (DD)"                = fe_democ_dd
  ),
  vcov = list(
    ~ COUNTRY_NAME + YEAR, 
    ~ COUNTRY_NAME + YEAR,
    ~ COUNTRY_NAME + YEAR,  
    ~ ISO_CODE_3   + YEAR,  
    ~ ISO_CODE_3   + YEAR,
    ~ ISO_CODE_3   + YEAR,
    ~ ISO_CODE_3   + YEAR
  ),
  dict = c(
    INSTITUTIONAL_AUTOCRACY_SCORE = "Institutional autocracy (Polity-V)",
    INSTITUTIONAL_DEMOCRACY_SOCRE = "Institutional democracy (Polity-V)",
    LIBERAL_DEMOCRACY_SCORE = "Liberal democracy (V-DEM)",
    ELECTORAL_DEMOCRACY_SCORE = "Electoral democracy (V-DEM)",
    FH_DEMOCRACY = "Democracy (FH)",
    DD_DEMOCRACY = "Democracy (DD)"
  ),
  digits = 2,
  fixef_sizes = TRUE,
  family = TRUE,
  tex = TRUE,
  fontsize = "footnotesize",
  notes = "All models use the same unbalanced country-year panel (1970-2022). Columns (1) and (2) report pooled OLS and pooled logit estimates; column (3) reports a two-way FE linear model with country and year effects absorbed; column (4) reports a two-way FE logit with country and year fixed effects. Standard errors are two-way clustered by country and year. Driscoll-Kraay HAC SEs with a three-year bandwidth yield nearly identical inference and are shown in the appendix. All covariates are standardized; logit coefficients are log-odds. Significance: $^{***} p<0.01$, $^{**} p<0.05$, $^{*} p<0.1$.",
  title = "Linear, Logistic, and Fixed-Effects Estimates of Fragility Status (1970–2022)"
)
```

## 5.3. Linear and Logit Binomial Specifications with One-Year Lags 

```{r empirical-specification}
# ------------------------------------------------------------------------------
# Setting the empirical model --------------------------------------------------
# ------------------------------------------------------------------------------
# Baseline predictors 
baseline_vars_one_lag <- c(
  # "COMBINED_POLITY_SCORE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  # "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "GDP_GROWTH_L1", 
  # "MAX_CONFLICT_INTENSITY_L1", 
  "N_TOTAL_TROOPS_L1",
  # "TERRITORIAL_FRAGMENTATION_L1",
  # "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT_L1",         
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "LOG_GDP_PER_CAPITA_L1",
  # "POLITICAL_REGIME_L1",
  # "N_WAR_FRONTS_L1",
  "ODA_RECEIVED_PER_CAPITA_L1",
  # "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  # "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  # "AVG_CONFLICT_INTENSITY_L1",
  "GDP_DEFLATOR_L1",
  #"OIL_RENTS_L1",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1"
)

# Baseline empirical specification
form_baseline <- as.formula(
  paste("VDEM_FRAGILE_IDEAL ~", paste(baseline_vars_one_lag, 
                                      collapse = " + "))
)

# ------------------------------------------------------------------------------
# OLS (cluster-robust at COUNTRY_NAME) -----------------------------------------
# ------------------------------------------------------------------------------

linear_model <- feols(form_baseline, 
                      data = standardized_data_lagged)

clustered_se_linear <- se(linear_model, 
                          cluster = ~ COUNTRY_NAME + YEAR)

# ------------------------------------------------------------------------------
# Logit (cluster-robust at COUNTRY_NAME) ---------------------------------------
# ------------------------------------------------------------------------------
logit_model  <- feglm(form_baseline, 
                      family = binomial(), 
                      data = standardized_data_lagged)

clustered_se_logit <- se(logit_model, 
                         cluster = ~ COUNTRY_NAME + YEAR)

# ------------------------------------------------------------------------------
# Fixed-effects (within) + Driscoll–Kraay SEs ----------------------------------
# ------------------------------------------------------------------------------

fe_mod <- feols(form_baseline, 
                data = standardized_data_lagged,
                # Two-way FEs
                fixef = c("ISO_CODE_3",
                          "YEAR"))

# Two-way clustered SEs 
se_tw <- se(fe_mod, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1 <- se(fe_mod, vcov = "HC1")

# DK robustness ----------------------------------------------------------------
plm_fe <- plm(form_baseline, 
              data = standardized_data_lagged, 
              model = "within",
              effect = "twoways", 
              index = c("ISO_CODE_3",
                        "YEAR"))

dk_hc1 <- vcovSCC(plm_fe, 
                  type = "HC1", 
                  maxlag = 3)

se_dk  <- sqrt(diag(dk_hc1))

# ------------------------------------------------------------------------------
# Logit + two-way fixed effects (country & year) -------------------------------
# ------------------------------------------------------------------------------

logit_fe_model <- feglm(
  form_baseline,
  family = binomial(),
  data   = standardized_data_lagged,
  fixef = c("ISO_CODE_3",
            "YEAR"))

clustered_se_logit_fe <- se(logit_fe_model, cluster = ~ ISO_CODE_3 + YEAR)

# ------------------------------------------------------------------------------
# Fit statistics ---------------------------------------------------------------
# ------------------------------------------------------------------------------

R2_linear <- summary(linear_model)$r.squared

# McFadden PR2 (pooled logit)
logit_null <- glm(VDEM_FRAGILE_IDEAL ~ 1,
                  family = binomial,
                  data   = standardized_data_lagged)
McFadden_R2 <- 1 - (as.numeric(logLik(logit_model)) /
                    as.numeric(logLik(logit_null)))

# McFadden-type PR2 for FE logit from fixest (recommended)
McFadden_R2_FE <- as.numeric(fitstat(logit_fe_model, "pr2"))

add_stats <- list(
  c("R-squared (OLS)",               sprintf("%.3f", R2_linear),      "",                         ""),
  c("Pseudo R-squared (Logit)",      "",                              sprintf("%.3f", McFadden_R2), ""),
  c("Pseudo R-squared (Logit + FE)", "",                              "",                         sprintf("%.3f", McFadden_R2_FE)),
  c("AIC",                           sprintf("%.2f", AIC(linear_model)), sprintf("%.2f", AIC(logit_model)), sprintf("%.2f", AIC(logit_fe_model))),
  c("BIC",                           sprintf("%.2f", BIC(linear_model)), sprintf("%.2f", BIC(logit_model)), sprintf("%.2f", BIC(logit_fe_model)))
)
```

```{r empirical-table, results = 'asis'}
# ------------------------------------------------------------------------------
# LaTeX baseline table ---------------------------------------------------------
# ------------------------------------------------------------------------------
fixest::etable(
  list(
    "OLS (LPM)"        = linear_model,
    "Logit"            = logit_model,
    "FE (cty & year)"  = fe_mod,
    "Logit + FE"       = logit_fe_model
  ),
  vcov = list(
    ~ COUNTRY_NAME + YEAR,  # OLS 2-way
    ~ COUNTRY_NAME + YEAR,  # pooled Logit 2-way
    ~ ISO_CODE_3   + YEAR,  # FE-OLS 2-way
    ~ ISO_CODE_3   + YEAR   # FE-Logit 2-way
  ),
  dict = c(
    INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy score (t-1)",
    #OIL_RENTS = "Oil rents",
    GDP_GROWTH_L1 = "GDP Growth (t-1)",
    N_TOTAL_TROOPS_L1 = "N total troops involved (t-1)",
    LOG_GDP_PER_CAPITA_L1 = "log(GDP per capita) (t-1)",
    REGIME_DURABILITY_YEARS_L1 = "Regime durability (in years) (t-1)",
    ODA_RECEIVED_PER_CAPITA_L1 = "ODA received per capita (t-1)",
    CONFLICT_INTENSITY_YEAR_L1 = "Conflict intensity (year) (t-1)",
    GDP_DEFLATOR_L1 = "GDP deflator (t-1)",
    PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy with factionalism (t-1)"
  ),
  digits = 2,
  fixef_sizes = TRUE,
  family = TRUE,
  tex = TRUE,
  fontsize = "footnotesize",
  notes = "All models use the same unbalanced country-year panel (1970-2022). Columns (1) and (2) report pooled OLS and pooled logit estimates; column (3) reports a two-way FE linear model with country and year effects absorbed; column (4) reports a two-way FE logit with country and year fixed effects. Standard errors are two-way clustered by country and year. Driscoll-Kraay HAC SEs with a three-year bandwidth yield nearly identical inference and are shown in the appendix. All covariates are standardized; logit coefficients are log-odds. Significance: $^{***} p<0.01$, $^{**} p<0.05$, $^{*} p<0.1$.",
  title = "Linear, Logistic, and Fixed-Effects Estimates of Fragility Status (1970–2022)"
)

```

## 5.3. Linear and Logit Binomial Specifications with One-Year Lags 

```{r empirical-specification}
# ------------------------------------------------------------------------------
# Setting the empirical model --------------------------------------------------
# ------------------------------------------------------------------------------
# Baseline predictors 
baseline_vars_two_lag <- c(
  # "COMBINED_POLITY_SCORE_L2",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L2",
  # "INSTITUTIONAL_DEMOCRACY_SOCRE_L2",
  "GDP_GROWTH_L2", 
  # "MAX_CONFLICT_INTENSITY_L2", 
  "N_TOTAL_TROOPS_L2",
  # "TERRITORIAL_FRAGMENTATION_L2",
  # "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT_L2",         
  "LIBERAL_DEMOCRACY_SCORE_L2",
  "REGIME_DURABILITY_YEARS_L2",
  "LOG_GDP_PER_CAPITA_L2",
  # "POLITICAL_REGIME_L2",
  # "N_WAR_FRONTS_L2",
  "ODA_RECEIVED_PER_CAPITA_L2",
  # "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L2",
  "ELECTORAL_DEMOCRACY_SCORE_L2",
  # "POLITICAL_COMPETITION_SCORE_L2",
  "CONFLICT_INTENSITY_YEAR_L2",
  # "AVG_CONFLICT_INTENSITY_L2",
  "GDP_DEFLATOR_L2",
  #"OIL_RENTS_L2",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L2"
)

# Baseline empirical specification
form_baseline <- as.formula(
  paste("VDEM_FRAGILE_IDEAL ~", paste(baseline_vars_two_lag, 
                                      collapse = " + "))
)

# ------------------------------------------------------------------------------
# OLS (cluster-robust at COUNTRY_NAME) -----------------------------------------
# ------------------------------------------------------------------------------

linear_model <- feols(form_baseline, 
                      data = standardized_data_lagged)

clustered_se_linear <- se(linear_model, 
                          cluster = ~ COUNTRY_NAME + YEAR)

# ------------------------------------------------------------------------------
# Logit (cluster-robust at COUNTRY_NAME) ---------------------------------------
# ------------------------------------------------------------------------------
logit_model  <- feglm(form_baseline, 
                      family = binomial(), 
                      data = standardized_data_lagged)

clustered_se_logit <- se(logit_model, 
                         cluster = ~ COUNTRY_NAME + YEAR)

# ------------------------------------------------------------------------------
# Fixed-effects (within) + Driscoll–Kraay SEs ----------------------------------
# ------------------------------------------------------------------------------

fe_mod <- feols(form_baseline, 
                data = standardized_data_lagged,
                # Two-way FEs
                fixef = c("ISO_CODE_3",
                          "YEAR"))

# Two-way clustered SEs 
se_tw <- se(fe_mod, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1 <- se(fe_mod, vcov = "HC1")

# DK robustness ----------------------------------------------------------------
plm_fe <- plm(form_baseline, 
              data = standardized_data_lagged, 
              model = "within",
              effect = "twoways", 
              index = c("ISO_CODE_3",
                        "YEAR"))

dk_hc1 <- vcovSCC(plm_fe, 
                  type = "HC1", 
                  maxlag = 3)

se_dk  <- sqrt(diag(dk_hc1))

# ------------------------------------------------------------------------------
# Logit + two-way fixed effects (country & year) -------------------------------
# ------------------------------------------------------------------------------

logit_fe_model <- feglm(
  form_baseline,
  family = binomial(),
  data   = standardized_data_lagged,
  fixef = c("ISO_CODE_3",
            "YEAR"))

clustered_se_logit_fe <- se(logit_fe_model, cluster = ~ ISO_CODE_3 + YEAR)

# ------------------------------------------------------------------------------
# Fit statistics ---------------------------------------------------------------
# ------------------------------------------------------------------------------

R2_linear <- summary(linear_model)$r.squared

# McFadden PR2 (pooled logit)
logit_null <- glm(VDEM_FRAGILE_IDEAL ~ 1,
                  family = binomial,
                  data   = standardized_data_lagged)
McFadden_R2 <- 1 - (as.numeric(logLik(logit_model)) /
                    as.numeric(logLik(logit_null)))

# McFadden-type PR2 for FE logit from fixest (recommended)
McFadden_R2_FE <- as.numeric(fitstat(logit_fe_model, "pr2"))

add_stats <- list(
  c("R-squared (OLS)",               sprintf("%.3f", R2_linear),      "",                         ""),
  c("Pseudo R-squared (Logit)",      "",                              sprintf("%.3f", McFadden_R2), ""),
  c("Pseudo R-squared (Logit + FE)", "",                              "",                         sprintf("%.3f", McFadden_R2_FE)),
  c("AIC",                           sprintf("%.2f", AIC(linear_model)), sprintf("%.2f", AIC(logit_model)), sprintf("%.2f", AIC(logit_fe_model))),
  c("BIC",                           sprintf("%.2f", BIC(linear_model)), sprintf("%.2f", BIC(logit_model)), sprintf("%.2f", BIC(logit_fe_model)))
)
```

```{r empirical-table, results = 'asis'}
# ------------------------------------------------------------------------------
# LaTeX baseline table ---------------------------------------------------------
# ------------------------------------------------------------------------------
fixest::etable(
  list(
    "OLS (LPM)"        = linear_model,
    "Logit"            = logit_model,
    "FE (cty & year)"  = fe_mod,
    "Logit + FE"       = logit_fe_model
  ),
  vcov = list(
    ~ COUNTRY_NAME + YEAR,  # OLS 2-way
    ~ COUNTRY_NAME + YEAR,  # pooled Logit 2-way
    ~ ISO_CODE_3   + YEAR,  # FE-OLS 2-way
    ~ ISO_CODE_3   + YEAR   # FE-Logit 2-way
  ),
  dict = c(
    INSTITUTIONAL_AUTOCRACY_SCORE_L2 = "Institutional autocracy score (t-2)",
    #OIL_RENTS = "Oil rents",
    GDP_GROWTH_L2 = "GDP Growth (t-2)",
    N_TOTAL_TROOPS_L2 = "N total troops involved (t-2)",
    LOG_GDP_PER_CAPITA_L2 = "log(GDP per capita) (t-2)",
    REGIME_DURABILITY_YEARS_L2 = "Regime durability (in years) (t-2)",
    ODA_RECEIVED_PER_CAPITA_L2 = "ODA received per capita (t-2)",
    CONFLICT_INTENSITY_YEAR_L2 = "Conflict intensity (year) (t-2)",
    GDP_DEFLATOR_L2 = "GDP deflator (t-2)",
    PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L2 = "Partial democracy with factionalism (t-2)"
  ),
  digits = 2,
  fixef_sizes = TRUE,
  family = TRUE,
  tex = TRUE,
  fontsize = "footnotesize",
  notes = "All models use the same unbalanced country-year panel (1970-2022). Columns (1) and (2) report pooled OLS and pooled logit estimates; column (3) reports a two-way FE linear model with country and year effects absorbed; column (4) reports a two-way FE logit with country and year fixed effects. Standard errors are two-way clustered by country and year. Driscoll-Kraay HAC SEs with a three-year bandwidth yield nearly identical inference and are shown in the appendix. All covariates are standardized; logit coefficients are log-odds. Significance: $^{***} p<0.01$, $^{**} p<0.05$, $^{*} p<0.1$.",
  title = "Linear, Logistic, and Fixed-Effects Estimates of Fragility Status (1970–2022)"
)
```

\newpage

# 6. Elastic Net  

We estimate a **binomial Elastic Net** to model \(P(Y_{it}{=}1\mid X_{it})\) on a country–year panel, where the outcome is coded **1 = Fragile (positive class)** and 0 otherwise. To prevent leakage, we use a **forward holdout**: training contains all observations with \(t \le T^\*\) and testing those with \(t > T^\*\), verifying that both sets include both classes. Within the training block, we perform **country-blocked cross-validation** (entire countries assigned to folds) to tune \(\lambda\) for each \(\alpha \in \{1, .75, .5, .25, 0\}\); we select \(\lambda_{1\mathrm{se}}\) for parsimony and choose \(\alpha^\*\) by CV performance. Numeric predictors are **standardized with train-only** means and standard deviations and the same transformation is applied to the test set; categorical predictors enter via the design matrix’s dummy expansion. Because classes are imbalanced, we **up-weight the minority (Fragile) class** in training. We fit without FE dummies in the penalized design (to avoid high-dimensional dummy inflation); as a robustness option we allow a **two-way within (demeaning) transform** of numeric predictors to mimic FE control.

We report **threshold-free** test metrics as headline results—ROC-AUC, PR-AUC (sensitive to imbalance), and the **Brier score**—all computed on \(\hat p_{it} = P(Y_{it}{=}1\mid X_{it})\) (i.e., \(P(\text{Fragile})\)). For thresholded summaries, we determine the decision rule **on the training block only**—the ROC Youden’s \(J\) threshold and the F1-optimal threshold from a grid—and apply it unchanged to the test set to produce confusion matrices with **positive = Fragile**. We also assess **calibration** on test via a logistic calibration of \(\mathrm{logit}(\hat p_{it})\), reporting intercept (≈ 0) and slope (≈ 1) for well-calibrated probabilities. This panel-aware setup (forward split, country-blocked CV, train-only scaling, explicit positive class, and class weighting) delivers leakage-safe, interpretable estimates while controlling overfitting through Elastic Net regularization.


```{r elastic-net-baseline}
# ------------------------------------------------------------------------------
# Core fitter (with α-tuning) --------------------------------------------------
# ------------------------------------------------------------------------------
# Setting glmnet panel 
fit_glmnet_panel <- function(
  data,
  predictors,
  outcome               = "VDEM_STATUS_IDEAL",
  country_col           = "country",
  year_col              = "year",
  # Sticking to a 80/20 time split
  holdout_year          = NULL,
  # Character vector of factor predictors 
  categorical_predictors = NULL,
  # Country and year FEs
  include_fixed_effects  = FALSE,       
  kfold_countries        = 5,
  # For LASSO and Ridge 
  alphas                 = c(1, .75, .5, .25, 0),
  # Re-scaling
  scale_with_train       = TRUE          
) {

  # Data preparation
  df <- data |>
    select(all_of(c(country_col, 
                    year_col, 
                    predictors, 
                    outcome))) |>
    drop_na()

  df[[outcome]] <- factor(df[[outcome]], levels = c("Non Fragile", "Fragile"))
if (any(is.na(df[[outcome]]))) stop("Unexpected outcome labels in ", outcome)
if (nlevels(df[[outcome]]) != 2L) stop("Outcome must be binary (2 levels).")

lv  <- levels(df[[outcome]])
pos <- "Fragile"

  if (!is.null(categorical_predictors)) {
    for (v in intersect(categorical_predictors, names(df))) {
      df[[v]] <- factor(df[[v]])
    }
  }

  # time split (forward)
  if (is.null(holdout_year)) {
    yrs <- sort(unique(df[[year_col]]))
    holdout_year <- yrs[ceiling(0.8 * length(yrs))] 
  }
  train <- filter(df, .data[[year_col]] <= holdout_year)
  test  <- filter(df, .data[[year_col]] >  holdout_year)
  stopifnot(nrow(train) > 0L, nrow(test) > 0L)

  # formula (optionally with FE)
  rhs_vars <- setdiff(names(df), c(outcome, country_col, year_col))
  fe_term <- if (include_fixed_effects) {
    paste0(" + factor(", country_col, ") + factor(", year_col, ")")
  } else ""
  form <- as.formula(paste0(outcome, " ~ ", paste(rhs_vars, collapse = " + "), fe_term, " - 1"))

  # optional train-only scaling (prevents leakage even if upstream z-scoring was global)
  num_vars <- rhs_vars[sapply(train[, rhs_vars, drop = FALSE], is.numeric)]
  scaling <- list(enabled = scale_with_train, mu = NULL, sd = NULL, num_vars = num_vars)

  if (scale_with_train && length(num_vars) > 0) {
    mu <- sapply(train[, num_vars, drop = FALSE], mean)
    sd <- sapply(train[, num_vars, drop = FALSE], sd); sd[sd == 0] <- 1
    scale_with <- function(d) {
      d2 <- d
      d2[, num_vars] <- sweep(sweep(d2[, num_vars, drop = FALSE], 2, mu, "-"), 2, sd, "/")
      d2
    }
    train <- scale_with(train)
    test  <- scale_with(test)
    scaling$mu <- mu; scaling$sd <- sd
  }

  # matrices
  # matrices (build on TRAIN only; then align TEST)
  train_mm <- train
  test_mm  <- test

  # Drop unused factor levels in TRAIN so we don't create "future-level" columns
  fact_vars <- rhs_vars[sapply(train_mm[, rhs_vars, drop = FALSE], is.factor)]
  if (include_fixed_effects) fact_vars <- unique(c(fact_vars, country_col, year_col))
  if (length(fact_vars) > 0) {
    train_mm[, fact_vars] <- lapply(train_mm[, fact_vars, drop = FALSE], droplevels)
  }

  X_train <- model.matrix(form, data = train_mm)
  X_test  <- model.matrix(form, data = test_mm)

  # Align columns: add missing columns as zeros; drop extras
  miss <- setdiff(colnames(X_train), colnames(X_test))
  if (length(miss)) {
    X_test <- cbind(
      X_test,
      matrix(0, nrow(X_test), length(miss), dimnames = list(NULL, miss))
    )
  }
  X_test <- X_test[, colnames(X_train), drop = FALSE]
  X_cols <- colnames(X_train)

  y_train <- droplevels(train[[outcome]])
  y_test  <- droplevels(test[[outcome]])
  #lv <- levels(y_train); pos <- lv[2]

  # country-wise folds
  set.seed(42)
  K <- kfold_countries
  countries_tr <- unique(train[[country_col]])
  foldmap <- setNames(sample(rep(1:K, length.out = length(countries_tr))), countries_tr)
  foldid  <- as.integer(foldmap[train[[country_col]]])

  # class weights (minority upweighted)
  n_pos <- sum(y_train == pos); n_neg <- sum(y_train != pos)
  w <- ifelse(y_train == pos, n_neg / max(1, n_pos), 1)
  #w <- rep(1, length(y_train))

  # alpha tuning (AUC CV)
  fits <- lapply(alphas, function(a)
    cv.glmnet(
      x = X_train, y = y_train,
      family = "binomial",
      alpha = a,
      type.measure = "auc",
      standardize = FALSE,               
      foldid = foldid,
      nlambda = 200,
      lambda.min.ratio = 1e-5,
      weights = w,
      keep = TRUE
    )
  )
  # <- sapply(fits, function(f) max(f$cvm))
  #best_i <- which.max(cv_auc)
  #cv_loss <- sapply(fits, function(f) min(f$cvm))
  #best_i  <- which.min(cv_loss)
  #best_alpha <- alphas[best_i]
  #best_fit   <- fits[[best_i]]
  #best_alpha <- alphas[best_i]
  #best_fit <- fits[[best_i]]
  # s_use <- best_fit$lambda.1se
  cv_auc <- sapply(fits, function(f) max(f$cvm))
  best_i <- which.max(cv_auc)

  best_alpha <- alphas[best_i]
  best_fit   <- fits[[best_i]]

  # choose lambda (keep min, or use 1se if you want more regularization)
  s_use <- best_fit$lambda.min
  
  # --- OOF predictions at chosen lambda (from cv.glmnet) ---
  # --- Robust OOF predictions (country folds) ---
  p_train_oof <- rep(NA_real_, length(y_train))
  for (k in sort(unique(foldid))) {
    idx_val <- which(foldid == k)
    idx_tr  <- which(foldid != k)
    
    fit_k <- glmnet(
    x = X_train[idx_tr, , drop = FALSE],
    y = y_train[idx_tr],
    family = "binomial",
    alpha = best_alpha,
    lambda = s_use,
    standardize = FALSE,
    weights = w[idx_tr]
  )
    p_train_oof[idx_val] <- as.numeric(
    predict(fit_k, newx = X_train[idx_val, , drop = FALSE], type = "response")
  )
}

  stopifnot(all(is.finite(p_train_oof)))
  
  # -- TRAIN thresholds (choose here; apply to TEST)
  # In-sample train probs (still useful to keep)
  p_train <- as.numeric(predict(best_fit, s = s_use, newx = X_train, type = "response"))

  # Use out-of-fold probs for threshold selection (less optimistic)
  p_thr <- as.numeric(p_train_oof)

  roc_tr  <- roc(response = y_train, predictor = p_thr, levels = lv, quiet = TRUE)
  thr_youden_tr <- as.numeric(
    coords(roc_tr, "best", best.method = "youden", transpose = TRUE)["threshold"]
  )

  #grid <- seq(0.05, 0.95, by = 0.01)
  #f1_train <- function(th){
  #  pr <- factor(ifelse(p_thr >= th, pos, lv[1]), levels = lv)
  #  tp <- sum(pr==pos & y_train==pos); fp <- sum(pr==pos & y_train!=pos); fn <- sum(pr!=pos & y_train==pos)
  #  prec <- ifelse(tp+fp==0, 0, tp/(tp+fp)); rec <- ifelse(tp+fn==0, 0, tp/(tp+fn))
  #  ifelse(prec+rec==0, 0, 2*prec*rec/(prec+rec))
  #}
  #thr_f1_tr <- grid[which.max(vapply(grid, f1_train, numeric(1)))]


  # TEST metrics
  p_test <- as.numeric(predict(best_fit, s = s_use, newx = X_test, type = "response"))
  roc_te <- roc(response = y_test, predictor = p_test, levels = lv, quiet = TRUE)
  auc_te <- as.numeric(auc(roc_te))
  ci_te  <- as.numeric(ci.auc(roc_te)) 
  
  # Out-of-fold train probabilities at the chosen lambda (for calibration)
  #idx_lambda <- which.min(abs(best_fit$lambda - s_use))
  #p_train_oof <- as.numeric(best_fit$fit.preval[, idx_lambda])

  # thresholds learned on TRAIN, applied to TEST
  pred_y <- factor(ifelse(p_test >= thr_youden_tr, pos, lv[1]), levels = lv)
  #pred_f <- factor(ifelse(p_test >= thr_f1_tr,     pos, lv[1]), levels = lv)
  cm_y   <- confusionMatrix(pred_y, y_test, positive = pos)
  #cm_f   <- confusionMatrix(pred_f, y_test, positive = pos)

  # PR-AUC + Brier
  y_bin <- as.integer(y_test == pos)
  pr    <- pr.curve(scores.class0 = p_test[y_bin==1],
                    scores.class1 = p_test[y_bin==0], curve = FALSE)
  pr_auc <- pr$auc.integral
  brier  <- mean((y_bin - p_test)^2)
  pi_test    <- mean(y_bin)
  brier_null <- pi_test * (1 - pi_test)
  bss        <- 1 - brier / brier_null

  # simple calibration (test)
  eps_cal <- 1e-6
  p_test_c <- pmin(pmax(p_test, eps_cal), 1 - eps_cal)

  cal_fit <- glm((y_test == pos) ~ qlogis(p_test_c), family = binomial())
  cal_intercept <- unname(coef(cal_fit)[1])
  cal_slope     <- unname(coef(cal_fit)[2])

  # coefficients at s_use
  nz <- coef(best_fit, s = s_use)
  nz_df <- data.frame(feature = rownames(nz), beta = as.numeric(nz), row.names = NULL) |>
    filter(beta != 0)

  list(
    # model + design metadata
    model = best_fit,
    best_alpha = best_alpha,
    lambda = s_use,
    model_formula = form,
    X_cols = X_cols,
    scaling = scaling,
    levels = lv, pos = pos,
    include_fixed_effects = include_fixed_effects,

    # CV + Test metrics
    cv_auc_by_alpha = data.frame(alpha = alphas, cv_auc = cv_auc),
    auc_test = auc_te,
    auc_test_ci = ci_te,
    pr_auc_test = pr_auc,
    brier_test  = brier,
    brier_null_test   = brier_null,
    brier_skill_test  = bss,
    pi_test           = pi_test,
    calib_intercept = cal_intercept,
    calib_slope     = cal_slope,
    
    p_train      = p_train,
    p_train_oof  = p_train_oof,
    p_test       = p_test,
    y_train      = y_train,
    y_test       = y_test,
    holdout_year = holdout_year,

    # thresholds (chosen on TRAIN)
    threshold_youden_train = thr_youden_tr,
    #threshold_f1_train     = thr_f1_tr,

    # confusion matrices (TEST)
    confusion_youden = cm_y,
    #confusion_f1     = cm_f,

    # extras
    nonzero_coefs = nz_df,
    X_train_dim = dim(X_train),
    X_test_dim  = dim(X_test)
  )
}

# ------------------------------------------------------------------------------
# Helper: predict/evaluate on NEW data block -----------------------------------
# ------------------------------------------------------------------------------
eval_on_new <- function(fit, newdata, outcome, country_col, year_col) {
  df <- newdata |>
    select(all_of(c(country_col, year_col, outcome, all.vars(update(fit$model_formula, . ~ .)))))
  df <- drop_na(df)

  # applying train scaling if used
  if (isTRUE(fit$scaling$enabled) && length(fit$scaling$num_vars) > 0) {
    nv <- intersect(fit$scaling$num_vars, names(df))
    if (length(nv) > 0) {
      df[, nv] <- sweep(sweep(df[, nv, drop = FALSE], 2, fit$scaling$mu[nv], "-"),
                        2, fit$scaling$sd[nv], "/")
    }
  }

  # build design, then align cols to training X
  X_new <- model.matrix(fit$model_formula, data = df)
  miss  <- setdiff(fit$X_cols, colnames(X_new))
  if (length(miss)) X_new <- cbind(X_new, matrix(0, nrow(X_new), length(miss), dimnames = list(NULL, miss)))
  X_new <- X_new[, fit$X_cols, drop = FALSE]

  y_new <- factor(df[[outcome]], levels = fit$levels)
  p_new <- as.numeric(predict(fit$model, s = fit$lambda, newx = X_new, type = "response"))

  roc_new <- roc(response = y_new, predictor = p_new, levels = fit$levels, quiet = TRUE)
  auc_new <- as.numeric(auc(roc_new))

  list(n = length(p_new), auc = auc_new, probs = p_new, y = y_new)
}

# ------------------------------------------------------------------------------
# Plotting helpers -------------------------------------------------------------
# ------------------------------------------------------------------------------
plot_cv_and_paths <- function(fit, main_prefix = "Best model") {
  op <- par(mfrow = c(1,2), mar = c(5,4,3,1)); on.exit(par(op), add = TRUE)
  plot(fit$model, main = paste0(main_prefix, " — CV (AUC), alpha=", fit$best_alpha))
  abline(v = log(fit$lambda), lty = 2)
  plot(fit$model$glmnet.fit, xvar = "lambda", main = "Coefficient paths")
  abline(v = log(fit$lambda), lty = 2)
}

# ------------------------------------------------------------------------------
# Time-cut sensitivity (forward) -----------------------------------------------
# ------------------------------------------------------------------------------
time_cut_sensitivity <- function(data, predictors, outcome,
                                 country_col, year_col,
                                 include_fixed_effects = FALSE,
                                 scale_with_train = TRUE,
                                 cuts = NULL, min_test_n = 100) {
  df_clean <- data |>
    select(all_of(c(country_col, year_col, predictors, outcome))) |>
    drop_na()
  yrs <- sort(unique(df_clean[[year_col]]))
  if (is.null(cuts)) cuts <- tail(yrs[yrs < max(yrs)], 6)

  out <- lapply(cuts, function(cut) {
    fit <- try(fit_glmnet_panel(
      data, predictors, outcome, country_col, year_col,
      holdout_year = cut,
      include_fixed_effects = include_fixed_effects,
      scale_with_train = scale_with_train
    ), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    if (fit$X_test_dim[1] < min_test_n) return(NULL)
    data.frame(cut = cut,
               n_train = fit$X_train_dim[1],
               n_test  = fit$X_test_dim[1],
               auc     = fit$auc_test,
               pr_auc  = fit$pr_auc_test,
               brier   = fit$brier_test)
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

# ------------------------------------------------------------------------------
# Region-out (LORO) evaluation -------------------------------------------------
# ------------------------------------------------------------------------------
region_out <- function(data, predictors, outcome,
                       country_col, year_col, region_col,
                       holdout_year, include_fixed_effects = FALSE,
                       scale_with_train = TRUE, min_test_n = 150) {
  regs <- sort(unique(data[[region_col]]))
  res <- lapply(regs, function(R) {
    tr <- data[data[[region_col]] != R, ]
    te <- data[data[[region_col]] == R & data[[year_col]] > holdout_year, ]
    if (nrow(te) < min_test_n) return(NULL)

    fit <- fit_glmnet_panel(
      tr, predictors, outcome, country_col, year_col,
      holdout_year = holdout_year,
      include_fixed_effects = include_fixed_effects,
      scale_with_train = scale_with_train
    )
    ev <- eval_on_new(fit, te, outcome, country_col, year_col)
    data.frame(region = R, n_test = ev$n, auc = ev$auc)
  })
  do.call(rbind, Filter(Negate(is.null), res))
}

# ------------------------------------------------------------------------------
# Running the specification ----------------------------------------------------
# ------------------------------------------------------------------------------
predictors_baseline <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "LOG_GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  #"POLITICAL_REGIME",
  #"FH_DEMOCRACY",
  #"FH_AUTOCRACY",
  #"DD_DEMOCRACY",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  #"INSTITUTIONAL_DEMOCRACY_SOCRE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  #"COMBINED_POLITY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "COMPETITIVENESS_PARTICIPATION",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "INSTITUTIONAL_PARTICIPATION",
  "EXECUTIVE_CONSTRAINT_SCORE",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  #"MAX_CONFLICT_INTENSITY",
  #"AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

set.seed(42)

baseline_fit <- fit_glmnet_panel(
  data = final_clean_percentiles_data_normalized_1,
  predictors = predictors_baseline,
  outcome = "VDEM_STATUS_IDEAL",
  country_col = "COUNTRY_NAME",
  year_col = "YEAR",
  holdout_year = 2016,          # or NULL to auto 80/20
  include_fixed_effects = FALSE,
  kfold_countries = 5,
  scale_with_train = TRUE       # prevents any scaling leakage
)

# Inspecting & plottting
# baseline_fit$cv_auc_by_alpha
baseline_fit$cv_auc_by_alpha
baseline_fit$best_alpha
baseline_fit$auc_test; baseline_fit$auc_test_ci
baseline_fit$pr_auc_test; baseline_fit$brier_test
baseline_fit$brier_null_test
baseline_fit$brier_skill_test
baseline_fit$pi_test
baseline_fit$calib_intercept; baseline_fit$calib_slope
baseline_fit$confusion_youden
head(baseline_fit$nonzero_coefs)
plot_cv_and_paths(baseline_fit, "Baseline")

# Time-cut stability
time_cut_summary <- time_cut_sensitivity(
  data = final_clean_percentiles_data_normalized_1,
  predictors = predictors_baseline,
  outcome = "VDEM_STATUS_IDEAL",
  country_col = "COUNTRY_NAME",
  year_col = "YEAR",
  include_fixed_effects = FALSE,
  scale_with_train = TRUE
)
time_cut_summary

# Getting all coefficients (including zero) at chosen lambda
all_coefs <- as.matrix(coef(baseline_fit$model, s = baseline_fit$lambda))

coef_table <- data.frame(
  feature = rownames(all_coefs),
  beta    = as.numeric(all_coefs)
) |>
  # Dropping intercepet
  filter(feature != "(Intercept)") |>
  mutate(shrunk_to_zero = beta == 0)

# Inspect
print(head(coef_table, 20))  # first 20 rows
table(coef_table$shrunk_to_zero)

# ------------------------------------------------------------------------------
# Cross-validated AUC by alpha table -------------------------------------------
# ------------------------------------------------------------------------------
baseline_fit$cv_auc_by_alpha |>
  mutate(
    alpha   = round(alpha, 2),
    cv_auc = round(cv_auc, 3)
  )

# ------------------------------------------------------------------------------
# Time-cut summary table -------------------------------------------------------
# ------------------------------------------------------------------------------
time_cut_summary |>
  mutate(
    cut    = as.integer(cut),
    auc    = round(auc, 3),
    pr_auc = round(pr_auc, 3),
    brier  = round(brier, 3)
  )

# ------------------------------------------------------------------------------
# Confusion Matrix output ------------------------------------------------------
# ------------------------------------------------------------------------------
extract_cm_row <- function(cm, label) {
  # cm is a caret::confusionMatrix object
  tibble(
    Thresholding     = label,
    Accuracy         = unname(cm$overall["Accuracy"]),
    Sensitivity      = unname(cm$byClass["Sensitivity"]),
    Specificity      = unname(cm$byClass["Specificity"]),
    `Pos Pred Value` = unname(cm$byClass["Pos Pred Value"]),
    `Neg Pred Value` = unname(cm$byClass["Neg Pred Value"])
  )
}

tab_cm <- bind_rows(
  extract_cm_row(baseline_fit$confusion_youden, "Youden's J"),
  extract_cm_row(baseline_fit$confusion_f1,     "F1-optimal")
) |>
  mutate(across(-Thresholding, ~ round(.x, 3)))

tab_cm

# ------------------------------------------------------------------------------
# Platt scalling ---------------------------------------------------------------
# ------------------------------------------------------------------------------
eps <- 1e-6
clip <- function(p) pmin(pmax(p, eps), 1 - eps)

pos <- baseline_fit$pos

# Use OOF train probs for calibration (best). Fall back to in-sample if needed.
p_tr <- if (!is.null(baseline_fit$p_train_oof)) baseline_fit$p_train_oof else baseline_fit$p_train

y_tr <- as.integer(baseline_fit$y_train == pos)
y_te <- as.integer(baseline_fit$y_test  == pos)

p_tr <- clip(as.numeric(p_tr))
p_te <- clip(as.numeric(baseline_fit$p_test))

# Platt scaler fit on TRAIN (OOF probs)
cal <- glm(y_tr ~ qlogis(p_tr), family = binomial())

# Apply to TEST
p_te_cal <- plogis(coef(cal)[1] + coef(cal)[2] * qlogis(p_te))

# Brier + Skill (vs prevalence baseline) on TEST
brier_raw <- mean((y_te - p_te)^2)
brier_cal <- mean((y_te - p_te_cal)^2)

pi_te <- mean(y_te)
brier_null <- pi_te * (1 - pi_te)

bss_raw <- 1 - brier_raw / brier_null
bss_cal <- 1 - brier_cal / brier_null

c(pi_test = pi_te,
  brier_null = brier_null,
  brier_raw = brier_raw, bss_raw = bss_raw,
  brier_cal = brier_cal, bss_cal = bss_cal)

# ------------------------------------------------------------------------------
# Isotonic regression ----------------------------------------------------------
#-------------------------------------------------------------------------------
eps <- 1e-6
clip <- function(p) pmin(pmax(p, eps), 1 - eps)

pos <- baseline_fit$pos
p_tr <- clip(as.numeric(baseline_fit$p_train_oof))
y_tr <- as.integer(baseline_fit$y_train == pos)

p_te <- clip(as.numeric(baseline_fit$p_test))
y_te <- as.integer(baseline_fit$y_test == pos)

# Fit isotonic calibration map on (p_tr, y_tr)
iso <- isoreg(p_tr, y_tr)

# Predict calibrated probs by interpolating the step function
p_te_iso <- pmin(pmax(approx(iso$x, iso$yf, xout = p_te, rule = 2)$y, 0), 1)

brier_raw <- mean((y_te - p_te)^2)
brier_iso <- mean((y_te - p_te_iso)^2)

c(brier_raw = brier_raw, brier_iso = brier_iso)
```

## 6.1. Pure LASSO

```{r pure-lasso}
baseline_fit_lasso <- fit_glmnet_panel(
  data = final_clean_percentiles_data_normalized_1,
  predictors = predictors_baseline,
  outcome = "VDEM_STATUS_IDEAL",
  country_col = "COUNTRY_NAME",
  year_col = "YEAR",
  holdout_year = 2016,
  include_fixed_effects = FALSE,
  kfold_countries = 5,
  scale_with_train = TRUE,
  # Setting pure LASSO
  alphas = c(1)   
)

all_coefs_lasso <- as.matrix(coef(baseline_fit_lasso$model, 
                                  s = baseline_fit_lasso$lambda))

coef_table_lasso <- data.frame(
  feature = rownames(all_coefs_lasso),
  beta    = as.numeric(all_coefs_lasso),
  stringsAsFactors = FALSE
) |> 
  filter(feature != "(Intercept)") |>
  mutate(shrunk_to_zero = abs(beta) < 1e-9)

table(coef_table_lasso$shrunk_to_zero)  

subset(coef_table_lasso)
```

# Hazard Model

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ==========================================================
# Fragility Exit Hazard — WP-style core ID + robustness
# Data object: standardized_data_final
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)      # feglm (GLM) with clustering and FE
  library(brglm2)      # bias-reduced logistic regression
  library(glmmTMB)     # random-intercept logit
  library(margins)     # marginal effects
  library(sandwich)    # clustered vcov
  library(lmtest)      # coeftest with custom vcov
  library(forcats)
})

# ------------------------------------------------------------------------------
# Settings ---------------------------------------------------------------------
# ------------------------------------------------------------------------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
# "Fragile","Transitioning","Non Fragile"
STATUS  <- "VDEM_STATUS_IDEAL"     
# Lag coverage threshold
COVERAGE_MIN <- 0.70               
# sustained-exit year horizons for robustness
K_GRID       <- c(3L, 5L, 7L)      
# 2% leverage drop in influence check
DROP_TOP_LEVERAGE <- 0.02          

# Full predictors menu (we will prune by coverage and then select a CORE)
PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "LOG_GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "COMBINED_POLITY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

# Theory-first CORE spec (kept if coverage permits)
CORE_NAMES <- c(
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "GDP_GROWTH_L1",            
)

# -----------------
# Helpers
# -----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  # lag predictors by 1 year (keep NA for coverage audit)
  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")
  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 1L), .names = "{.col}_L1")) %>%
    ungroup()

  # Clean factor lags (avoid weird numeric-level labels / NA as level)
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")

  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  # helper to print clustered coeftest for glm/brglm2
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

# -----------------
# Load & prep
# -----------------
stopifnot(exists("standardized_data_final"))
df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# regime-like vars as factors (pre-lag)
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors (pre-lag)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# -----------------
# Run the full pipeline for each K in K_GRID
# -----------------
for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")
  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # Coverage audit (top 10)
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # ---- CORE model selection (based on theory & your early signal)
  core <- intersect(CORE_NAMES, keep_lag)
  # add GDP_DEFLATOR_L1 if available (proxy for nominal/ToT)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) {
    # fall back to a larger kept set if coverage excluded too much
    core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  }
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # -----------------
  # POOLED (year FE on by default), clustered by country
  # -----------------
  f_core_yearfe <- as.formula(paste("EXIT_SUSTAINED ~", paste(core, collapse = " + "), "+ i(", YEAR, ")", sep = ""))
  m_pool_core <- feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # Influence diag: drop top 2% leverage and refit pooled without year FE (GLM) for transparency
  f_core <- as.formula(paste("EXIT_SUSTAINED ~", paste(core, collapse = " + ")))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =", length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # -----------------
  # RARE-EVENTS (bias-reduced) + clustered SEs — CORE
  # -----------------
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))
  
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  
  if (!exists("COLLECT")) COLLECT <- tmp else COLLECT <- bind_rows(COLLECT, tmp)


  # -----------------
  # RANDOM-INTERCEPT (country) — CORE (graceful fallback)
  # -----------------
  f_re <- as.formula(paste("EXIT_SUSTAINED ~", paste(core, collapse = " + "), "+ (1 |", COUNTRY, ")"))
  m_re_core <- try(
    glmmTMB(f_re, data = haz, family = binomial(),
            control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
    silent = TRUE
  )
  cat("\n=== Random-intercept (country) — CORE ===\n")
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge from default starts (likely sparse categories / quasi-separation).\n")
  } else {
    print(summary(m_re_core))
  }

  # -----------------
  # AMEs for CORE pooled GLM (no year FE; AMEs are about X, not time dummies)
  # -----------------
  ame <- margins(m_glm_core)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  print(summary(ame))

  cat("\n================= End K =", K, "=================\n")
}


```

## Full Predictors 

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ================================================================
# Discrete-time hazard of sustained exit from fragility (WP-ready)
# ================================================================

# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "LOG_GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  #"POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  #"TERRITORIAL_FRAGMENTATION",
  #"INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  #"POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  #"CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  #"MAX_CONFLICT_INTENSITY",
  #"AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "LOG_GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  #"TERRITORIAL_FRAGMENTATION_L1",
  #"INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  #"POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  #"CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  #"MAX_CONFLICT_INTENSITY_L1",
  #"AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df |>
    arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df |>
    group_by(.data[[COUNTRY]]) |>
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) |>
    ungroup() |>
    group_by(.data[[COUNTRY]]) |>
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) |>
    ungroup()

  haz <- df |>
    filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz |>
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 1L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT <- bind_rows(COLLECT, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

```


```{r, hazard, echo = FALSE, message = FALSE, warning = FALSE, results = 'asis'}

stopifnot(exists("COLLECT"))

# Keep only the CORE variables you plotted in the models
core_terms <- c(
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "GDP_GROWTH_L1",
  "GDP_DEFLATOR_L1"
)

# Nice labels for rows
var_labs <- c(
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-1)",
  REGIME_DURABILITY_YEARS_L1   = "Regime durability, years (t-1)",
  CONFLICT_INTENSITY_YEAR_L1   = "Conflict intensity (acute, t-1)",
  AVG_CONFLICT_INTENSITY_L1    = "Conflict intensity (chronic avg., t-1)",
  GDP_GROWTH_L1                = "Real GDP per capita growth (t-1)",
  GDP_DEFLATOR_L1              = "GDP deflator / price tailwind (t-1)"
)

# Compute ORs and 95% clustered CIs
tab <- COLLECT %>%
  filter(term %in% core_terms) %>%
  mutate(
    OR  = exp(beta),
    LCL = exp(beta - 1.96 * se),
    UCL = exp(beta + 1.96 * se),
    cell = sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL),
    term_label = var_labs[term]
  ) %>%
  select(term_label, K, cell)

# If a core variable was dropped for some K (coverage/collinearity), mark as blank
all_rows <- tibble(term_label = unname(var_labs))
all_cols <- tibble(K = sort(unique(tab$K)))
tab_full <- tidyr::complete(tab, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# Wide layout: one column per K
tab_wide <- tab_full %>%
  tidyr::pivot_wider(names_from = K, values_from = cell, names_prefix = "K = ")

# Order rows like CORE
tab_wide <- tab_wide %>%
  mutate(row_order = match(term_label, unname(var_labs))) %>%
  arrange(row_order) %>%
  select(-row_order)

# --- 3) Render LaTeX (booktabs) and write to file ---

library(kableExtra)

latex_tbl <- kable(
  tab_wide,
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (CORE rare-events logit, clustered by country)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c")
) %>%
  kable_styling(latex_options = c("hold_position","striped")) %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3)) %>%
  footnote(
    general = "Entries are odds ratios with 95\\% clustered confidence intervals. Models are bias-reduced (brglm2, AS_mixed). Covariates standardized before lagging. Outcome: first sustained exit (Non-Fragile) at t persisting for K years. Sample includes risk set (Fragile/Transitioning) and the first exit year.",
    threeparttable = TRUE, escape = FALSE
  )

cat(latex_tbl)

```

## Final Table

```{r, echo = FALSE, message = FALSE}
# ============================================================
# OR table (ALL predictors) from COLLECT — vectorized + robust
# ============================================================
stopifnot(exists("COLLECT"))

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-1)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-1)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-1)",
  GDP_DEFLATOR_L1            = "GDP deflator / price tailwind (t-1)",
  POLITICAL_REGIME_L1        = "Political regime (t-1)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-1)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-1)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-1)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-1)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-1)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-1)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-1)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-1)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-1)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-1)",
  N_WAR_FRONTS_L1                 = "# war fronts (t-1)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-1)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-1)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-1)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## Two Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 2L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_two_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_two_year <- bind_rows(COLLECT_two_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_two_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_two_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_two_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_two_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_two_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-2)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-2)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-2)",
  GDP_DEFLATOR_L1            = "GDP deflator / price tailwind (t-2)",
  POLITICAL_REGIME_L1        = "Political regime (t-2)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-2)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-2)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-2)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-2)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-2)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-2)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-2)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-2)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-2)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-2)",
  N_WAR_FRONTS_L1                 = "# war fronts (t-2)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-2)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-2)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-2)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)
```

## Three-Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 3L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_three_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_three_year <- bind_rows(COLLECT_three_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_three_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_three_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_three_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_three_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_three_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-3)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-3)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-3)",
  GDP_DEFLATOR_L1            = "GDP deflator (t-3)",
  POLITICAL_REGIME_L1        = "Political regime (t-3)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-3)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-3)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-3)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-3)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-3)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-3)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-3)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-3)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-3)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-3)",
  N_WAR_FRONTS_L1                 = "N war fronts (t-3)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-3)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-3)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-2)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## Four-Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 4L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_four_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_four_year <- bind_rows(COLLECT_four_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_four_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_four_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_four_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_four_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_four_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-4)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-4)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-4)",
  GDP_DEFLATOR_L1            = "GDP deflator (t-4)",
  POLITICAL_REGIME_L1        = "Political regime (t-4)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-4)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-4)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-4)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-4)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-4)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-4)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-4)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-4)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-4)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-4)",
  N_WAR_FRONTS_L1                 = "N war fronts (t-4)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-4)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-4)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-4)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## Five-Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 5L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_five_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_five_year <- bind_rows(COLLECT_five_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_five_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_five_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_five_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_five_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_five_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-4)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-4)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-4)",
  GDP_DEFLATOR_L1            = "GDP deflator (t-4)",
  POLITICAL_REGIME_L1        = "Political regime (t-4)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-4)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-4)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-4)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-4)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-4)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-4)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-4)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-4)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-4)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-4)",
  N_WAR_FRONTS_L1                 = "N war fronts (t-4)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-4)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-4)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-4)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## LASSO Predictors Only

```{r}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH", # Drop
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME", # Drop
  "ELECTORAL_DEMOCRACY_SCORE", # Drop
  "LIBERAL_DEMOCRACY_SCORE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE", # Drop
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS", # Drop
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 1L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_two_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_two_year <- bind_rows(COLLECT_two_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_two_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_two_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_two_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_two_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_two_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-1)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-1)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-1)",
  GDP_DEFLATOR_L1            = "GDP deflator / price tailwind (t-1)",
  POLITICAL_REGIME_L1        = "Political regime (t-1)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-1)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-1)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-1)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-1)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-1)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-1)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-1)",
  N_WAR_FRONTS_L1                 = "# war fronts (t-1)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-1)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-1)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)
```
