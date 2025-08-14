---
title: "FRAGILITY EXIT - Final VDEM Models and Findings"
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

# Cleaning my environment
rm(list = ls())

# Managing memory
gc()

# Load packages ---------------------------------------------------------------
library(here)
library(coefplot)
library(pROC)
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
library(caret)
library(stargazer)
library(VGAM)
library(factoextra)
library(cluster)
library(ggplot2)
library(vip)
library(themis)
library(janitor)
library(WDI)

# Setting seed for replication -------------------------------------------------
set.seed(247)

# Load data, functions, and scripts --------------------------------------------
load("final_new_clean_output_data.RDS")
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

# Performing the join ----------------------------------------------------------
final_clean_percentiles_data <- left_join(final_clean_percentiles_data,
                                          gdp_per_capita_data,
                                          by = c("ISO_CODE_3", 
                                                 "YEAR"),
                                          keep = FALSE)

# Final cleanning --------------------------------------------------------------
final_clean_percentiles_data <- final_clean_percentiles_data |>
  # Dropping duplicate columns 
  rename(COUNTRY_NAME = COUNTRY_NAME.x) |>
  select(-COUNTRY_NAME.y) |>
  # Allocating GDP per capita along with economic indicators 
  relocate(82, .after = 3)

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
    across(35:55, as.character)  # If needed for export or plotting
  )

# Computing means and SDs for numeric columns ----------------------------------
means_data <- sapply(final_data, function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else NA)
sd_data    <- sapply(final_data, function(x) if (is.numeric(x)) sd(x, na.rm = TRUE) else NA)

# Initializing standardized data
standardized_data <- final_data

# Defining dummies to EXCLUDE from standardization
dummy_vars <- c(
  "VDEM_FRAGILE_BASELINE", 
  "VDEM_NON_FRAGILE_BASELINE",
  "VDEM_FRAGILE_IDEAL", 
  "VDEM_NON_FRAGILE_IDEAL",
  "VDEM_FRAGILE_ENDLINE", 
  "VDEM_NON_FRAGILE_ENDLINE"
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

\newpage

# 1. Executive Summary

This brief summary estimates:

* Country-year linear and logistic coefficients of 20 numeric predictors on a dichotomous outcome: whether a country is `Fragile` or `Non-Fragile` in a given year. The strategy is repeated with two-year lags of the same 20 numeric predictors. 

* The performance of the linear model (which returns the marginal effect of of one unit increase in the predictor over the target variable) and of the logistic model (which estimates the log-likelihood that the dichotomous outcome is true for a unit increase in the predictor). Predictors are normalized, so coefficients are directly comparable across the two models: the coefficients capture the change in the target variable for one standard deviation (positive or negative) in the predictor variable. The strategy is repeated with two-year lags of the same 20 numeric predictors.

* The LASSO regression parameters, tackling overfitting and multicollinaerity in identifying less relevant predictors. 

* The confusion matrix, which evaluates the accuracy of trained data in correctly identifying the class a country falls in for a given year (whether `Fragile` or `Non-Fragile`)

\newpage

# 2. Descriptive Statistics 

```{r, echo = FALSE, message = FALSE, Wwarning = FALSE}

descriptive_data <- final_clean_percentiles_data |>
  # Selecting only relevant predictors
  select(GDP_PER_CAPITA,
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

descriptive_data <- descriptive_data |>
  mutate(GDP_PER_CAPITA_LOG = log(GDP_PER_CAPITA)) |>
  relocate(23, .after = 1)
```

```{r, descriptive-stats, echo = FALSE, message = FALSE, warning = FALSE, results = 'asis'}

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

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# Time trends ------------------------------------------------------------------

# Count number of countries per category per year
trend_data <- standardized_data_final |>
  group_by(YEAR, 
           VDEM_STATUS_IDEAL) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(YEAR = as.numeric(YEAR)) |>
  arrange(VDEM_STATUS_IDEAL,
          YEAR)

# Plotting the time trends
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

# Plot
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
  facet_wrap(~ gdp_pc_quartile, ncol = 2) +  # rectangle headers you liked
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "top"
  )
```

\newpage 

# 3. Success Episodes 

## 3.1. ODA: Pre vs. Post-Escape Pooled Estimates

The scatterplot below shows, for each case of successful breakthrough, the relationship between the average ODA per capita in the pre-escape period (x-axis, log1p-transformed) and the change in ODA per capita from pre- to post-escape (y-axis). Each dot represents one country. The black line is an OLS fit with a 95% confidence interval, summarizing the average association across countries.

```{r, echo = FALSE, message = FALSE, warning = FALSE}

# Setting the critical junctures
standardized_data_final <- standardized_data_final |>
  mutate(VDEM_ESCAPE_YEAR = case_when(
    (COUNTRY_NAME == "Brazil" & YEAR == 1988) ~ 1,
    (COUNTRY_NAME == "Gabon" & YEAR == 1990) ~ 1,
    (COUNTRY_NAME == "Gambia, The" & YEAR == 2017) ~ 1,
    (COUNTRY_NAME == "Kenya" & YEAR == 2010) ~ 1,
    (COUNTRY_NAME == "Mauritania" & YEAR == 2006) ~ 1,
    (COUNTRY_NAME == "Sierra Leone" & YEAR == 2018) ~ 1,
    (COUNTRY_NAME == "Tunisia" & YEAR == 2010) ~ 1,
    (COUNTRY_NAME == "Uganda" & YEAR == 1987) ~ 1,
    (COUNTRY_NAME == "Iran, Islamic Rep." & YEAR == 1997) ~ 1,
    (COUNTRY_NAME == "Nepal" & YEAR == 2007) ~ 1,
    (COUNTRY_NAME == "Philippines" & YEAR == 1987) ~ 1,
    (COUNTRY_NAME == "Sri Lanka" & YEAR == 2009) ~ 1,
    (COUNTRY_NAME == "Thailand" & YEAR == 1991) ~ 1,
    (COUNTRY_NAME == "Bosnia and Herzegovina" & YEAR == 1995) ~ 1,
    (COUNTRY_NAME == "Serbia" & YEAR == 2001) ~ 1,
    (COUNTRY_NAME == "Georgia" & YEAR == 2004) ~ 1,
    (COUNTRY_NAME == "Kyrgyz Republic" & YEAR == 2011) ~ 1,
    (COUNTRY_NAME == "Armenia" & YEAR == 2018) ~ 1,
    (COUNTRY_NAME == "Argentina" & YEAR == 1983) ~ 1,
    (COUNTRY_NAME == "Colombia" & YEAR == 2002) ~ 1,
    (COUNTRY_NAME == "Ecuador" & YEAR == 1976) ~ 1,
    (COUNTRY_NAME == "Dominican Republic" & YEAR == 1997) ~ 1,
    (COUNTRY_NAME == "El Salvador" & YEAR == 2009) ~ 1,
    (COUNTRY_NAME == "Mexico" & YEAR == 1994) ~ 1,
    (COUNTRY_NAME == "Nicaragua" & YEAR == 1990) ~ 1,
    (COUNTRY_NAME == "Panama" & YEAR == 1990) ~ 1,
    (COUNTRY_NAME == "Peru" & YEAR == 2001) ~ 1,
    TRUE ~ 0
  )) |>
  mutate(VDEM_ESCAPE_YEAR = as.numeric(VDEM_ESCAPE_YEAR))

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
```

Most countries experienced a positive change in ODA per capita after their escape year (dots above zero), although the magnitude varies widely. The regression line has a slight upward slope, suggesting that countries with higher baseline ODA per capita tended, on average, to see somewhat larger increases afterward. We recognize the broad confidence interval. 

## 3.2. ODA: Pre vs. Post-Escape Country-Estimates

```{r, echo = FALSE, message = FALSE, warning = FALSE}

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
  tidyr::pivot_wider(names_from = period, 
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
```

\newpage

# 4. A Classification Problem 

To validate the internal consistency and empirical grounding of our fragility classification (`Fragile`, `Transitioning` or `Non-Fragile`), we implemented a supervised machine learning exercise using two distinct algorithms: multinomial logistic regression (logit) and random forest. We trained both models on a set of 20 theoretically motivated structural predictors, including regime type, democracy indicators, institutional durability, economic performance, and conflict intensity (i.e., all predictors we employ in our empirical specifications). The outcome variable was our three-tier classification (\texttt{VDEM\_STATUS\_IDEAL}), derived independently of the predictors via percentile cutoffs on a normalized V-Dem latent factor. 

We label our country cases according to the subsequent taxonomy:

| Classification Label     | Fragile Threshold            | Transitioning Range              | Non-Fragile Threshold           |
|--------------------------|------------------------------|----------------------------------|---------------------------------|
| `RELAXED`                | Below 20th percentile        | 20th to 25th percentile          | Above 25th percentile           |
| `IDEAL`                  | Below 30th percentile        | 30th to 35th percentile          | Above 35th percentile           |
| `STRICT`                 | Below 20th percentile        | 20th to 40th percentile          | Above 40th percentile           |


The two supervised learning algorithms serve complementary purposes to validate our classification:

- **Multinomial Logistic Regression (Logit)** estimates the probability that a given observation belongs to one of three mutually exclusive classes (`Fragile`, `Transitioning` or `Non-Fragile`), based on a linear combination of structural predictors. It is fully interpretable, assumes additive relationships, and provides insight into the marginal effects of each variable. Because it is transparent and theory-aligned, it serves as an ideal benchmark for testing the recoverability of our classification from exogenous structural features.

- **Random Forest** is a non-parametric ensemble method that builds hundreds of decision trees on random subsets of predictors and observations. It excels at modeling complex, nonlinear interactions and is robust to multicollinearity and overfitting. In our context, it offers a high-performing, assumption-light alternative that tests whether the classification can still be recovered without relying on linearity or additivity assumptions.

\newpage 

## 4.1. Ideal Classification

Model performance was evaluated using repeated 10-fold cross-validation (3 repeats, 30 folds total). The multinomial logit model achieved an average classification accuracy of 91.3\% and a multiclass log-loss of 0.48, with low standard errors (0.0038 and 0.0242, respectively). The random forest model, while slightly less accurate at 89.5\%, yielded a lower log-loss of 0.31, indicating superior calibration of predicted probabilities. Confusion matrices reveal that most misclassifications occurred along conceptual boundaries—e.g., Fragile vs.\ Transitioning, or Transitioning vs.\ Non-Fragile — rather than across extreme categories (i.e., Fragile vs. Non Fragile). Both models exhibited high stability and low variance across cross-validation folds, suggesting that the classification scheme aligns closely with observable patterns in structural data.

These results support the validity of our classification strategy: not only is it recoverable from exogenous institutional, economic, and conflict indicators, but its internal logic also reflects meaningful distinctions between country-years across the fragility spectrum. In sum, the exercise confirms that our tiered labels are empirically grounded and not artifacts of arbitrary thresholds.

```{r, echo = FALSE, message = FALSE, warning = FALSE}

# Step 1: Defining our predictors predictors
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

# Building our dataset
df <- standardized_data_final |>
  select(all_of(predictors_20), 
         VDEM_STATUS_IDEAL) |>
  filter(!is.na(VDEM_STATUS_IDEAL)) |>
  drop_na() |>
  mutate(VDEM_STATUS_IDEAL = as.factor(VDEM_STATUS_IDEAL))

# Train/test split
data_split <- initial_split(df, strata = VDEM_STATUS_IDEAL)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Recipe
rec <- recipe(VDEM_STATUS_IDEAL ~ ., data = train_data) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# Model specifications ----------------------------------------------------------

# Logit Model 
logit_spec <- multinom_reg(mode = "classification") |> 
  set_engine("nnet")

# Random Forest Model 
rf_spec <- rand_forest(mode = "classification", 
                       trees = 500) |>
  set_engine("ranger", 
             importance = "permutation")

# Step 6: Workflows
logit_wf <- workflow() |> add_model(logit_spec) |> add_recipe(rec)
rf_wf    <- workflow() |> add_model(rf_spec)    |> add_recipe(rec)

# Step 7: CV setup
folds <- vfold_cv(train_data, v = 10, repeats = 3, strata = VDEM_STATUS_IDEAL)

# Step 8: CV model evaluation
logit_res <- fit_resamples(logit_wf, folds, metrics = metric_set(accuracy, mn_log_loss))
rf_res    <- fit_resamples(rf_wf, folds, metrics = metric_set(accuracy, mn_log_loss))

# Step 9: Final model fits
logit_fit <- fit(logit_wf, data = train_data)
rf_fit    <- fit(rf_wf, data = train_data)

# Step 10: Predict
logit_preds <- predict(logit_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

rf_preds <- predict(rf_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

# Step 11: Compute confusion matrices + accuracy
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds, truth = truth, estimate = .pred_class)

accuracy_logit <- accuracy(logit_preds, truth = truth, estimate = .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds, truth = truth, estimate = .pred_class)$.estimate

# Step 12: Plot confusion matrix with custom colors + accuracy in title
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(
    title = paste0("Logit Confusion Matrix (Accuracy = ", round(accuracy_logit * 100, 1), "%)"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(
    title = paste0("RE Confusion Matrix (Accuracy = ", round(accuracy_rf * 100, 1), "%)"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

# Step 13: Optional - Variable importance
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "RE Variable Importance") +
  theme_classic(base_size = 14)

# Step 14: Print summary metrics
collect_metrics(logit_res)
collect_metrics(rf_res) 
```

\newpage

## 4.2. Relaxed Classification

As a robustness check, we re-estimated our supervised learning validation using an alternative classification scheme — `VDEM_STATUS_BASELINE` — which adopts a narrower definition of the "Transitioning" category by labeling as `Fragile` all country-years below the 20th percentile of the V-Dem fragility index and as `Non-Fragile` those above the 25th percentile. Observations between these thresholds are coded as `Transitioning`. This specification provides a more relaxed test of the classification’s internal validity by increasing the conceptual and empirical distance between adjacent categories.

Under this baseline setup, both the multinomial logistic regression and the random forest models continued to perform exceptionally well. The logit model achieved an out-of-sample classification accuracy of **93.7%**, with a multiclass log-loss of **0.40** and a standard error of **0.0032**, indicating strong predictive power and low variance across folds. The random forest model followed closely with an accuracy of **93.3%** and an even lower log-loss of **0.24**, suggesting better-calibrated probabilistic predictions. As with the ideal specification, the confusion matrices reveal that most misclassifications are concentrated along adjacent categories — particularly between `Fragile` and `Transitioning`, and between `Transitioning` and `Non-Fragile` — with virtually no confusion between `Fragile` and `Non-Fragile`. This pattern reinforces the ordinal logic and substantive coherence of the classification.

Taken together, the results suggest that our classification thresholds remain empirically robust even under narrower definitions of fragility and transition. The models' ability to recover the hand-coded labels using only structural predictors confirms that the classification captures meaningful differences in regime characteristics, conflict exposure, and institutional performance across country-years.

```{r, echo = FALSE, message = FALSE, warning = FALSE}

# Set seed for reproducibility
set.seed(247)

# Step 1: Define predictors
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

# Step 2: Build dataset
df <- standardized_data_final |>
  select(all_of(predictors_20), VDEM_STATUS_BASELINE) |>
  filter(!is.na(VDEM_STATUS_BASELINE)) |>
  drop_na() |>
  mutate(VDEM_STATUS_BASELINE = as.factor(VDEM_STATUS_BASELINE))

# Step 3: Train/test split
data_split <- initial_split(df, strata = VDEM_STATUS_BASELINE)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Step 4: Recipe
rec <- recipe(VDEM_STATUS_BASELINE ~ ., data = train_data) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# Step 5: Model specs
logit_spec <- multinom_reg(mode = "classification") |> set_engine("nnet")
rf_spec <- rand_forest(mode = "classification", trees = 500) |>
  set_engine("ranger", importance = "permutation")

# Step 6: Workflows
logit_wf <- workflow() |> add_model(logit_spec) |> add_recipe(rec)
rf_wf    <- workflow() |> add_model(rf_spec)    |> add_recipe(rec)

# Step 7: CV setup
folds <- vfold_cv(train_data, v = 10, repeats = 3, strata = VDEM_STATUS_BASELINE)

# Step 8: CV model evaluation
logit_res <- fit_resamples(logit_wf, folds, metrics = metric_set(accuracy, mn_log_loss))
rf_res    <- fit_resamples(rf_wf, folds, metrics = metric_set(accuracy, mn_log_loss))

# Step 9: Final model fits
logit_fit <- fit(logit_wf, data = train_data)
rf_fit    <- fit(rf_wf, data = train_data)

# Step 10: Predict
logit_preds <- predict(logit_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_BASELINE)

rf_preds <- predict(rf_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_BASELINE)

# Step 11: Compute confusion matrices + accuracy
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds, truth = truth, estimate = .pred_class)

accuracy_logit <- accuracy(logit_preds, truth = truth, estimate = .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds, truth = truth, estimate = .pred_class)$.estimate

# Step 12: Plot confusion matrix with custom colors + accuracy in title
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(
    title = paste0("Logit Confusion Matrix (Accuracy = ", round(accuracy_logit * 100, 1), "%)"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(
    title = paste0("RE Confusion Matrix (Accuracy = ", round(accuracy_rf * 100, 1), "%)"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

# Step 13: Optional - Variable importance
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "RE Variable Importance") +
  theme_classic(base_size = 14)

# Step 14: Print summary metrics
collect_metrics(logit_res)
collect_metrics(rf_res)
```

\newpage

## 4.3. Strict Classification

Lastly, we re-ran the validation exercise using an alternative, stricter scheme: countries are classified as `Fragile` if they fall below the 20th percentile of the normalized V-Dem latent factor, as `Non-Fragile` if above the 40th percentile, and as `Transitioning` in between (20th–40th percentiles). As before, we trained two supervised machine learning models — a multinomial logistic regression and a random forest — using the same 20 structural predictors capturing regime type, democratic institutions, economic performance, and conflict intensity.

Both models demonstrate high predictive performance under this revised definition. The **multinomial logit model** achieved an average classification accuracy of **91.8%** with a multiclass log-loss of **0.36**, and low standard errors (0.0029 and 0.0263, respectively) across 30 cross-validation folds. The **random forest model** performed similarly, with slightly lower accuracy at **90.4%**, but a significantly improved log-loss of **0.31**, suggesting more calibrated probability estimates. Confusion matrices for both models show that most misclassifications continue to occur between adjacent conceptual categories (e.g., `Fragile` vs. `Transitioning`), with few cases where `Fragile` is mistaken for `Non-Fragile`, and vice versa.

These findings reinforce the internal validity of our classification strategy. Even under a more demanding definition, the models are able to accurately recover the three-tiered labels based purely on structural indicators, without relying on the latent V-Dem factor. This underscores that our classification is not arbitrary or circular, but rather captures empirically discernible patterns in institutional and socio-political fragility.

```{r, echo = FALSE, message = FALSE, warning = FALSE}

# Set seed for reproducibility
set.seed(247)

# Step 1: Define predictors
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

# Building our dataset
df <- standardized_data_final |>
  select(all_of(predictors_20), 
         VDEM_STATUS_IDEAL) |>
  filter(!is.na(VDEM_STATUS_IDEAL)) |>
  drop_na() |>
  mutate(VDEM_STATUS_IDEAL = as.factor(VDEM_STATUS_IDEAL))

# Train/test split
data_split <- initial_split(df, strata = VDEM_STATUS_IDEAL)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Recipe
rec <- recipe(VDEM_STATUS_IDEAL ~ ., data = train_data) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# Model specifications ----------------------------------------------------------

# Logit Model 
logit_spec <- multinom_reg(mode = "classification") |> 
  set_engine("nnet")

# Random Forest Model 
rf_spec <- rand_forest(mode = "classification", 
                       trees = 500) |>
  set_engine("ranger", 
             importance = "permutation")

# Step 6: Workflows
logit_wf <- workflow() |> add_model(logit_spec) |> add_recipe(rec)
rf_wf    <- workflow() |> add_model(rf_spec)    |> add_recipe(rec)

# Step 7: CV setup
folds <- vfold_cv(train_data, v = 10, repeats = 3, strata = VDEM_STATUS_IDEAL)

# Step 8: CV model evaluation
logit_res <- fit_resamples(logit_wf, folds, metrics = metric_set(accuracy, mn_log_loss))
rf_res    <- fit_resamples(rf_wf, folds, metrics = metric_set(accuracy, mn_log_loss))

# Step 9: Final model fits
logit_fit <- fit(logit_wf, data = train_data)
rf_fit    <- fit(rf_wf, data = train_data)

# Step 10: Predict
logit_preds <- predict(logit_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

rf_preds <- predict(rf_fit, test_data, type = "class") |>
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

# Step 11: Compute confusion matrices + accuracy
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds, truth = truth, estimate = .pred_class)

accuracy_logit <- accuracy(logit_preds, truth = truth, estimate = .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds, truth = truth, estimate = .pred_class)$.estimate

# Step 12: Plot confusion matrix with custom colors + accuracy in title
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(
    title = paste0("Logit Confusion Matrix (Accuracy = ", round(accuracy_logit * 100, 1), "%)"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(
    title = paste0("RE Confusion Matrix (Accuracy = ", round(accuracy_rf * 100, 1), "%)"),
    fill = "Count"
  ) +
  theme_minimal(base_size = 14)

# Step 13: Optional - Variable importance
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "RE Variable Importance") +
  theme_classic(base_size = 14)

# Step 14: Print summary metrics
collect_metrics(logit_res)
collect_metrics(rf_res) 
```

## 4.4. Conclusion: Validating Fragility Classifications

The analyses presented in this section provide strong empirical validation for our proposed fragility classification strategy. Across all three specifications — *ideal* (30–35th percentiles), *relaxed* (20–25th), and *strict* (20–40th) — the use of supervised machine learning algorithms demonstrated that our tiered classification (`Fragile`, `Transitioning` or `Non-Fragile`) is not only conceptually sound but also empirically recoverable from a set of 20 structural predictors reflecting regime type, democratic quality, institutional durability, economic performance, and conflict exposure.

Both multinomial logistic regression and random forest models performed consistently well across specifications. The ideal model achieved over 91% accuracy in logit and 89% in random forest, with misclassifications concentrated along adjacent conceptual categories (e.g., between `Fragile` and `Transitioning`). Even when the classification thresholds were narrowed (relaxed) or broadened (stricted), predictive performance remained stable: accuracy ranged from 90% to 94%, and log-loss values decreased as the spacing between categories widened—indicating better probabilistic calibration in more polarized cases.

The key takeaway is that our classification system is **not arbitrary**. The fact that simple models trained only on exogenous structural variables can consistently reproduce the labels we defined using V-Dem latent scores suggests that those scores capture real, observable dimensions of state fragility. Moreover, the consistent structure of misclassifications — rarely jumping from `Fragile` to `Non-Fragile` or vice versa — confirms that the Transitioning category reflects a meaningful, intermediate empirical space.

These results support the continued use of this supervised classification framework. At this stage, the use of **unsupervised machine learning methods** (e.g., clustering, dimensionality reduction) does not appear necessary to re-derive the classification from scratch, as our labels already reflect an interpretable and externally valid mapping of fragility status across country-years. We are thus well-positioned to use these classifications as inputs into downstream empirical analyses.

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

```{r, echo = FALSE, message = FALSE, warning = FALSE}

###############################################################################
# Fragility-tier prediction with *country-blocked* CV
# --------------------------------------------------

# ---------------------------------------------------------------------------
# 1. Predictor list (20 structural & conflict variables)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# 2. Build modelling data  (country & year kept as IDs)
# ---------------------------------------------------------------------------
df <- standardized_data_final |>
  select(COUNTRY_NAME, 
         YEAR,
         all_of(predictors_20), 
         VDEM_STATUS_IDEAL) |> 
  filter(!is.na(VDEM_STATUS_IDEAL)) |> 
  drop_na() |>
  mutate(VDEM_STATUS_IDEAL = as_factor(VDEM_STATUS_IDEAL)) |>
  # a constant, country-level label for stratification
  group_by(COUNTRY_NAME) |>
  mutate(country_mode = names(which.max(table(VDEM_STATUS_IDEAL)))) |> 
  ungroup()

# ---------------------------------------------------------------------------
# 3. Country-blocked 80/20 train-test split  (stratify on country_mode)
# ---------------------------------------------------------------------------
data_split <- group_initial_split(
  df,
  group  = COUNTRY_NAME,
  prop   = 0.80,
  strata = country_mode            # safe: constant within group
)

train_data <- training(data_split)

test_data  <- testing(data_split)

# ---------------------------------------------------------------------------
# 4. Pre-processing recipe  (upsample minority tier *inside each resample*)
# ---------------------------------------------------------------------------
rec <- recipe(VDEM_STATUS_IDEAL ~ ., data = train_data) |> 
  update_role(COUNTRY_NAME, YEAR, country_mode, new_role = "id") %>% 
  step_zv(all_predictors()) %>% 
  # -------- imbalance fix --------
  step_upsample(VDEM_STATUS_IDEAL, over_ratio = 1) %>%  # OR step_smote()
  # --------------------------------
  step_normalize(all_numeric_predictors())

# ---------------------------------------------------------------------------
# 5. Model specifications  (class-weighted RF)
# ---------------------------------------------------------------------------
## 5a. Multinomial logit  — learns on the upsampled data
logit_spec <- multinom_reg(mode = "classification") %>% 
  set_engine("nnet")

## 5b. Random forest  — add inverse-frequency weights
cls_freq    <- table(train_data$VDEM_STATUS_IDEAL)
rf_weights  <- as.numeric(max(cls_freq) / cls_freq)  # named numeric
names(rf_weights) <- names(cls_freq)

rf_spec <- rand_forest(mode = "classification", trees = 500) %>% 
  set_engine("ranger",
             importance    = "permutation",
             class.weights = rf_weights)

# ---------------------------------------------------------------------------
# 6. Workflows
# ---------------------------------------------------------------------------
logit_wf <- workflow() %>% add_model(logit_spec) %>% add_recipe(rec)
rf_wf    <- workflow() %>% add_model(rf_spec)    %>% add_recipe(rec)

# ---------------------------------------------------------------------------
# 7. 10-fold, 3-repeat grouped CV  (countries intact; class balance via mode)
# ---------------------------------------------------------------------------
folds <- group_vfold_cv(
  train_data,
  group   = COUNTRY_NAME,
  v       = 10,
  repeats = 3,
  strata  = country_mode
)

# ---------------------------------------------------------------------------
# 8. Metric set  – include macro recall & macro F1
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 8. Metric set  – include macro recall & macro F1
# ---------------------------------------------------------------------------
library(tidymodels)   # loads yardstick
# detach caret if it’s still masking recall:
if ("package:caret" %in% search()) detach("package:caret", unload = TRUE)

metrics_set <- yardstick::metric_set(
  yardstick::accuracy,        # unambiguous
  yardstick::mn_log_loss,     # unambiguous

  # ----- macro recall -----
  yardstick::metric_tweak(
    "recall_macro",           # ← .name
    yardstick::recall,        # ← .fn (the metric function)
    estimator = "macro"       # override default
  ),

  # ----- macro F₁ -----
  yardstick::metric_tweak(
    "f_meas_macro",           # .name
    yardstick::f_meas,        # .fn
    estimator = "macro"
  )
)

# ---------------------------------------------------------------------------
# 9. Cross-validated training & evaluation
# ---------------------------------------------------------------------------
logit_res <- fit_resamples(
  logit_wf,
  resamples = folds,
  metrics   = metrics_set,
  control   = control_resamples(save_pred = TRUE)
)
rf_res <- fit_resamples(
  rf_wf,
  resamples = folds,
  metrics   = metrics_set,
  control   = control_resamples(save_pred = TRUE)
)

# ---------------------------------------------------------------------------
# 10. Final fits on *training* data
# ---------------------------------------------------------------------------
logit_fit <- fit(logit_wf, data = train_data)
rf_fit    <- fit(rf_wf,    data = train_data)

# ---------------------------------------------------------------------------
# 11. Predict on the *held-out* test countries
# ---------------------------------------------------------------------------
logit_preds <- predict(logit_fit, test_data, type = "class") %>% 
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)
rf_preds <- predict(rf_fit, test_data, type = "class") %>% 
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

# ---------------------------------------------------------------------------
# 12. Confusion matrices & overall accuracy
# ---------------------------------------------------------------------------
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds,    truth = truth, estimate = .pred_class)

accuracy_logit <- accuracy(logit_preds, truth, .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds,    truth, .pred_class)$.estimate

# visualise
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(title = paste0("ML – Held-out Countries (Acc = ",
                      round(accuracy_logit*100,1), "%)"),
       fill  = "Count") +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(title = paste0("RF – Held-out Countries (Acc = ",
                      round(accuracy_rf*100,1), "%)"),
       fill  = "Count") +
  theme_minimal(base_size = 14)

# ---------------------------------------------------------------------------
# 13. Variable-importance plot (RF)
# ---------------------------------------------------------------------------
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "Random Forest – Permutation Importance") +
  theme_classic(base_size = 14)

# ---------------------------------------------------------------------------
# 14. Cross-validated summary metrics
# ---------------------------------------------------------------------------
collect_metrics(logit_res)
collect_metrics(rf_res)
```

### 4.5.a. Model with Imputation and without Highly Correlated Predictors 

```{r, echo = FALSE, message = FALSE, warning = FALSE}

# ---------------------------------------------------------------------------
# 2. Build modelling data  (country & year kept as IDs)
# ---------------------------------------------------------------------------
df <- standardized_data_final |>
  select(COUNTRY_NAME, 
         YEAR,
         all_of(predictors_20), 
         VDEM_STATUS_IDEAL) |> 
  filter(!is.na(VDEM_STATUS_IDEAL)) |> 
  drop_na() |>
  mutate(VDEM_STATUS_IDEAL = as_factor(VDEM_STATUS_IDEAL)) |>
  # a constant, country-level label for stratification
  group_by(COUNTRY_NAME) |>
  mutate(country_mode = names(which.max(table(VDEM_STATUS_IDEAL)))) |> 
  ungroup()

# ---------------------------------------------------------------------------
# 3. Country-blocked 80/20 train-test split  (stratify on country_mode)
# ---------------------------------------------------------------------------
data_split <- group_initial_split(
  df,
  group  = COUNTRY_NAME,
  prop   = 0.80,
  strata = country_mode            # safe: constant within group
)

train_data <- training(data_split)

test_data  <- testing(data_split)

# ---------------------------------------------------------------------------
# 4. Pre-processing recipe  (upsample minority tier *inside each resample*)
# ---------------------------------------------------------------------------
rec <- recipe(VDEM_STATUS_IDEAL ~ ., data = train_data) |> 
  update_role(COUNTRY_NAME, YEAR, country_mode, new_role = "id") %>% 
  step_zv(all_predictors()) %>% 
  # -------- imbalance fix --------
  step_upsample(VDEM_STATUS_IDEAL, over_ratio = 1) %>%  # OR step_smote()
  # --------------------------------
  step_normalize(all_numeric_predictors())

# ---------------------------------------------------------------------------
# 5. Model specifications  (class-weighted RF)
# ---------------------------------------------------------------------------
## 5a. Multinomial logit  — learns on the upsampled data
logit_spec <- multinom_reg(mode = "classification") %>% 
  set_engine("nnet")

## 5b. Random forest  — add inverse-frequency weights
cls_freq    <- table(train_data$VDEM_STATUS_IDEAL)
rf_weights  <- as.numeric(max(cls_freq) / cls_freq)  # named numeric
names(rf_weights) <- names(cls_freq)

rf_spec <- rand_forest(mode = "classification", trees = 500) %>% 
  set_engine("ranger",
             importance    = "permutation",
             class.weights = rf_weights)

# ---------------------------------------------------------------------------
# 6. Workflows
# ---------------------------------------------------------------------------
logit_wf <- workflow() %>% add_model(logit_spec) %>% add_recipe(rec)
rf_wf    <- workflow() %>% add_model(rf_spec)    %>% add_recipe(rec)

# ---------------------------------------------------------------------------
# 7. 10-fold, 3-repeat grouped CV  (countries intact; class balance via mode)
# ---------------------------------------------------------------------------
folds <- group_vfold_cv(
  train_data,
  group   = COUNTRY_NAME,
  v       = 10,
  repeats = 3,
  strata  = country_mode
)

# ---------------------------------------------------------------------------
# 8. Metric set  – include macro recall & macro F1
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 8. Metric set  – include macro recall & macro F1
# ---------------------------------------------------------------------------
library(tidymodels)   # loads yardstick
# detach caret if it’s still masking recall:
if ("package:caret" %in% search()) detach("package:caret", unload = TRUE)

metrics_set <- yardstick::metric_set(
  yardstick::accuracy,        # unambiguous
  yardstick::mn_log_loss,     # unambiguous

  # ----- macro recall -----
  yardstick::metric_tweak(
    "recall_macro",           # ← .name
    yardstick::recall,        # ← .fn (the metric function)
    estimator = "macro"       # override default
  ),

  # ----- macro F₁ -----
  yardstick::metric_tweak(
    "f_meas_macro",           # .name
    yardstick::f_meas,        # .fn
    estimator = "macro"
  )
)

# ---------------------------------------------------------------------------
# 9. Cross-validated training & evaluation
# ---------------------------------------------------------------------------
logit_res <- fit_resamples(
  logit_wf,
  resamples = folds,
  metrics   = metrics_set,
  control   = control_resamples(save_pred = TRUE)
)
rf_res <- fit_resamples(
  rf_wf,
  resamples = folds,
  metrics   = metrics_set,
  control   = control_resamples(save_pred = TRUE)
)

# ---------------------------------------------------------------------------
# 10. Final fits on *training* data
# ---------------------------------------------------------------------------
logit_fit <- fit(logit_wf, data = train_data)
rf_fit    <- fit(rf_wf,    data = train_data)

# ---------------------------------------------------------------------------
# 11. Predict on the *held-out* test countries
# ---------------------------------------------------------------------------
logit_preds <- predict(logit_fit, test_data, type = "class") %>% 
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)
rf_preds <- predict(rf_fit, test_data, type = "class") %>% 
  bind_cols(truth = test_data$VDEM_STATUS_IDEAL)

# ---------------------------------------------------------------------------
# 12. Confusion matrices & overall accuracy
# ---------------------------------------------------------------------------
conf_mat_logit <- conf_mat(logit_preds, truth = truth, estimate = .pred_class)
conf_mat_rf    <- conf_mat(rf_preds,    truth = truth, estimate = .pred_class)

accuracy_logit <- accuracy(logit_preds, truth, .pred_class)$.estimate
accuracy_rf    <- accuracy(rf_preds,    truth, .pred_class)$.estimate

# visualise
autoplot(conf_mat_logit, type = "heatmap") +
  scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
  labs(title = paste0("ML – Held-out Countries (Acc = ",
                      round(accuracy_logit*100,1), "%)"),
       fill  = "Count") +
  theme_minimal(base_size = 14)

autoplot(conf_mat_rf, type = "heatmap") +
  scale_fill_gradient(low = "#fee5d9", high = "#cb181d") +
  labs(title = paste0("RF – Held-out Countries (Acc = ",
                      round(accuracy_rf*100,1), "%)"),
       fill  = "Count") +
  theme_minimal(base_size = 14)

# ---------------------------------------------------------------------------
# 13. Variable-importance plot (RF)
# ---------------------------------------------------------------------------
vip(extract_fit_parsnip(rf_fit), num_features = 20, geom = "col") +
  labs(title = "Random Forest – Permutation Importance") +
  theme_classic(base_size = 14)

# ---------------------------------------------------------------------------
# 14. Cross-validated summary metrics
# ---------------------------------------------------------------------------
collect_metrics(logit_res)
collect_metrics(rf_res)
```


# 5. Linear and Logit Binomial Models 

## 5.1. Linear and Logit Binomial Coefficients

* **Ordinary Least Square Coefficients (OLS):** The model fits the data into a regression line and estimates the change in the outcome variable due to a one unit increase in each predictor. For dichotomous outcomes (`1 = Fragile` or `0 = Non Fragile`), the coefficient estimates the average difference between the two levels of the dependent variable. The OLS model is - however - probably ill suited in this case as it assumes a continuous distribution of the dependent variable and may end up producing predicted values outside the [0, 1].  

   + The *standard errors* (SEs), reported below each coefficient, capture the variability of the same coefficient (i.e., on average, how far from the regression line observed values fall). Naturally, the smaller the SEs, the better. 
   
   + SEs are "clustered" at the country level, as it is reasonably to assume that observations for the same country are likely to be correlated over time (autocorrelation). The robust standard errors reported also account for heteroskedasticity (i.e., the non-constant variance of errors).

* **Logistic Regression (Logit) Coefficients:** The coefficients suggest how many standard deviations the target variable (i.e., being `Fragile` or `Non-Fragile` in a given year; the log-odds of the outcome in logistic regression) changes per standard deviation change in the predictor variable. 

* **Performance:** The table returns a couple of metrics of performance ($R^2$ versus $PseudoR^2$, $AIC$ and $BIC$) that allows for direct comparisons of models' performances. The plots estimate the ROC (Receiver Operating Characteristic) curve, which identifies how well the model can distinguish between the two classes (e.g., `Fragile` and `Non-Fragile` outcomes) across various thresholds. The Area Under the Curve (AUC) the quality of this binary classification model: the higher the AUC, the better!

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ------------------------------------------------------------
# 1. Build RHS string (typos intact) + year FE ---------------
# ------------------------------------------------------------
baseline_vars <- c(
  "COMBINED_POLITY_SCORE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "GDP_GROWTH",
  "MAX_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",   # <-- typo kept intentionally
  "INSTITUTIONAL_DEMOCRACY_SOCRE",         # <-- typo kept intentionally
  "LIBERAL_DEMOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "GDP_PER_CAPITA_GROWTH",
  "POLITICAL_REGIME",
  "N_WAR_FRONTS",
  "ODA_RECEIVED_PER_CAPITA",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "ELECTORAL_DEMOCRACY_SCORE",
  "POLITICAL_COMPETITION_SCORE",
  "CONFLICT_INTENSITY_YEAR",
  "AVG_CONFLICT_INTENSITY",
  "GDP_DEFLATOR",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM"
)

form_baseline <- as.formula(
  paste("VDEM_FRAGILE_IDEAL ~", paste(baseline_vars, collapse = " + "))
)

# ------------------------------------------------------------
# 2. OLS (cluster-robust at COUNTRY_NAME) --------------------
# ------------------------------------------------------------
linear_model <- lm(form_baseline, data = standardized_data_final)
clustered_se_linear <- sqrt(diag(vcovCL(linear_model, cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 3. Logit (cluster-robust at COUNTRY_NAME) ------------------
# ------------------------------------------------------------
logit_model  <- glm(form_baseline,
                    family = binomial(link = "logit"),
                    data   = standardized_data_final)
clustered_se_logit <- sqrt(diag(vcovCL(logit_model, cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 4. Fixed-effects (within) + Driscoll–Kraay SEs -------------
# ------------------------------------------------------------
fixed_effects_model <- plm(form_baseline,
                           data  = standardized_data_final,
                           model = "within",
                           index = c("ISO_CODE_3", "YEAR"))
fe_dk_se <- sqrt(diag(vcovSCC(fixed_effects_model,
                              maxlag = 3, type = "HC1")))

# ------------------------------------------------------------
# 5. Fit statistics ------------------------------------------
# ------------------------------------------------------------
R2_linear <- summary(linear_model)$r.squared

logit_null <- glm(VDEM_FRAGILE_BASELINE ~ 1,
                  family = binomial,
                  data   = standardized_data_final)
McFadden_R2 <- 1 - (as.numeric(logLik(logit_model)) /
                    as.numeric(logLik(logit_null)))

add_stats <- list(
  c("R-squared (OLS)",          sprintf("%.3f", R2_linear),      ""),
  c("Pseudo R-squared (Logit)", "",                              sprintf("%.3f", McFadden_R2)),
  c("AIC",                      sprintf("%.2f", AIC(linear_model)), sprintf("%.2f", AIC(logit_model))),
  c("BIC",                      sprintf("%.2f", BIC(linear_model)), sprintf("%.2f", BIC(logit_model)))
)
```

```{r, linear-binomial, echo = FALSE, message = FALSE, warning = FALSE, results = 'asis'}
# ------------------------------------------------------------
# 6. LaTeX table (labels unchanged, year FE omitted) ---------
# ------------------------------------------------------------
stargazer(linear_model,
          logit_model,
          fixed_effects_model,
          se = list(clustered_se_linear,
                    clustered_se_logit,
                    fe_dk_se),
          covariate.labels = c("Combined polity score",
                               "Institutional autocracy score",
                               "GDP growth",
                               "MAX conflict intensity",
                               "N total troops involved",
                               "Territorial fragmentation",
                               "Executive recruitment",
                               "Institutional democracy score",
                               "Liberal democracy score",
                               "Regime durability in years",
                               "GDP per capita growth",
                               "Political regime",
                               "N war of fronts",
                               "ODA received per capita",
                               "Conflitc cumul. intensity across years",
                               "Electoral democracy score",
                               "Political competition score",
                               "Conflict intensity year",
                               "AVG conflict intensity",
                               "GDP deflator",
                               "Partial democracy with factionalism",
                               "Constant"),
          omit = "factor\\(YEAR\\)",     # hide year dummies
          header = FALSE,
          no.space = TRUE,
          font.size = "small",
          add.lines = add_stats,
          title = "Linear, Logistic and FE Estimates of Fragility Status",
          dep.var.caption = "Fragility Status",
          dep.var.labels = "")
```

\newpage

## 5.2 Multicolinearity 

Multicollinearity refers to a situation in which two or more predictor variables in a regression model are highly correlated, meaning they provide overlapping or redundant information. While multicollinearity does not bias coefficient estimates per se, it inflates standard errors, reduces statistical power, and complicates interpretation — especially in models aiming to isolate the marginal effect of individual predictors.

This problem is particularly relevant in our context, where many institutional and democracy indicators are theoretically and empirically interrelated (e.g., different components of executive recruitment, regime durability, and democracy scores often move together).

Multicollinearity matters differently across models:

- **Ordinary Least Squares (OLS):**
  - Inflated standard errors make it difficult to determine which predictors are truly significant;
  - Coefficient signs and magnitudes may become unstable or counterintuitive;
  - The overall model fit ($R^2$) may remain high even when individual predictors appear insignificant.

- **Logistic Regression (Logit):**
  - The same inflation of standard errors occurs, obscuring the relationship between predictors and outcome;
  - Convergence issues may arise during model fitting if predictors are highly collinear;
  - Interpretation of odds ratios or marginal effects becomes less reliable when predictors are entangled.

- **Country-Year Fixed Effects Models:**
  - With fixed effects absorbing between-country and/or between-year variation, there's less "within" variation left to estimate the effects of collinear predictors;
  - Perfect or near-perfect collinearity can lead to dropped variables or singularity warnings;
  - It exacerbates overfitting and can reduce the degrees of freedom substantially in already high-dimensional panel settings.

In sum, diagnosing and addressing multicollinearity is essential to ensure that our estimated effects are stable, interpretable, and not artifacts of redundancy in the predictor set.

```{r, echo = FALSE, message = FALSE, warning = FALSE}
library(car)
library(dplyr)
library(knitr)
library(kableExtra)

# ---- 1. Re-estimate the LPM WITHOUT the year dummies just for VIF ----------
formula_no_year <- update(form_baseline, . ~ . - factor(YEAR))
lpm_no_year     <- lm(formula_no_year, data = standardized_data_final)

# ---- 2. Compute VIFs -------------------------------------------------------
raw_vif <- car::vif(lpm_no_year)

# Convert to a tidy data frame, handling both vector and matrix cases
vif_tbl <- if (is.matrix(raw_vif)) {
  data.frame(
    Variable = rownames(raw_vif),
    VIF      = raw_vif[, "GVIF^(1/(2*Df))", drop = TRUE]
  )
} else {
  data.frame(
    Variable = names(raw_vif),
    VIF      = as.numeric(raw_vif)
  )
}

vif_tbl <- vif_tbl %>% arrange(desc(VIF))

# ---- 3. Print, highlighting VIF > 10 --------------------------------------
kable(vif_tbl, digits = 2,
      caption = "Variance–Inflation Factors for Baseline Predictors (Year FE removed for diagnostics)") %>%
  kable_styling(latex_options = c("hold_position", "striped")) %>%
  row_spec(which(vif_tbl$VIF > 10),
           bold = TRUE, color = "white", background = "#D7261E")
```

## 5.3. Colinearity Heatmap

```{r, echo = FALSE, message = FALSE, warning=FALSE}

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 1. Predictor whitelist (typos intact) ---------------------------------------
predictor_vars <- c(
  "COMBINED_POLITY_SCORE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "GDP_GROWTH",
  "MAX_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",
  "LIBERAL_DEMOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "GDP_PER_CAPITA_GROWTH",
  "POLITICAL_REGIME",
  "N_WAR_FRONTS",
  "ODA_RECEIVED_PER_CAPITA",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "ELECTORAL_DEMOCRACY_SCORE",
  "POLITICAL_COMPETITION_SCORE",
  "CONFLICT_INTENSITY_YEAR",
  "AVG_CONFLICT_INTENSITY",
  "GDP_DEFLATOR",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM"
)

# 2. Identify which of those columns are numeric ------------------------------
numeric_vars <- predictor_vars[
  sapply(standardized_data_final[predictor_vars], is.numeric)
]

non_numeric <- setdiff(predictor_vars, numeric_vars)
if (length(non_numeric) > 0) {
  message("Dropped non-numeric variables from heat-map: ",
          paste(non_numeric, collapse = ", "))
}

# 3. Compute correlation matrix (ungroup first!) ------------------------------
corr_mat <- standardized_data_final %>%
  ungroup() %>%                         # remove COUNTRY_NAME grouping
  select(all_of(numeric_vars)) %>%
  cor(use = "pairwise.complete.obs")

# 4. Re-order with hierarchical clustering for clearer blocks -----------------
hc <- hclust(as.dist(1 - abs(corr_mat)))
corr_mat <- corr_mat[hc$order, hc$order]

# 5. Long format for ggplot ----------------------------------------------------
corr_long <- as.data.frame(as.table(corr_mat)) %>%
  rename(Var1 = Var1, Var2 = Var2, r = Freq)

# 6. Heat-map ------------------------------------------------------------------
ggplot(corr_long, aes(Var1, Var2, fill = r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#8B0000", mid = "white", high = "#00468B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Correlation") +
  coord_fixed() +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank()) +
  labs(title = "Pairwise Correlations among Baseline Numeric Predictors",
       x = "", y = "")
```

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# Load dplyr just in case
library(dplyr)

# ---------------------------
# 1. PCA: Democracy indicators
# ---------------------------
democracy_vars <- c("INSTITUTIONAL_DEMOCRACY_SOCRE",
                    "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
                    "LIBERAL_DEMOCRACY_SCORE")

# Extract clean matrix (force numeric and drop malformed)
democracy_mat <- standardized_data_final %>%
  ungroup() %>%
  select(all_of(democracy_vars)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))

# Drop rows with any NA or Inf
democracy_clean <- democracy_mat %>%
  filter(if_all(everything(), ~ !is.na(.) & is.finite(.)))

# Run PCA
pca_democracy <- prcomp(democracy_clean, scale. = TRUE)

# Project PC1 scores back
PC1_demo <- predict(pca_democracy, newdata = democracy_mat)[, 1]
standardized_data_final$PC_DEMOCRACY <- as.numeric(scale(PC1_demo))

# Optional diagnostics
print("PC_DEMOCRACY loadings:")
print(round(pca_democracy$rotation[, 1], 3))
print("Explained variance:")
print(summary(pca_democracy))


# ---------------------------
# 2. PCA: Conflict indicators
# ---------------------------
# 1. Pull and flatten conflict variables
conflict_vars <- c(
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS"
)

# Extract clean matrix (force numeric and drop malformed)
conflict_mat <- standardized_data_final %>%
  ungroup() %>%
  select(all_of(conflict_vars)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))

# Drop rows with any NA or Inf
conflict_clean <- conflict_mat %>%
  filter(if_all(everything(), ~ !is.na(.) & is.finite(.)))

# Run PCA
pca_conflict <- prcomp(conflict_clean, scale. = TRUE)

# Project PC1 scores back
PC1_conf <- predict(pca_conflict, newdata = conflict_mat)[, 1]
standardized_data_final$PC_CONFLICT <- as.numeric(scale(PC1_conf))

# Optional diagnostics
print("PC_CONFLICT loadings:")
print(round(pca_conflict$rotation[, 1], 3))
print("Explained variance:")
print(summary(pca_conflict))



```

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ------------------------------------------------------------
# 1. Build final predictor list (typos intact + PCA + drop GDP_GROWTH)
# ------------------------------------------------------------
baseline_vars_pca <- c(
  "COMBINED_POLITY_SCORE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "N_TOTAL_TROOPS",
  "TERRITORIAL_FRAGMENTATION",
  "REGIME_DURABILITY_YEARS",
  "GDP_PER_CAPITA_GROWTH",
  "POLITICAL_REGIME",
  "N_WAR_FRONTS",
  "ODA_RECEIVED_PER_CAPITA",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "ELECTORAL_DEMOCRACY_SCORE",
  "POLITICAL_COMPETITION_SCORE",
  "GDP_DEFLATOR",
  # Replacements via PCA
  "PC_DEMOCRACY",
  "PC_CONFLICT"
)

form_baseline_pca <- as.formula(
  paste("VDEM_FRAGILE_BASELINE ~", paste(baseline_vars_pca, collapse = " + "), "+ factor(YEAR)")
)

# ------------------------------------------------------------
# 2. OLS (cluster-robust at COUNTRY_NAME)
# ------------------------------------------------------------
linear_model_pca <- lm(form_baseline_pca, data = standardized_data_final)
clustered_se_linear_pca <- sqrt(diag(vcovCL(linear_model_pca, cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 3. Logit (cluster-robust at COUNTRY_NAME)
# ------------------------------------------------------------
logit_model_pca <- glm(form_baseline_pca,
                       family = binomial(link = "logit"),
                       data   = standardized_data_final)
clustered_se_logit_pca <- sqrt(diag(vcovCL(logit_model_pca, cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 4. Fixed-effects (within) + Driscoll–Kraay SEs
# ------------------------------------------------------------
fixed_effects_model_pca <- plm(form_baseline_pca,
                               data  = standardized_data_final,
                               model = "within",
                               index = c("ISO_CODE_3", 
                                         "YEAR"))

fe_dk_se_pca <- sqrt(diag(vcovSCC(fixed_effects_model_pca, maxlag = 3, type = "HC1")))

# ------------------------------------------------------------
# 5. Fit statistics
# ------------------------------------------------------------
R2_linear_pca <- summary(linear_model_pca)$r.squared

logit_null_pca <- glm(VDEM_FRAGILE_BASELINE ~ 1,
                      family = binomial,
                      data   = standardized_data_final)

McFadden_R2_pca <- 1 - (as.numeric(logLik(logit_model_pca)) /
                        as.numeric(logLik(logit_null_pca)))

add_stats <- list(
  c("R-squared (OLS)",          sprintf("%.3f", R2_linear_pca),      ""),
  c("Pseudo R-squared (Logit)", "",                              sprintf("%.3f", McFadden_R2_pca)),
  c("AIC",                      sprintf("%.2f", AIC(linear_model_pca)), sprintf("%.2f", AIC(logit_model_pca))),
  c("BIC",                      sprintf("%.2f", BIC(linear_model_pca)), sprintf("%.2f", BIC(logit_model_pca)))
)
```


```{r, linear-binomial-pca, echo = FALSE, message = FALSE, warning = FALSE, results = 'asis'}
stargazer(linear_model_pca,
          logit_model_pca,
          fixed_effects_model_pca,
          se = list(clustered_se_linear_pca,
                    clustered_se_logit_pca,
                    fe_dk_se_pca),
          covariate.labels = c(
            "Combined polity score",
            "Institutional autocracy score",
            "N total troops involved",
            "Territorial fragmentation",
            "Regime durability in years",
            "GDP per capita growth",
            "Political regime",
            "N war of fronts",
            "ODA received per capita",
            "Partial democracy with factionalism",
            "Electoral democracy score",
            "Political competition score",
            "GDP deflator",
            "Institutional + liberal democracy (PC)",
            "Conflict intensity composite (PC)",
            "Constant"
          ),
          omit = "factor\\(YEAR\\)",   # keep in model but omit from table
          header = FALSE,
          no.space = TRUE,
          font.size = "small",
          add.lines = add_stats,
          title = "Linear, Logistic and FE Estimates of Fragility Status wit PCAs",
          dep.var.caption = "Fragility Status",
          dep.var.labels = "")
```

# Appendix A

## A.1. Logit and Binomial Coefficients: Relaxed Classification

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ------------------------------------------------------------
# 1. Build RHS string (typos intact) + year FE ---------------
# ------------------------------------------------------------

form_baseline_relaxed <- as.formula(
  paste("VDEM_FRAGILE_BASELINE ~", paste(baseline_vars, 
                                         collapse = " + "))
)

# ------------------------------------------------------------
# 2. OLS (cluster-robust at COUNTRY_NAME) --------------------
# ------------------------------------------------------------
linear_model_relaxed <- lm(form_baseline, 
                           data = standardized_data_final)

clustered_se_linear_relaxed <- sqrt(diag(vcovCL(linear_model_relaxed, 
                                                cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 3. Logit (cluster-robust at COUNTRY_NAME) ------------------
# ------------------------------------------------------------
logit_model_relaxed  <- glm(form_baseline_relaxed,
                            family = binomial(link = "logit"),
                            data   = standardized_data_final)

clustered_se_logit_relaxed <- sqrt(diag(vcovCL(logit_model_relaxed, 
                                               cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 4. Fixed-effects (within) + Driscoll–Kraay SEs -------------
# ------------------------------------------------------------
fixed_effects_model_relaxed <- plm(form_baseline_relaxed,
                                   data  = standardized_data_final,
                                   model = "within",
                                   index = c("ISO_CODE_3", "YEAR"))

fe_dk_se_relaxed <- sqrt(diag(vcovSCC(fixed_effects_model_relaxed,
                                      maxlag = 3,
                                      type = "HC1")))

# ------------------------------------------------------------
# 5. Fit statistics ------------------------------------------
# ------------------------------------------------------------
R2_linear_relaxed <- summary(linear_model_relaxed)$r.squared

logit_null_relaxed <- glm(VDEM_FRAGILE_BASELINE ~ 1,
                          family = binomial,
                          data   = standardized_data_final)

McFadden_R2_relaxed <- 1 - (as.numeric(logLik(logit_model_relaxed)) /
                            as.numeric(logLik(logit_null_relaxed)))

add_stats_relaxed <- list(
  c("R-squared (OLS)",          sprintf("%.3f", R2_linear_relaxed),      ""),
  c("Pseudo R-squared (Logit)", "",                              sprintf("%.3f", McFadden_R2_relaxed)),
  c("AIC",                      sprintf("%.2f", AIC(linear_model_relaxed)), sprintf("%.2f", AIC(logit_model_relaxed))),
  c("BIC",                      sprintf("%.2f", BIC(linear_model_relaxed)), sprintf("%.2f", BIC(logit_model_relaxed)))
)
```

```{r, linear-binomial-relaxed, echo = FALSE, message = FALSE, warning = FALSE, results = 'asis'}
# ------------------------------------------------------------
# 6. LaTeX table (labels unchanged, year FE omitted) ---------
# ------------------------------------------------------------
stargazer(linear_model_relaxed,
          logit_model_relaxed,
          fixed_effects_model_relaxed,
          se = list(clustered_se_linear_relaxed,
                    clustered_se_logit_relaxed,
                    fe_dk_se_relaxed),
          covariate.labels = c("Combined polity score",
                               "Institutional autocracy score",
                               "GDP growth",
                               "MAX conflict intensity",
                               "N total troops involved",
                               "Territorial fragmentation",
                               "Executive recruitment",
                               "Institutional democracy score",
                               "Liberal democracy score",
                               "Regime durability in years",
                               "GDP per capita growth",
                               "Political regime",
                               "N war of fronts",
                               "ODA received per capita",
                               "Conflitc cumul. intensity across years",
                               "Electoral democracy score",
                               "Political competition score",
                               "Conflict intensity year",
                               "AVG conflict intensity",
                               "GDP deflator",
                               "Partial democracy with factionalism",
                               "Constant"),
          omit = "factor\\(YEAR\\)",     # hide year dummies
          header = FALSE,
          no.space = TRUE,
          font.size = "small",
          add.lines = add_stats_relaxed,
          title = "Linear, Logistic and FE Estimates of Fragility Status (Relaxed Classification)",
          dep.var.caption = "Fragility Status",
          dep.var.labels = "")
```

## A.2. Logit and Binomial Coefficients: Strict Classification

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ------------------------------------------------------------
# 1. Build RHS string (typos intact) + year FE ---------------
# ------------------------------------------------------------

form_baseline_strict <- as.formula(
  paste("VDEM_FRAGILE_ENDLINE ~", paste(baseline_vars, 
                                         collapse = " + "))
)

# ------------------------------------------------------------
# 2. OLS (cluster-robust at COUNTRY_NAME) --------------------
# ------------------------------------------------------------
linear_model_strict <- lm(form_baseline, 
                          data = standardized_data_final)

clustered_se_linear_strict <- sqrt(diag(vcovCL(linear_model_strict, 
                                               cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 3. Logit (cluster-robust at COUNTRY_NAME) ------------------
# ------------------------------------------------------------
logit_model_strict  <- glm(form_baseline_strict,
                           family = binomial(link = "logit"),
                           data   = standardized_data_final)

clustered_se_logit_strict <- sqrt(diag(vcovCL(logit_model_strict, 
                                              cluster = ~ COUNTRY_NAME)))

# ------------------------------------------------------------
# 4. Fixed-effects (within) + Driscoll–Kraay SEs -------------
# ------------------------------------------------------------
fixed_effects_model_strict <- plm(form_baseline_strict,
                                  data  = standardized_data_final,
                                  model = "within",
                                  index = c("ISO_CODE_3", "YEAR"))

fe_dk_se_strict <- sqrt(diag(vcovSCC(fixed_effects_model_strict,
                                     maxlag = 3,
                                     type = "HC1")))

# ------------------------------------------------------------
# 5. Fit statistics ------------------------------------------
# ------------------------------------------------------------
R2_linear_strict <- summary(linear_model_strict)$r.squared

logit_null_strict <- glm(VDEM_FRAGILE_BASELINE ~ 1,
                         family = binomial,
                         data   = standardized_data_final)

McFadden_R2_strict <- 1 - (as.numeric(logLik(logit_model_strict)) /
                           as.numeric(logLik(logit_null_strict)))

add_stats_strict <- list(
  c("R-squared (OLS)",          sprintf("%.3f", R2_linear_strict),      ""),
  c("Pseudo R-squared (Logit)", "",                              sprintf("%.3f", McFadden_R2_strict)),
  c("AIC",                      sprintf("%.2f", AIC(linear_model_strict)), sprintf("%.2f", AIC(logit_model_strict))),
  c("BIC",                      sprintf("%.2f", BIC(linear_model_strict)), sprintf("%.2f", BIC(logit_model_strict)))
)
```

```{r, linear-binomial-strict, echo = FALSE, message = FALSE, warning = FALSE, results = 'asis'}
# ------------------------------------------------------------
# 6. LaTeX table (labels unchanged, year FE omitted) ---------
# ------------------------------------------------------------
stargazer(linear_model_strict,
          logit_model_strict,
          fixed_effects_model_strict,
          se = list(clustered_se_linear_strict,
                    clustered_se_logit_strict,
                    fe_dk_se_strict),
          covariate.labels = c("Combined polity score",
                               "Institutional autocracy score",
                               "GDP growth",
                               "MAX conflict intensity",
                               "N total troops involved",
                               "Territorial fragmentation",
                               "Executive recruitment",
                               "Institutional democracy score",
                               "Liberal democracy score",
                               "Regime durability in years",
                               "GDP per capita growth",
                               "Political regime",
                               "N war of fronts",
                               "ODA received per capita",
                               "Conflitc cumul. intensity across years",
                               "Electoral democracy score",
                               "Political competition score",
                               "Conflict intensity year",
                               "AVG conflict intensity",
                               "GDP deflator",
                               "Partial democracy with factionalism",
                               "Constant"),
          omit = "factor\\(YEAR\\)",     
          header = FALSE,
          no.space = TRUE,
          font.size = "small",
          add.lines = add_stats_strict,
          title = "Linear, Logistic and FE Estimates of Fragility Status (Strict Classification)",
          dep.var.caption = "Fragility Status",
          dep.var.labels = "")
```
