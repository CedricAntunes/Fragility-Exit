# AUTHOR: Cedric Antunes -------------------------------------------------------
# Date: February, 2026 ---------------------------------------------------------

# Cleaning my environment 
rm(list = ls())

# Managing memory
gc()

# Required packages ------------------------------------------------------------
library(readxl)
library(janitor)
library(countrycode)

# ------------------------------------------------------------------------------
# Loading the data -------------------------------------------------------------
# ------------------------------------------------------------------------------
# Country-year loading factors -------------------------------------------------
wgi_lf <- read_excel("D:/Users/cedric/Downloads/wgi - PCA data.xlsx") |>
  clean_names() |>
  select(country_name,
         year,
         loading_factor_1_normalized)
  
# Percentiles ------------------------------------------------------------------
wgi_pc <- read_excel("D:/Users/cedric/Downloads/wgi - percentiles.xlsx") |>
  select(Year, 
         P30,
         P35) |>
  rename_with(tolower)

# ------------------------------------------------------------------------------
# Cleaning the data ------------------------------------------------------------
# ------------------------------------------------------------------------------
# Setting iso-codes 3 ----------------------------------------------------------
wgi_lf$iso_code3 <- countrycode(sourcevar = wgi_lf$country_name, 
                                origin = "country.name", 
                                destination = "iso3c")

wgi_lf$iso_code3[wgi_lf$country_name == "Kosovo"] <- "XKX"

wgi_lf$iso_code3[wgi_lf$country_name == "Netherlands Antilles (former)"] <- "ANT"

# ------------------------------------------------------------------------------
# Join -------------------------------------------------------------------------
# ------------------------------------------------------------------------------
wgi_revised <- wgi_lf |>
  left_join(wgi_pc,
            by = "year")

# Final cleaning ---------------------------------------------------------------
wgi_revised_clean <- wgi_revised |>
  mutate(year = as.character(year)) |>
  rename(ISO_CODE_3 = iso_code3,
         YEAR = year,
         WGI_30TH_PERCENTILE = p30,
         WGI_35TH_PERCENTILE = p35,
         WGI_LOADING_FACTOR_1_NORMALIZED = loading_factor_1_normalized) |>
  select(-country_name)

# Saving revised data ----------------------------------------------------------
save(wgi_revised_clean,
     file = "wgi_revised.RDS")
