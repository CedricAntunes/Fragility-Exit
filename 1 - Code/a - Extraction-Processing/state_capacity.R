# AUTHOR: Cedric Antunes (FGV) -------------------------------------------------
# DATE: April, 2026 ------------------------------------------------------------

# Cleaning my environment
rm(list = ls())

# Managing memory
gc()

# Required packages ------------------------------------------------------------
library(dplyr)
library(haven)

# Loading and cleaning data ----------------------------------------------------
sc <- read_dta("D:/Users/cedric/Downloads/StateCapacityDataset_v1.dta/StateCapacityDataset_v1.dta") |>
  select(year, 
         iso3,
         Capacity,
         Capacity_sd,
         StateHist50s,
         taxrev_gdp) |>
  rename(YEAR = year,
         ISO_CODE_3 = iso3,
         STATE_CAPACITY = Capacity,
         STATE_CAPACITY_SD = Capacity_sd,
         TAX_REVENUE_GDP = taxrev_gdp) |>
  mutate(YEAR = as.character(YEAR))

# Saving the data --------------------------------------------------------------
save(sc,
     file = "state_capacity.RDS")
