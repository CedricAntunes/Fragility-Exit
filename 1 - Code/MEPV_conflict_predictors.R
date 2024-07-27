# Attaching my packages ---------------------------------------------------
library(here)
library(tidyverse)
library(data.table)
library(countrycode)
library(readr)
library(stringr)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

# Cleaning my environment
rm(list = ls())

# Loading data ------------------------------------------------------------
mepv2018 <- read_excel("MEPVv2018.xls") |>
  # Filtering for years of interest
  filter(year >= 1970 & year <= 2023) |>
  # Renaming variable for perfect match 
  rename("iso_code3" = "scode")

# Also adding Iso-code2 for perfect match 
mepv2018$iso_code2 <- countrycode(sourcevar = mepv2018$country, 
                                  origin = "country.name", 
                                  destination = "iso2c")

# Checking for NAs: cases in which setting the right ISO-code did not work
cbind(lapply(lapply(mepv2018, is.na), sum)) 

missing_iso <- mepv2018[is.na(mepv2018$iso_code2), ] # We have some noise here

# Selecting irrelevant vars: international conflict
irrelevant_vars <- c("intviol", "intwar", "inttot", "actotal")

# Removing irrelevant vars: keeping only domestic/civil conflict
mepv2018 <- mepv2018[, -which(names(mepv2018) %in% irrelevant_vars)]

# Saving data ------------------------------------------------------------------
save(mepv2018, file = "mepv2018.RDS")
write.csv(mepv2018, file = "mep2018.csv")
