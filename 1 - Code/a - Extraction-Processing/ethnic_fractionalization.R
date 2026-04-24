# Author: Cedric Antunes (FGV-CEPESP) ------------------------------------------
# Date: April, 2026 ------------------------------------------------------------

# Cleaning my environment
rm(list = ls())

# Managing memory
gc()

# Required packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(countrycode)
})

# ------------------------------------------------------------------------------
# Loading data -----------------------------------------------------------------
# ------------------------------------------------------------------------------
ef <- read.csv("D:/Users/cedric/Downloads/dataverse_files/HIEF_data.csv")

# ------------------------------------------------------------------------------
# Cleaning data ----------------------------------------------------------------
# ------------------------------------------------------------------------------
ef <- ef |>
  # Filtering for years of interest
  filter(Year >= 1970 & Year <= 2023) 

# Setting historical/ambiguous ISO codes manually
historical_matches <- c(
  "Czechoslovakia"              = "CSK",
  "German Democratic Republic"  = "DDR",
  "Republic of Vietnam"         = "VNM",
  "Yemen Arab Republic"         = "YEM",
  "Yemen PDR"                   = "YMD",
  "Yugoslavia"                  = "YUG"
)

# Adding ISO-CODE3 for matching
ef <- ef |>
  filter(Year >= 1970 & Year <= 2023) |>
  mutate(
    ISO_CODE_3 = countrycode(
      sourcevar = Country,
      origin = "country.name",
      destination = "iso3c",
      custom_match = historical_matches,
      warn = FALSE
    )
  )

# Final clean ------------------------------------------------------------------
ef <- ef |>
  select(-Country) |>
  mutate(Year = as.character(Year)) |>
  rename(YEAR = Year,
         EF_INDEX = EFindex)

# ------------------------------------------------------------------------------
# Saving the data --------------------------------------------------------------
# ------------------------------------------------------------------------------
save(ef, file = "ef_clean.RDS")
