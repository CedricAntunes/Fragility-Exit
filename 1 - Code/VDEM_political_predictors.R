# Attaching my packages ---------------------------------------------------
library(devtools)
library(countrycode)

# Setting V-DEM API
devtools::install_github("vdeminstitute/vdemdata")

# Loading V-DEM data
library(vdemdata)

# Cleaning my environment
rm(list = ls())

# Generating data ---------------------------------------------------------

# V-DEM data 
vdem_data <- vdem |>
  # Compiling data frame
  as.data.frame() |>
  # Filtering for years of interest
  filter(year >= 1970 & year <= 2023) |>
  # Selecting only variables of interest
  select(country_name, country_text_id, year, v2x_regime, 
         v2x_polyarchy, v2x_liberal) |>
  rename(iso_code3 = country_text_id)

# Adding iso-code2 for matching 
vdem_data$iso_code2 <- countrycode(sourcevar = vdem$country_name, 
                                   origin = "country.name", 
                                   destination = "iso2c")

# Checking for NAs: cases in which ISO-code did not work 
cbind(lapply(lapply(vdem_data, is.na), sum)) # Missing 1.7% of cases

# Noisy observations 
missing_iso <- vdem[is.na(vdem_data$iso_code2), ]

# Saving data ------------------------------------------------------------------

save(vdem_data, file = "vdem_data.RDS")
write.csv(vdem_data, "vdem_data.csv")
