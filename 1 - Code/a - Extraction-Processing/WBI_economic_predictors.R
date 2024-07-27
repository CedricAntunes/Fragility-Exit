# Attaching my packages ---------------------------------------------------
library(here)
library(readr)
library(stringr)
library(dplyr)
library(WDI) # World Bank API
library(tidyr)
library(purrr)

# Creating local direcotry for storing data
dir.create("data")

rm(list = ls())

# Generating WB data -----------------------------------------------------------

# ODA indicators
oda <- WDIsearch("ODA received",
                 short = FALSE)

# ODA data 
oda_data <- WDI(
  country = "all",
  indicator = c("DT.ODA.DACD.CD", "DT.ODA.DACD.CD.PC", "DT.ODA.DACD.KD",
                "DT.ODA.MULT.CD", "DT.ODA.MULT.CD.PC", "DT.ODA.NDAC.CD",
                "DT.ODA.ODAT.PC.ZS"),
  start = 1970,
  end = 2023,
  extra = FALSE,
  cache = NULL,
  latest = NULL,
  language = "en"
)

# GDP indicators
GDP <- WDIsearch("GDP",
                 short = FALSE)

# GDP data
gdp_data <- WDI(
  country = "all",
  indicator = c("NY.GDP.PCAP.KD.ZG", "NY.GDP.MKTP.KD.ZG", "NY.GDP.DEFL.KD.ZG"),
  start = 1970,
  end = 2023,
  extra = FALSE,
  cache = NULL,
  latest = NULL,
  language = "en"
)

# Function to generate simple vector of WB indicators
cat("c(", paste0("\"", per_capita_vars, "\"", collapse = ", "), ")\n")

# Preparing and saving data ----------------------------------------------------

# ODA data
oda_data <- oda_data |>
  # Renaming vars for merging 
  rename(iso_code2 = iso2c,
         iso_code3 = iso3c)

# Saving ODA data
save(oda_data, file = "wb_oda_data.RDS")
write.csv(oda_data, "wb_oda_data.csv")

# ODA data
gdp_data <- gdp_data |>
  # Renaming vars for merging 
  rename(iso_code2 = iso2c,
         iso_code3 = iso3c)

# GDP data 
save(gdp_data, file = "gdp_data.RDS")
write.csv(gdp_data, file = "gdp_data.csv")
write.csv(GDP, "gdp_indicators.csv")
