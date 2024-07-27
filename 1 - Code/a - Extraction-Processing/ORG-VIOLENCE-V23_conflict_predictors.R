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
orgviolence_v23 <- read_excel("data/organizedviolencecy_v23_1.xlsx") |>
  rename(year = year_cy)

# Adding iso-code2 for matching 
orgviolence_v23$iso_code2 <- countrycode(sourcevar = orgviolence_v23$country_cy, 
                                         origin = "country.name", 
                                         destination = "iso2c")

# Adding iso-code3 for matching 
orgviolence_v23$iso_code3 <- countrycode(sourcevar = orgviolence_v23$country_cy, 
                                         origin = "country.name", 
                                         destination = "iso3c")

# Checking for NAs: cases in which setting the ISO-code did not work
cbind(lapply(lapply(orgviolence_v23, is.na), sum)) # Missing ~0.8% of cases

# Some noise again
missing_iso <- orgviolence_v23[is.na(orgviolence_v23$iso_code2), ] 

# Assigning ISO-CODE3 to the problematic cases 
replacement_iso3 <- c("Czechoslovakia" = "CZE", 
                      "German Democratic Republic" = "GDR", 
                      "Kosovo" = "KOS", 
                      "South Vietnam" = "RVN",
                      "Yemen (North Yemen)" = "YAR", 
                      "Yemen (South Yemen)" = "YPR",
                      "Yugoslavia" = "YUG")

# Preparing the final dataset 
orgviolence_v23$iso_code3 <- ifelse(is.na(orgviolence_v23$iso_code3),
                                    # Condition for replacement 
                                    replacement_iso3[orgviolence_v23$country_cy], 
                                    orgviolence_v23$iso_code3)

# Checking for NAs: we are ready to go!
cbind(lapply(lapply(orgviolence_v23, is.na), sum))

# Saving data -------------------------------------------------------------
save(orgviolence_v23, file = "orgviolence_v23.RDS")
write.csv(orgviolence_v23, file = "orgviolence_v23.csv")
