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

# Loading data -----------------------------------------------------------------
p5v2018 <- read_excel("p5v2018.xls") |>
  # Filtering for years of interest
  filter(year >= 1970 & year <= 2023)

p5v2018$iso_code2 <- countrycode(sourcevar = p5v2018$country, 
                                 origin = "country.name", 
                                 destination = "iso2c")

p5v2018$iso_code3 <- countrycode(sourcevar = p5v2018$country, 
                                 origin = "country.name", 
                                 destination = "iso3c")

# Checking for NAs: cases in which setting the ISO-code did not work
# ISO-code did not work for ~0.8% of cases
cbind(lapply(lapply(p5v2018, is.na), sum))

# Some noise again
missing_iso <- p5v2018[is.na(p5v2018$iso_code2), ] 

glimpse(p5v2018)

# Goldstone et al. (2010) ------------------------------------------------------
p5v2018_final <- p5v2018 |>
  # Creating a new (dichotomous) var for partial democracy with factionalism
  mutate(PARTIAL_DEMOCRACY_WITH_FACTIONALISM = 
           ifelse(exrec %in% c(6:8) & parcomp == 3, 1, 0))

# Saving data ------------------------------------------------------------------
save(p5v2018_final, file = "polV.RDS")
write.csv(p5v2018_final, "polv.csv")
