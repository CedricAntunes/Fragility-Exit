# Attaching my packages --------------------------------------------------------
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

# Unzipping UCDP/PRIO Armed Conflict files
unzip(
  here("data", "ucdp-prio-acd-231-xlsx.zip"),
  exdir = here("data"))

# Cleaning data ----------------------------------------------------------------

# Extracting only years of interest (starting with 1970 and ending with 2018)
ucdp <- read_excel("data/UcdpPrioConflict_v23_1.xlsx")|>
  # Setting year var as numeric
  mutate(year = as.numeric(year)) |>
  # Filtering for years of interest
  filter(year >= 1970 & year <= 2023)

# Adding Isoc2c (iso-code) variable to the dataset
ucdp$iso_code2 <- countrycode(sourcevar = ucdp$location, 
                              origin = "country.name", 
                              destination = "iso2c")
  

## For sure, this is abit problematic

# Checking for NAs: cases in which setting the right ISO-code did not work
cbind(lapply(lapply(ucdp, is.na), sum)) # ISO-code did not work for ~5% of cases

# Checking those cases
missing_iso <- ucdp[is.na(ucdp$iso_code2), ] # The ´location´ consist of more
# than one country

# Let's fix this: create a string
missing_countries <- missing_iso$location

# Splitting countries: creating a single observation for 
# each individual country
split_countries <- strsplit(missing_countries, ", ")

# Fixing this headeach: create a single row for each individual country
# and repeat the data entries
missing_iso <- data.frame(
  location = rep(missing_iso$location, sapply(split_countries, length)),
  # Repeat all other variables (for each indivudla country) with lapply
  do.call(cbind, lapply(missing_iso[, -c(2)], 
                        function(x) rep(x, 
                                        sapply(split_countries, length)))),
  # Creating a single observation for each individual country
  individual_country = unlist(split_countries),
  stringsAsFactors = FALSE
)

# Assignning ISO-code to the new split observations
missing_iso$iso_code2 <- countrycode(sourcevar = missing_iso$individual_country, 
                                     origin = "country.name", 
                                     destination = "iso2c")

cbind(lapply(lapply(missing_iso, is.na), sum)) # We still have a problem for
# 23 observations


different_vars <- setdiff(names(ucdp), names(missing_iso))

missing_iso_still <- missing_iso[is.na(missing_iso$iso_code2), ]
# The problem is Yemen!

# Assigning the Yemen ISO-code
missing_iso$iso_code2[is.na(missing_iso$iso_code2)] <- "YE"

# We are finally ready to go!
cbind(lapply(lapply(missing_iso, is.na), sum))

ucdp_clean <- ucdp |> 
  filter(!is.na(iso_code2)) |>
  mutate(individual_country = location)

# Final df 
ucdp_clean <- rbind(ucdp, missing_iso)

## Adding iso-code 3 to the final df
ucdp_clean$iso_code3 <- countrycode(sourcevar = ucdp_clean$individual_country, 
                                    origin = "country.name", 
                                    destination = "iso3c")

ucdp_clean$iso_code3[is.na(ucdp_clean$iso_code3)] <- "YEM"


# Countries that appear more than once -----------------------------------------

country_year_counts <- ucdp_clean |>
  # Grouping by country-year
  group_by(individual_country, year) |>
  # Summarizing frequency
  summarise(count = n())

# Filtering countries engaged in more than one war/front in a single year
duplicates <- country_year_counts |>
  filter(count > 1)

# Creating new df only with countries engaged in more than one war/front in a 
# single year
more_than_one_front_year_ucdp <- ucdp_clean %>%
                                 semi_join(duplicates, 
                                 by = c("individual_country", "year"))

ucdp_clean <- ucdp_clean |>
  # Merging the count of duplicates back into the original df
  left_join(country_year_counts, 
            by = c("individual_country", "year"))

ucdp_clean <- ucdp_clean |>
  group_by(individual_country, year) |>
  mutate(MAX_CONFLICT_INTENSITY = max(intensity_level),
         MIN_CONFLICT_INTENSITY = min(intensity_level),
         AVG_CONFLICT_INTENSITY = mean(as.numeric(intensity_level)),
         N_STATES_SUPPORT_SIDE_A_WITH_TROOPS = str_count(side_a_2nd, ",") + 1,
         N_STATES_SUPPORT_SIDE_B_WITH_TROOPS = str_count(side_b_2nd, ",") + 1) |>
  # Ungroup
  ungroup() |>
  # Filtering only for (1) extrasystemic conflicts; (3) intrastate conflicts;
  # and (4) internationalized intrastate conflicts
  filter(type_of_conflict %in% c("1", "3", "4"))

# Generating the final df
ucdp_final_clean <- ucdp_clean |>
  distinct(individual_country, year,
           .keep_all = TRUE)

# Saving data ------------------------------------------------------------------

# Saving files
save(ucdp_final_clean, file = "ucdp_final_clean.RDS")
save(more_than_one_front_year_ucdp, file = "more_than_one_front_year_ucdp.RDS")
write.csv(ucdp_final_clean, "ucdp_final_clean.csv")
write.csv(more_than_one_front_year_ucdp, "more_than_one_front_year_ucdp.csv")
