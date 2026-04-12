# AUTHOR: Cedric Antunes (FGV) -------------------------------------------------
# DATE: March, 2026 ------------------------------------------------------------

# Cleaning my enrvironment 
rm(list = ls())

# Managing memory
gc()

# Required packages ------------------------------------------------------------
library(dplyr)

# Loading data -----------------------------------------------------------------
load("D:/Users/cedric/Downloads/episodes.rda")

# Cleaning & preparing data ----------------------------------------------------
ert <- episodes |>
  select(country_text_id,
         year, 
         reg_trans,
         aut_ep) |>
  filter(year %in% 1971:2022) |>
  rename(ISO_CODE_3 = country_text_id,
         YEAR = year, 
         ERT_TRANSITION = reg_trans,
         AUTOCRATIZATION_EPISODE = aut_ep) |>
  mutate(ERT_DEMOCRATIZATION = if_else(ERT_TRANSITION == 1, 1, 0),
         ERT_DEM_BREAKDOWN = if_else(ERT_TRANSITION == -1, 1, 0)) 

# Setting autocratization episodes ---------------------------------------------
ert <- ert |>
  arrange(ISO_CODE_3, YEAR) |>
  group_by(ISO_CODE_3) |>
  mutate(
    IN_AUTOCRATIZATION = AUTOCRATIZATION_EPISODE != 0 & !is.na(AUTOCRATIZATION_EPISODE),
    ERT_AUTOCRATIZATION = as.integer(IN_AUTOCRATIZATION & c(TRUE, !head(IN_AUTOCRATIZATION, -1)))
  ) |>
  ungroup() |>
  select(-IN_AUTOCRATIZATION) |>
  mutate(YEAR = as.character(YEAR))

# Saving the data --------------------------------------------------------------
save(ert,
     file = "ert_transitions.RDS")
