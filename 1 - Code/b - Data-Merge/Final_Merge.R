# Attaching my packages --------------------------------------------------------
library(here)
library(tidyverse)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)

# Cleaning my environment
rm(list = ls())

# Loading data -----------------------------------------------------------------

# Input Vars
load("orgviolence_v23.RDS")
load("wb_oda_data.RDS")
load("polV.RDS")
load("gdp_data.RDS")
load("ucdp_final_clean.RDS")
load("vdem_data.RDS")

# Output vars
load("output-data/vdem_output_final.RDS")
load("output-data/vdem_output_percentiles_final.RDS")
load("output-data/wgi_output_final.RDS")
load("output-data/wgi_output_percentiles_final.RDS")

# Merging data -----------------------------------------------------------------

# Merging World Bank data 
wb_data <- merge(oda_data,
                 gdp_data,
                 by = c("iso_code2", "year"),
                 suffixes = c("", "_gdp"),
                 no.dups = TRUE) |>
  select(-iso_code3_gdp)

# Checking for NAs
cbind(lapply(lapply(wb_vdem, is.na), sum))

# Merging World Bank data with V-DEM
wb_vdem <- merge(wb_data,
                 vdem_data,
                 by = c("iso_code3", "year"),
                 all.x = TRUE,
                 no.dups = TRUE)

# Iso-codes to drop
regions_code <- c("1A", # Arab World 
                  "1W", # World
                  "4E", # East Asia and Pacific
                  "7E", # Europe and Central Asia,
                  "8S", # South Asia
                  "B8", # Central Europe and Baltics
                  "EU", # European Union
                  "F1", # Fragile and conflict affected situations
                  "OE", # OECD members
                  "S1", # Small states
                  "S2", # Pacific island small states
                  "S3", # Caribbean small states
                  "S4", # Other small states
                  "T2", # Latin America & the Caribbean (IDA & IBRD countries)
                  "T3", # Middle East & North Africa (IDA & IBRD countries
                  "T4", # East Asia & Pacific (IDA & IBRD countries
                  "T5", # South Asia (IDA & IBRD)
                  "T6", # Sub-Saharan Africa (IDA & IBRD countries)
                  "T7", # Europe & Central Asia (IDA & IBRD countries)
                  "V1", # Pre-demographic dividend
                  "V2", # Early-demographic dividend
                  "V3", # Late-demographic dividend
                  "V4", # Post-demographic dividend
                  "XC", # Euro area
                  "XD", # High income
                  "XE", # Heavily indebted poor countries (HIPC)
                  "XF", # IBRD only
                  "XG", # IDA total
                  "XH", # IDA blend
                  "XI", # IDA only
                  "XJ", # Latin America & Caribbean (excluding high income)
                  "XL", # Least developed countries: UN classification
                  "XM", # Low income
                  "XN", # Lower middle income
                  "XO", # Low & middle income
                  "XP", # Middle income
                  "XQ", # Middle East & North Africa (excluding high income)
                  "XT", # Upper middle income
                  "XU", # North America
                  "XY", # Not classified
                  "Z4", # East Asia & Pacific
                  "Z7", # Europe & Central Asia
                  "ZF", # Sub-Saharan Africa (excluding high income)
                  "ZG", # Sub-Saharan Africa
                  "ZH", # Africa Eastern and Southern
                  "ZI", # Africa Western and Central
                  "ZJ", # Latin America & Caribbean
                  "ZQ", # Middle East & North Africa
                  "ZT" # IDA & IBRD total
)

wb_vdem_clean <- wb_vdem |>
  # Dropping all regions; keeping only countries
  filter(!iso_code2 %in% regions_code) |>
  # Keeping variables type consistent
  mutate(year = as.character(year))

# Merging Polity-V data
wb_vdem_pv <- merge(wb_vdem_clean,
                    p5v2018_final,
                    by = c("iso_code3", "year"),
                    all.x = TRUE,
                    no.dups = TRUE)

# Merging UCDP data
wb_vdem_pv_ucdp <- merge(wb_vdem_pv,
                         ucdp_final_clean,
                         by = c("iso_code3", "year"),
                         all.x = TRUE,
                         no.dups = TRUE) |>
  # Dropping irrelevant vars
  select(-iso_code2.x, -iso_code2.y, -country_gdp, -country_name, 
         -country.y) |>
  # Keeping labeling consistency
  rename(country_name = country.x)

clean_data <- wb_vdem_pv_ucdp |>
  # Dropping irrelevant variables
  select(-iso_code2, -p5, -cyear, -ccode, -scode,
         -side_a_id, -side_b_id, -gwno_a, -gwno_a_2nd,
         -gwno_b, -gwno_b_2nd, -gwno_loc, -region, -version,
         -individual_country, -exconst) |>
  # Renaming variables 
  rename(
    # World Bank indicators
    ISO_CODE_3 = iso_code3,
    YEAR = year, 
    COUNTRY_NAME = country_name, 
    ODA_RECEIVED_PER_CAPITA = DT.ODA.ODAT.PC.ZS,
    GDP_PER_CAPITA_GROWTH = NY.GDP.PCAP.KD.ZG,
    GDP_GROWTH = NY.GDP.MKTP.KD.ZG,
    GDP_DEFLATOR = NY.GDP.DEFL.KD.ZG,
    # V-Dem indicators
    POLITICAL_REGIME = v2x_regime,
    ELECTORAL_DEMOCRACY_SCORE = v2x_polyarchy,
    LIBERAL_DEMOCRACY_SCORE = v2x_liberal,
    # Polity-V confidence index
    POLITYV_CONFIDENCE_INDEX = flag,
    TERRITORIAL_FRAGMENTATION = fragment,
    INSTITUTIONAL_DEMOCRACY_SOCRE = democ,
    INSTITUTIONAL_AUTOCRACY_SCORE = autoc, 
    COMBINED_POLITY_SCORE = polity,
    REVISED_COMBINED_POLITY_SCORE = polity2,
    REGIME_DURABILITY_YEARS = durable,
    INSTITUTIONAL_EXECUTIVE_RECRUTIMENT = xrreg,
    COMPETITIVENESS_EXECUTIVE_RECRUITMENT = xrcomp,
    OPENESS_EXECUTIVE_RECRUTIMENT = xropen,
    EXECUTIVE_CONSTRAINT_SCORE = xconst,
    INSTITUTIONAL_PARTICIPATION = parreg,
    COMPETITIVENESS_PARTICIPATION = parcomp,
    PRIOR_COMBINED_POLITY_SCORE = prior,
    POLITY_END_MONTH = emonth,
    POLITY_END_DAY = eday,
    POLITY_END_YEAR = eyear,
    END_DATE_PRECISION = eprec,
    INTERIM_POLITY_SCORE = interim,
    POLITY_BEGIN_MONTH = bmonth,
    POLITY_BEGIN_DAY = bday,
    POLITY_BEGIN_YEAR = byear,
    BEGIN_DATE_PRECISION = bprec,
    POST_POLITY_SCORE = post,
    TOTAL_CHANGE_POLITY_SCORE = change,
    REGIME_TRANSITION_COMPLETED = d5,
    STATE_FAILURE = sf,
    REGIME_TRANSITION = regtrans,
    EXECUTIVE_RECRUTIMENT_SCORE = exrec,
    POLITICAL_COMPETITION_SCORE = polcomp,
    # UCDP data 
    CONFLICT_ID = conflict_id,
    CONFLICT_LOCATION = location,
    COUNTRY_CONFLICT_SIDE_A = side_a,
    COUNTRY_OR_ACTOR_CONFLICT_SIDE_B = side_b,
    CONFLICT_CAUSE = incompatibility,
    TERRITORY_UNDER_DISPUTE = territory_name,
    CONFLICT_INTENSITY_YEAR = intensity_level,
    CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS = cumulative_intensity,
    CONFLICT_TYPE = type_of_conflict,
    CONFLICT_START_DATE = start_date,
    PRECISION_OF_CONFLICT_START_DATE = start_prec,
    CONFLICT_DEATHS_THRESHOLD_DATE = start_date2,
    PRECISION_OF_CONFLICT_DEATHS_THRESHOLD_DATE = start_prec2,
    CONFLICT_INACTIVE = ep_end,
    CONFLICT_END_DATE = ep_end_date,
    PRECISION_OF_CONFLICT_END_DATE = ep_end_prec,
    N_WAR_FRONTS = count,
    STATES_SUPPORT_SIDE_A_WITH_TROOPS = side_a_2nd,
    STATES_SUPPORT_SIDE_B_WITH_TROOPS = side_b_2nd) |>
  # Preparing Vars for analysis
  mutate(across(c(N_STATES_SUPPORT_SIDE_A_WITH_TROOPS, 
                  N_STATES_SUPPORT_SIDE_B_WITH_TROOPS, 
                  CONFLICT_INTENSITY_YEAR, 
                  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS,
                  N_WAR_FRONTS,
                  MAX_CONFLICT_INTENSITY,
                  MIN_CONFLICT_INTENSITY,
                  AVG_CONFLICT_INTENSITY), ~ replace(., is.na(.), as.numeric(0))),
         CONFLICT_INTENSITY_YEAR = as.numeric(CONFLICT_INTENSITY_YEAR),
         CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS = as.numeric(CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS),
         N_WAR_FRONTS = as.numeric(N_WAR_FRONTS),
         MAX_CONFLICT_INTENSITY = as.numeric(MAX_CONFLICT_INTENSITY),
         MIN_CONFLICT_INTENSITY = as.numeric(MAX_CONFLICT_INTENSITY),
         YEAR = as.numeric(YEAR),
         N_TOTAL_TROOPS = N_STATES_SUPPORT_SIDE_A_WITH_TROOPS + N_STATES_SUPPORT_SIDE_B_WITH_TROOPS)

          
# Checking for NAs
cbind(lapply(lapply(clean_data, is.na), sum))

# Keeping only UN members ------------------------------------------------------
un_countries <- c(
  "AFG", "ALB", "DZA", "AND", "AGO", "ATG", "ARG", "ARM", "AUS", "AUT", "AZE", "BHS",
  "BHR", "BGD", "BRB", "BLR", "BEL", "BLZ", "BEN", "BTN", "BOL", "BIH", "BWA", "BRA",
  "BRN", "BGR", "BFA", "BDI", "CPV", "KHM", "CMR", "CAN", "CAF", "TCD", "CHL", "CHN",
  "COL", "COM", "COG", "COK", "CRI", "CIV", "HRV", "CUB", "CYP", "CZE", "COD", "DNK",
  "DJI", "DMA", "DOM", "ECU", "EGY", "SLV", "GNQ", "ERI", "EST", "SWZ", "ETH", "FJI",
  "FIN", "FRA", "GAB", "GMB", "GEO", "DEU", "GHA", "GRC", "GRD", "GTM", "GIN", "GNB",
  "GUY", "HTI", "HND", "HUN", "ISL", "IND", "IDN", "IRN", "IRQ", "IRL", "ISR", "ITA",
  "JAM", "JPN", "JOR", "KAZ", "KEN", "KIR", "KWT", "KGZ", "LAO", "LVA", "LBN", "LSO",
  "LBR", "LBY", "LIE", "LTU", "LUX", "MDG", "MWI", "MYS", "MDV", "MLI", "MLT", "MHL",
  "MRT", "MUS", "MEX", "FSM", "MDA", "MCO", "MNG", "MNE", "MAR", "MOZ", "MMR", "NAM",
  "NRU", "NPL", "NLD", "NZL", "NIC", "NER", "NGA", "NIU", "PRK", "MKD", "NOR", "OMN",
  "PAK", "PLW", "PAN", "PNG", "PRY", "PER", "PHL", "POL", "PRT", "QAT", "ROU", "RUS",
  "RWA", "KNA", "LCA", "VCT", "WSM", "SMR", "STP", "SAU", "SEN", "SRB", "SYC", "SLE",
  "SGP", "SVK", "SVN", "SLB", "SOM", "ZAF", "KOR", "SSD", "ESP", "LKA", "SDN", "SUR",
  "SWE", "CHE", "SYR", "TWN", "TJK", "TZA", "THA", "TLS", "TGO", "TON", "TTO", "TUN",
  "TUR", "TKM", "TUV", "UGA", "UKR", "ARE", "GBR", "USA", "URY", "UZB", "VUT", "VEN",
  "VNM", "YEM", "ZMB", "ZWE"
)

clean_data_un_only <- clean_data |>
  filter(ISO_CODE_3 %in% un_countries)


# Imputing Success Cases Data --------------------------------------------------

# Merging V-Dem output data
final_clean_data <- merge(clean_data,
                          wgi_normalized_final,
                          by = c("ISO_CODE_3", "YEAR"),
                          all.x = TRUE,
                          no.dups = TRUE) |>
  distinct(ISO_CODE_3, 
           YEAR, 
           .keep_all = TRUE)

# Always double check for duplicates whenever merging 
duplicated_rows <- final_clean_data[duplicated(final_clean_data[c("ISO_CODE_3", "YEAR")]) | duplicated(final_clean_data[c("ISO_CODE_3", "YEAR")], fromLast = TRUE), ]

# Merging WGI output data
final_clean_data <- merge(final_clean_data,
                          vdem_normalized_final,
                          by = c("ISO_CODE_3", "YEAR"),
                          all.x = TRUE,
                          no.dups = TRUE) |>
  distinct(ISO_CODE_3, 
           YEAR, 
           .keep_all = TRUE)

# Always double check for duplicates whenever merging 
duplicated_rows <- final_clean_data[duplicated(final_clean_data[c("ISO_CODE_3", "YEAR")]) | duplicated(final_clean_data[c("ISO_CODE_3", "YEAR")], fromLast = TRUE), ]

# Matching V-Dem percentiles
final_clean_percentiles_data <- merge(final_clean_data, 
                                      vdem_percentiles, 
                                      by = "YEAR", 
                                      all.x = TRUE)

# Matching WGI percentiles
final_clean_percentiles_data <- merge(final_clean_percentiles_data, 
                                      wgi_percentiles, 
                                      by = "YEAR", 
                                      all.x = TRUE) |>
  # Dropping redundant vars due to matching
  select(-COUNTRY_NAME.y, -COUNTRY_NAME) |>
  rename(COUNTRY_NAME = COUNTRY_NAME.x)

# Saving data ------------------------------------------------------------------
save(clean_data, file = "final_data.RDS")
save(final_clean_percentiles_data, file = "final_clean_output_data.RDS")
save(clean_data_un_only, final = "final_un_only_data.RDS")
write.csv(clean_data, "final_data.csv")
write.csv(final_clean_percentiles_data, "final_clean_percentiles_data.csv")
write.csv(clean_data_un_only, "final_un_only_data.csv")
