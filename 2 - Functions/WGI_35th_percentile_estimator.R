# Function to fill missing values for specific WGI years: 1997, 1999, and 2001
estimating_WGI_35th_missing_values <- function(final_clean_percentiles_data) {
  # Looping through each country
  for (country in unique(final_clean_percentiles_data$COUNTRY_NAME)) {
    # Subsetting the data for each country
    countries <- final_clean_percentiles_data[final_clean_percentiles_data$COUNTRY_NAME == country, ]
    
    # Filling missing values for 1997
    if (any(is.na(final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 1997]))) {
      value_1996 <- final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 1996]
      value_1998 <- final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 1998]
      final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$COUNTRY_NAME == country & final_clean_percentiles_data$YEAR == 1997] <- (value_1996 + value_1998) / 2
    }
    
    # Filling missing values for 1999
    if (any(is.na(final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 1999]))) {
      value_1998 <- final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 1998]
      value_2000 <- final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 2000]
      final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$COUNTRY_NAME == country & final_clean_percentiles_data$YEAR == 1999] <- (value_1998 + value_2000) / 2
    }
  }
  
  # Filling missing values for 2001 after filling 1999
  for (country in unique(final_clean_percentiles_data$COUNTRY_NAME)) {
    # Subsetting the data frame for the current country
    countries <- final_clean_percentiles_data[final_clean_percentiles_data$COUNTRY_NAME == country, ]
    
    # Filling missing values for 2001
    if (any(is.na(final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 2001]))) {
      value_2000 <- final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 2000]
      value_2002 <- final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$YEAR == 2002]
      final_clean_percentiles_data$WGI_35TH_PERCENTILE[final_clean_percentiles_data$COUNTRY_NAME == country & final_clean_percentiles_data$YEAR == 2001] <- (value_2000 + value_2002) / 2
    }
  }
  
  return(final_clean_percentiles_data)
}
