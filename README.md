# Fragility Exit

> *Colonies do not cease to be colonies because they are indendependent.* <br>
> **Benjamin Disraeli**, British Prime Minister, circa 1880.
> 
> *Progress is not an illusion, it happens, but it is slow and ivariably disappointing.* <br>
> **George Orwell**, British political essayst, circa 1950.

## Working Abstract

# Project Overview 

## Repository Organization
This repository is (currently) organized accross X folders: 
* `1 - Data`: this folder is, by it's turn, organized across two sub-folders:
   + `a - Data Processing-Extraction`: which systematizes data extraction for each class of predictors from different datasources and prepare data for merging;
   + `b - Data-Merge`: sets the key for merging datasets, and clean and prepare final data for analysis.
* `2 - Functions`: set of auxiliary functions used in data visualization and analysis.

## Code Source 
### Data Extraction and Cleaning 
|Folder|Step|File Name|Predictor Class|File Description|
|---|---|---|---|---|
|World Bank data|Data extraction and cleaning|WB_economic_predictor.R|Economic|**(1)** Collecting data via API; <br> **(2)** selecting relevant indicators; <br> **(3)** pooling country-year data.|
