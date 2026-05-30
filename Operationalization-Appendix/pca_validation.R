# ==============================================================================
# Pairwise Correlations Among Fragility Components
# V-Dem and WGI versions
# ==============================================================================
# Author: Cedric Antunes (FGV-CEPESP) ------------------------------------------
# Date: May, 2026 --------------------------------------------------------------

# Assumes standardized_data_final is loaded

# Required packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(knitr)
  library(jsonlite)
  library(httr2)
  library(vdemdata)
  library(WDI)
})

# ------------------------------------------------------------------------------
# Loading V-Dem data -----------------------------------------------------------
# ------------------------------------------------------------------------------
vdem_raw <- vdemdata::vdem

# ------------------------------------------------------------------------------
# Required variables -----------------------------------------------------------
# ------------------------------------------------------------------------------
required_main_vars <- c(
  "COUNTRY_NAME",
  "ISO_CODE_3",
  "YEAR",
  "VDEM_FRAGILE_IDEAL",
  "WGI_FRAGILE_IDEAL"
)

required_vdem_vars <- c(
  "country_text_id",
  "country_name",
  "year",
  "v2x_corr",
  "v2caviol",
  "v2x_rule",
  "v2clrspct"
)

# ------------------------------------------------------------------------------
# Helpers ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
rescale01 <- function(x) {
  x_nonmiss <- x[!is.na(x)]
  
  if (length(x_nonmiss) == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  rng <- range(x_nonmiss)
  
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  (x - rng[1]) / diff(rng)
}

safe_cor <- function(x, y) {
  cc <- complete.cases(x, y)
  
  if (sum(cc) < 10) {
    return(NA_real_)
  }
  
  suppressWarnings(
    cor(x[cc], y[cc], use = "complete.obs")
  )
}

orient_to_fragility <- function(x, fragile_binary) {
  r <- safe_cor(x, fragile_binary)
  
  if (is.na(r)) return(x)
  
  # If negatively correlated with the fragile dummy,
  # flip so higher values mean more fragility.
  if (r < 0) -x else x
}

make_correlation_matrix <- function(data, component_cols) {
  data |>
    select(all_of(component_cols)) |>
    cor(
      use = "pairwise.complete.obs",
      method = "pearson"
    )
}

lower_triangle_for_latex <- function(cor_mat, panel_label) {
  cor_df <- as.data.frame(round(cor_mat, 3)) |>
    rownames_to_column("Component")
  
  component_names <- cor_df$Component
  
  for (j in seq_along(component_names)) {
    for (i in seq_len(nrow(cor_df))) {
      if (j > i) {
        cor_df[i, j + 1] <- NA
      }
    }
  }
  
  cor_df |>
    mutate(
      panel = panel_label,
      across(
        -c(Component, panel),
        ~ if_else(is.na(.x), "", sprintf("%.3f", as.numeric(.x)))
      )
    )
}

# ------------------------------------------------------------------------------
# Loading WGI data -------------------------------------------------------------
# ------------------------------------------------------------------------------
start_year <- min(standardized_data_final$YEAR, na.rm = TRUE)
end_year   <- max(standardized_data_final$YEAR, na.rm = TRUE)

# WGI covers 1996 onward.
start_year_wgi <- max(start_year, 1996)
end_year_wgi   <- min(end_year, 2024)

wgi_indicators <- c(
  wgi_control_corruption       = "CC.EST",
  wgi_political_stability      = "PV.EST",
  wgi_rule_of_law              = "RL.EST",
  wgi_government_effectiveness = "GE.EST"
)

try_wdi_source <- function(source_id) {
  tryCatch(
    {
      WDI::WDI(
        country = "all",
        indicator = wgi_indicators,
        start = start_year_wgi,
        end = end_year_wgi,
        source = source_id,
        extra = TRUE
      )
    },
    error = function(e) {
      message("WDI::WDI(source = ", source_id, ") failed.")
      message("Original error: ", conditionMessage(e))
      NULL
    }
  )
}

wgi_raw <- try_wdi_source(3)

if (is.null(wgi_raw) || nrow(wgi_raw) == 0) {
  wgi_raw <- try_wdi_source(75)
}

fetch_wb_indicator <- function(indicator_code, indicator_name, source_id,
                               start_year, end_year) {
  
  fetch_page <- function(page) {
    url <- paste0(
      "https://api.worldbank.org/v2/country/all/indicator/",
      indicator_code
    )
    
    resp <- request(url) |>
      req_url_query(
        format = "json",
        source = source_id,
        date = paste0(start_year, ":", end_year),
        per_page = 20000,
        page = page
      ) |>
      req_perform()
    
    parsed <- fromJSON(
      resp_body_string(resp),
      flatten = TRUE
    )
    
    if (length(parsed) < 2 || is.null(parsed[[2]])) {
      return(tibble())
    }
    
    as_tibble(parsed[[2]]) |>
      transmute(
        iso3c = as.character(countryiso3code),
        country = as.character(country.value),
        year = as.integer(date),
        value = as.numeric(value)
      ) |>
      filter(
        !is.na(iso3c),
        iso3c != "",
        !is.na(year)
      )
  }
  
  url <- paste0(
    "https://api.worldbank.org/v2/country/all/indicator/",
    indicator_code
  )
  
  first_resp <- request(url) |>
    req_url_query(
      format = "json",
      source = source_id,
      date = paste0(start_year, ":", end_year),
      per_page = 20000,
      page = 1
    ) |>
    req_perform()
  
  first_parsed <- fromJSON(
    resp_body_string(first_resp),
    flatten = TRUE
  )
  
  if (length(first_parsed) < 2 || is.null(first_parsed[[2]])) {
    warning(
      "No data returned for indicator ",
      indicator_code,
      " from source ",
      source_id
    )
    
    return(tibble(
      iso3c = character(),
      country = character(),
      year = integer(),
      !!indicator_name := numeric()
    ))
  }
  
  n_pages <- as.integer(first_parsed[[1]]$pages[1])
  
  map_dfr(seq_len(n_pages), fetch_page) |>
    rename(!!indicator_name := value)
}

fetch_wgi_api_source <- function(source_id) {
  message("Trying direct World Bank API fallback with source = ", source_id)
  
  wgi_list <- imap(
    wgi_indicators,
    ~ fetch_wb_indicator(
      indicator_code = .x,
      indicator_name = .y,
      source_id = source_id,
      start_year = start_year_wgi,
      end_year = end_year_wgi
    )
  )
  
  reduce(
    wgi_list,
    full_join,
    by = c("iso3c", "country", "year")
  )
}

if (is.null(wgi_raw) || nrow(wgi_raw) == 0) {
  wgi_raw <- fetch_wgi_api_source(3)
}

if (is.null(wgi_raw) || nrow(wgi_raw) == 0) {
  wgi_raw <- fetch_wgi_api_source(75)
}

wgi_raw <- wgi_raw |>
  filter(
    !is.na(iso3c),
    iso3c != "",
    year >= start_year_wgi,
    year <= end_year_wgi
  )

if (nrow(wgi_raw) == 0) {
  stop(
    "No WGI data were downloaded. The World Bank API may be unavailable, ",
    "or the WGI source/indicator codes may have changed. ",
    "Try again later or download the WGI CSV manually."
  )
}

print(head(wgi_raw))
print(summary(wgi_raw))

# ------------------------------------------------------------------------------
# Preparing analysis samples ---------------------------------------------------
# ------------------------------------------------------------------------------
analysis_sample <- standardized_data_final |>
  ungroup() |>
  transmute(
    country = COUNTRY_NAME,
    iso3 = as.character(ISO_CODE_3),
    year = as.integer(YEAR),
    vdem_fragile_binary = VDEM_FRAGILE_IDEAL,
    wgi_fragile_binary = WGI_FRAGILE_IDEAL
  )

analysis_sample_wgi <- analysis_sample |>
  filter(
    year >= start_year_wgi,
    year <= end_year_wgi
  )

# ------------------------------------------------------------------------------
# Extracting and preparing V-Dem components ------------------------------------
# ------------------------------------------------------------------------------
vdem_components <- vdem_raw |>
  transmute(
    iso3 = as.character(country_text_id),
    country_vdem = country_name,
    year = as.integer(year),
    
    v2x_corr = v2x_corr,
    v2caviol = v2caviol,
    v2x_rule = v2x_rule,
    v2clrspct = v2clrspct
  )

df_vdem_components_raw <- analysis_sample |>
  left_join(
    vdem_components,
    by = c("iso3", "year")
  )

unmatched_vdem_components <- df_vdem_components_raw |>
  filter(
    is.na(v2x_corr) &
      is.na(v2caviol) &
      is.na(v2x_rule) &
      is.na(v2clrspct)
  ) |>
  distinct(country, iso3, year)

print(unmatched_vdem_components)

df_vdem_components <- df_vdem_components_raw |>
  mutate(
    corruption_raw = v2x_corr,
    nonstate_violence_raw = v2caviol,
    weak_rule_of_law_raw = -v2x_rule,
    weak_public_admin_raw = -v2clrspct
  ) |>
  mutate(
    corruption_frag =
      orient_to_fragility(corruption_raw, vdem_fragile_binary),
    
    nonstate_violence_frag =
      orient_to_fragility(nonstate_violence_raw, vdem_fragile_binary),
    
    weak_rule_of_law_frag =
      orient_to_fragility(weak_rule_of_law_raw, vdem_fragile_binary),
    
    weak_public_admin_frag =
      orient_to_fragility(weak_public_admin_raw, vdem_fragile_binary)
  ) |>
  mutate(
    `Corruption / predation` =
      rescale01(corruption_frag),
    
    `Violence / instability` =
      rescale01(nonstate_violence_frag),
    
    `Weak rule of law` =
      rescale01(weak_rule_of_law_frag),
    
    `Weak public administration` =
      rescale01(weak_public_admin_frag)
  )

vdem_component_cols <- c(
  "Corruption / predation",
  "Violence / instability",
  "Weak rule of law",
  "Weak public administration"
)

vdem_orientation_diagnostics <- df_vdem_components |>
  summarise(
    corr_corruption_with_fragile =
      safe_cor(`Corruption / predation`, vdem_fragile_binary),
    
    corr_violence_with_fragile =
      safe_cor(`Violence / instability`, vdem_fragile_binary),
    
    corr_rule_law_with_fragile =
      safe_cor(`Weak rule of law`, vdem_fragile_binary),
    
    corr_admin_with_fragile =
      safe_cor(`Weak public administration`, vdem_fragile_binary)
  )

print(vdem_orientation_diagnostics)

vdem_cor_matrix <- make_correlation_matrix(
  data = df_vdem_components,
  component_cols = vdem_component_cols
)

print(round(vdem_cor_matrix, 3))

# ------------------------------------------------------------------------------
# Extracting and preparing WGI components --------------------------------------
# ------------------------------------------------------------------------------
wgi_components <- wgi_raw |>
  filter(!is.na(iso3c)) |>
  transmute(
    iso3 = as.character(iso3c),
    country_wgi = country,
    year = as.integer(year),
    
    wgi_control_corruption = wgi_control_corruption,
    wgi_political_stability = wgi_political_stability,
    wgi_rule_of_law = wgi_rule_of_law,
    wgi_government_effectiveness = wgi_government_effectiveness
  )

df_wgi_components_raw <- analysis_sample_wgi |>
  left_join(
    wgi_components,
    by = c("iso3", "year")
  )

unmatched_wgi_components <- df_wgi_components_raw |>
  filter(
    is.na(wgi_control_corruption) &
      is.na(wgi_political_stability) &
      is.na(wgi_rule_of_law) &
      is.na(wgi_government_effectiveness)
  ) |>
  distinct(country, iso3, year)

print(unmatched_wgi_components)

df_wgi_components <- df_wgi_components_raw |>
  mutate(
    # WGI indicators are oriented so higher = better governance.
    # Reversing them so higher = greater institutional fragility.
    corruption_raw = -wgi_control_corruption,
    violence_raw = -wgi_political_stability,
    weak_rule_of_law_raw = -wgi_rule_of_law,
    weak_public_admin_raw = -wgi_government_effectiveness
  ) |>
  mutate(
    corruption_frag =
      orient_to_fragility(corruption_raw, wgi_fragile_binary),
    
    violence_frag =
      orient_to_fragility(violence_raw, wgi_fragile_binary),
    
    weak_rule_of_law_frag =
      orient_to_fragility(weak_rule_of_law_raw, wgi_fragile_binary),
    
    weak_public_admin_frag =
      orient_to_fragility(weak_public_admin_raw, wgi_fragile_binary)
  ) |>
  mutate(
    `Corruption / predation` =
      rescale01(corruption_frag),
    
    `Violence / instability` =
      rescale01(violence_frag),
    
    `Weak rule of law` =
      rescale01(weak_rule_of_law_frag),
    
    `Weak public administration` =
      rescale01(weak_public_admin_frag)
  )

wgi_component_cols <- c(
  "Corruption / predation",
  "Violence / instability",
  "Weak rule of law",
  "Weak public administration"
)

wgi_orientation_diagnostics <- df_wgi_components |>
  summarise(
    corr_corruption_with_fragile =
      safe_cor(`Corruption / predation`, wgi_fragile_binary),
    
    corr_violence_with_fragile =
      safe_cor(`Violence / instability`, wgi_fragile_binary),
    
    corr_rule_law_with_fragile =
      safe_cor(`Weak rule of law`, wgi_fragile_binary),
    
    corr_admin_with_fragile =
      safe_cor(`Weak public administration`, wgi_fragile_binary)
  )

print(wgi_orientation_diagnostics)

wgi_cor_matrix <- make_correlation_matrix(
  data = df_wgi_components,
  component_cols = wgi_component_cols
)

print(round(wgi_cor_matrix, 3))

# ------------------------------------------------------------------------------
# Building LaTeX-ready matrices ------------------------------------------------
# ------------------------------------------------------------------------------
vdem_latex_matrix <- lower_triangle_for_latex(
  cor_mat = vdem_cor_matrix,
  panel_label = "Panel A: V-Dem components"
)

wgi_latex_matrix <- lower_triangle_for_latex(
  cor_mat = wgi_cor_matrix,
  panel_label = "Panel B: WGI components"
)

# ------------------------------------------------------------------------------
# Printing separate LaTeX tables -----------------------------------------------
# ------------------------------------------------------------------------------
knitr::kable(
  vdem_latex_matrix |> select(-panel),
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Pairwise Correlations among V-Dem Fragility Components",
  label = "fragility_component_correlation_matrix_vdem"
)

knitr::kable(
  wgi_latex_matrix |> select(-panel),
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Pairwise Correlations among WGI Fragility Components",
  label = "fragility_component_correlation_matrix_wgi"
)

# ------------------------------------------------------------------------------
# Manually assembling two-panel LaTeX table ------------------------------------
# ------------------------------------------------------------------------------
latex_two_panel <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Pairwise Correlations among Fragility Components}",
  "\\label{tab:fragility_component_correlation_matrix}",
  "\\begingroup",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\renewcommand{\\arraystretch}{1.15}",
  "",
  "\\begin{tabular}{lrrrr}",
  "\\toprule",
  "Component & \\makecell{Corruption /\\\\predation} & \\makecell{Violence /\\\\instability} & \\makecell{Weak rule\\\\of law} & \\makecell{Weak public\\\\administration} \\\\",
  "\\midrule",
  "\\multicolumn{5}{l}{\\textit{Panel A: V-Dem components}} \\\\",
  paste0(
    vdem_latex_matrix$Component, " & ",
    vdem_latex_matrix[[2]], " & ",
    vdem_latex_matrix[[3]], " & ",
    vdem_latex_matrix[[4]], " & ",
    vdem_latex_matrix[[5]], " \\\\"
  ),
  "\\addlinespace",
  "\\multicolumn{5}{l}{\\textit{Panel B: WGI components}} \\\\",
  paste0(
    wgi_latex_matrix$Component, " & ",
    wgi_latex_matrix[[2]], " & ",
    wgi_latex_matrix[[3]], " & ",
    wgi_latex_matrix[[4]], " & ",
    wgi_latex_matrix[[5]], " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{0.2cm}",
  "",
  "\\begin{minipage}{0.92\\textwidth}",
  "\\footnotesize",
  "\\emph{Notes:} The table reports Pearson correlations among reoriented fragility components. Panel A uses V-Dem components; Panel B uses WGI components. All variables are scaled so that higher values indicate greater institutional fragility. For V-Dem, corruption/predation is based on the political corruption index, violence/instability is based on the V-Dem political violence measure, weak rule of law reverses the rule-of-law index, and weak public administration reverses the rigorous and impartial public administration measure. For WGI, corruption/predation reverses control of corruption, violence/instability reverses political stability and absence of violence/terrorism, weak rule of law reverses rule of law, and weak public administration reverses government effectiveness. Correlations are computed using pairwise complete observations. The table is intended as a measurement-validity diagnostic: the components are positively correlated, consistent with a common institutional syndrome, but not perfectly collinear.",
  "\\end{minipage}",
  "",
  "\\endgroup",
  "\\end{table}"
)

writeLines(
  latex_two_panel,
  con = "fragility_component_correlation_matrix_two_panel.tex"
)

cat(paste(latex_two_panel, collapse = "\n"))

# ==============================================================================
# PCA Loadings and Variance Explained
# V-Dem and WGI fragility components
#
# Assumes the previous pipeline already created:
#   df_vdem_components
#   df_wgi_components
#   vdem_component_cols
#   wgi_component_cols
#
# Higher component values = greater institutional fragility
# ==============================================================================
# ------------------------------------------------------------------------------
# Helpers: estimating PCA and extract diagnostics ------------------------------
# ------------------------------------------------------------------------------
run_component_pca <- function(data, component_cols, source_label) {
  
  pca_data <- data |>
    select(all_of(component_cols)) |>
    drop_na()
  
  pca_fit <- prcomp(
    pca_data,
    center = TRUE,
    scale. = TRUE
  )
  
  eigenvalues <- pca_fit$sdev^2
  variance_share <- eigenvalues / sum(eigenvalues)
  
  pc1_loadings <- pca_fit$rotation[, 1]
  
  # PCA signs are arbitrary. Reorient PC1 so positive loadings
  # correspond to greater fragility on average.
  if (mean(pc1_loadings, na.rm = TRUE) < 0) {
    pc1_loadings <- -pc1_loadings
  }
  
  loadings_table <- tibble(
    source = source_label,
    component = names(pc1_loadings),
    pc1_loading = as.numeric(pc1_loadings)
  )
  
  summary_table <- tibble(
    source = source_label,
    n_complete_country_years = nrow(pca_data),
    pc1_eigenvalue = eigenvalues[1],
    pc1_variance_explained = variance_share[1],
    pc2_variance_explained = variance_share[2],
    pc3_variance_explained = variance_share[3],
    pc4_variance_explained = variance_share[4]
  )
  
  list(
    pca_fit = pca_fit,
    loadings = loadings_table,
    summary = summary_table
  )
}

# ------------------------------------------------------------------------------
# Running PCA separately for V-Dem and WGI -------------------------------------
# ------------------------------------------------------------------------------
pca_vdem <- run_component_pca(
  data = df_vdem_components,
  component_cols = vdem_component_cols,
  source_label = "V-Dem"
)

pca_wgi <- run_component_pca(
  data = df_wgi_components,
  component_cols = wgi_component_cols,
  source_label = "WGI"
)

# ------------------------------------------------------------------------------
# Combining outputs ------------------------------------------------------------
# ------------------------------------------------------------------------------
pca_loadings <- bind_rows(
  pca_vdem$loadings,
  pca_wgi$loadings
)

pca_summary <- bind_rows(
  pca_vdem$summary,
  pca_wgi$summary
)

print(pca_loadings)
print(pca_summary)

# ------------------------------------------------------------------------------
# LaTeX table ------------------------------------------------------------------
# ------------------------------------------------------------------------------
pca_loadings_wide <- pca_loadings |>
  mutate(
    pc1_loading_print = sprintf("%.3f", pc1_loading)
  ) |>
  select(source, component, pc1_loading_print) |>
  pivot_wider(
    names_from = source,
    values_from = pc1_loading_print
  )

pca_summary_print <- pca_summary |>
  transmute(
    source,
    n_complete_country_years,
    pc1_eigenvalue_print = sprintf("%.3f", pc1_eigenvalue),
    pc1_variance_print = percent(pc1_variance_explained, accuracy = 0.1),
    pc2_variance_print = percent(pc2_variance_explained, accuracy = 0.1),
    pc3_variance_print = percent(pc3_variance_explained, accuracy = 0.1),
    pc4_variance_print = percent(pc4_variance_explained, accuracy = 0.1)
  )

print(pca_loadings_wide)
print(pca_summary_print)

# ------------------------------------------------------------------------------
# Combined appendix table: loadings + PC1 variance explained -------------------
# ------------------------------------------------------------------------------
latex_pca_table <- pca_loadings |>
  mutate(
    component = factor(
      component,
      levels = c(
        "Corruption / predation",
        "Violence / instability",
        "Weak rule of law",
        "Weak public administration"
      )
    ),
    pc1_loading_print = sprintf("%.3f", pc1_loading)
  ) |>
  select(source, component, pc1_loading_print) |>
  pivot_wider(
    names_from = source,
    values_from = pc1_loading_print
  ) |>
  arrange(component) |>
  mutate(
    component = as.character(component)
  ) |>
  select(
    Component = component,
    `V-Dem PC1 loading` = `V-Dem`,
    `WGI PC1 loading` = WGI
  )

knitr::kable(
  latex_pca_table,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "PCA Loadings for Fragility Components",
  label = "tab:fragility_component_pca_loadings"
)

# ------------------------------------------------------------------------------
# Assembled LaTeX table with notes ---------------------------------------------
# ------------------------------------------------------------------------------
vdem_pc1_share <- pca_summary |>
  filter(source == "V-Dem") |>
  pull(pc1_variance_explained)

wgi_pc1_share <- pca_summary |>
  filter(source == "WGI") |>
  pull(pc1_variance_explained)

vdem_n <- pca_summary |>
  filter(source == "V-Dem") |>
  pull(n_complete_country_years)

wgi_n <- pca_summary |>
  filter(source == "WGI") |>
  pull(n_complete_country_years)

latex_pca_manual <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{PCA Loadings and Variance Explained by the First Component}",
  "\\label{tab:fragility_component_pca_loadings}",
  "\\begingroup",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{6pt}",
  "\\renewcommand{\\arraystretch}{1.15}",
  "",
  "\\begin{tabular}{lrr}",
  "\\toprule",
  "Component & \\makecell{V-Dem\\\\PC1 loading} & \\makecell{WGI\\\\PC1 loading} \\\\",
  "\\midrule",
  paste0(
    latex_pca_table$Component, " & ",
    latex_pca_table$`V-Dem PC1 loading`, " & ",
    latex_pca_table$`WGI PC1 loading`, " \\\\"
  ),
  "\\midrule",
  paste0(
    "Complete country-years & ",
    format(vdem_n, big.mark = ","), " & ",
    format(wgi_n, big.mark = ","), " \\\\"
  ),
  paste0(
    "Variance explained by PC1 & ",
    percent(vdem_pc1_share, accuracy = 0.1), " & ",
    percent(wgi_pc1_share, accuracy = 0.1), " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{0.2cm}",
  "",
  "\\begin{minipage}{0.92\\textwidth}",
  "\\footnotesize",
  "\\emph{Notes:} The table reports first-component PCA loadings for the four reoriented fragility components. All variables are scaled so that higher values indicate greater institutional fragility. PCA is estimated separately for the V-Dem and WGI operationalizations using standardized component scores. The sign of the first component is oriented so that positive loadings correspond to greater fragility. The final row reports the share of total variance explained by the first principal component.",
  "\\end{minipage}",
  "",
  "\\endgroup",
  "\\end{table}"
)

writeLines(
  latex_pca_manual,
  con = "fragility_component_pca_loadings_table.tex"
)

cat(paste(latex_pca_manual, collapse = "\n"))

# ==============================================================================
# Cross-Dimensional Range Diagnostic
# V-Dem and WGI fragility components
#
# Goal:
#   Assess whether country-years often combine very high fragility
#   on one dimension with very low fragility on another.
#
# Range:
#   max(component scores) - min(component scores)
#
# Higher range = more cross-dimensional imbalance
# Lower range  = components move more jointly within a country-year
# ==============================================================================
# ------------------------------------------------------------------------------
# Helpers ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
make_range_diagnostic <- function(data, component_cols, source_label) {
  
  d <- data |>
    select(
      country,
      iso3,
      year,
      all_of(component_cols)
    ) |>
    drop_na(all_of(component_cols))
  
  component_matrix <- as.matrix(d[, component_cols])
  
  max_idx <- apply(component_matrix, 1, which.max)
  min_idx <- apply(component_matrix, 1, which.min)
  
  d |>
    mutate(
      source = source_label,
      component_min = apply(component_matrix, 1, min),
      component_max = apply(component_matrix, 1, max),
      component_mean = rowMeans(component_matrix),
      component_sd = apply(component_matrix, 1, sd),
      cross_dimensional_range = component_max - component_min,
      most_fragile_component = component_cols[max_idx],
      least_fragile_component = component_cols[min_idx]
    )
}

# ------------------------------------------------------------------------------
# Diagnostics for V-Dem and WGI ------------------------------------------------
# ------------------------------------------------------------------------------
vdem_range <- make_range_diagnostic(
  data = df_vdem_components,
  component_cols = vdem_component_cols,
  source_label = "V-Dem"
)

wgi_range <- make_range_diagnostic(
  data = df_wgi_components,
  component_cols = wgi_component_cols,
  source_label = "WGI"
)

range_data <- bind_rows(
  vdem_range,
  wgi_range
)

print(head(range_data))

# ------------------------------------------------------------------------------
# Summary statistics -----------------------------------------------------------
# ------------------------------------------------------------------------------
range_summary <- range_data |>
  group_by(source) |>
  summarise(
    country_years = n(),
    mean_range = mean(cross_dimensional_range),
    sd_range = sd(cross_dimensional_range),
    median_range = median(cross_dimensional_range),
    p75_range = quantile(cross_dimensional_range, 0.75),
    p90_range = quantile(cross_dimensional_range, 0.90),
    p95_range = quantile(cross_dimensional_range, 0.95),
    p99_range = quantile(cross_dimensional_range, 0.99),
    share_range_le_025 = mean(cross_dimensional_range <= 0.25),
    share_range_le_050 = mean(cross_dimensional_range <= 0.50),
    share_range_ge_050 = mean(cross_dimensional_range >= 0.50),
    share_range_ge_075 = mean(cross_dimensional_range >= 0.75),
    .groups = "drop"
  )

print(range_summary)

# ------------------------------------------------------------------------------
# Binned distribution ----------------------------------------------------------
# ------------------------------------------------------------------------------
range_bins <- range_data |>
  mutate(
    range_bin = cut(
      cross_dimensional_range,
      breaks = c(0, 0.10, 0.25, 0.50, 0.75, 1.00),
      include.lowest = TRUE,
      right = TRUE,
      labels = c(
        "0.00--0.10",
        "0.10--0.25",
        "0.25--0.50",
        "0.50--0.75",
        "0.75--1.00"
      )
    )
  ) |>
  count(source, range_bin, name = "country_years") |>
  group_by(source) |>
  mutate(
    share = country_years / sum(country_years),
    share_print = percent(share, accuracy = 0.1)
  ) |>
  ungroup()

print(range_bins)

# ------------------------------------------------------------------------------
# Most discordant country-years ------------------------------------------------
# ------------------------------------------------------------------------------
top_discordant_cases <- range_data |>
  group_by(source) |>
  slice_max(
    order_by = cross_dimensional_range,
    n = 25,
    with_ties = FALSE
  ) |>
  ungroup() |>
  arrange(source, desc(cross_dimensional_range)) |>
  select(
    source,
    country,
    iso3,
    year,
    cross_dimensional_range,
    component_min,
    component_max,
    least_fragile_component,
    most_fragile_component,
    all_of(vdem_component_cols)
  )

print(top_discordant_cases)

# ------------------------------------------------------------------------------
# LaTeX summary table ----------------------------------------------------------
# ------------------------------------------------------------------------------
range_summary_latex <- range_summary |>
  transmute(
    Source = source,
    `Country-years` = format(country_years, big.mark = ","),
    `Mean range` = sprintf("%.3f", mean_range),
    `Median range` = sprintf("%.3f", median_range),
    `75th pct.` = sprintf("%.3f", p75_range),
    `90th pct.` = sprintf("%.3f", p90_range),
    `95th pct.` = sprintf("%.3f", p95_range),
    `Share range $\\leq$ 0.25` = percent(share_range_le_025, accuracy = 0.1),
    `Share range $\\leq$ 0.50` = percent(share_range_le_050, accuracy = 0.1),
    `Share range $\\geq$ 0.75` = percent(share_range_ge_075, accuracy = 0.1)
  )

print(range_summary_latex)

knitr::kable(
  range_summary_latex,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Distribution of Cross-Dimensional Ranges among Fragility Components",
  label = "fragility_component_cross_dimensional_range_summary"
)

# ------------------------------------------------------------------------------
# LaTeX binned distribution table ----------------------------------------------
# ------------------------------------------------------------------------------
range_bins_latex <- range_bins |>
  select(source, range_bin, country_years, share_print) |>
  pivot_wider(
    names_from = source,
    values_from = c(country_years, share_print),
    names_glue = "{source}_{.value}"
  ) |>
  transmute(
    `Cross-dimensional range` = as.character(range_bin),
    `V-Dem country-years` = format(`V-Dem_country_years`, big.mark = ","),
    `V-Dem share` = `V-Dem_share_print`,
    `WGI country-years` = format(WGI_country_years, big.mark = ","),
    `WGI share` = WGI_share_print
  )

print(range_bins_latex)

knitr::kable(
  range_bins_latex,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Binned Distribution of Cross-Dimensional Ranges among Fragility Components",
  label = "fragility_component_cross_dimensional_range_bins"
)

# ------------------------------------------------------------------------------
# Plot: distribution of cross-dimensional ranges -------------------------------
# ------------------------------------------------------------------------------
dark_red <- "#A12A2A"

range_medians <- range_data |>
  group_by(source) |>
  summarise(
    median_range = median(cross_dimensional_range),
    .groups = "drop"
  )

p_range_distribution <- ggplot(
  range_data,
  aes(x = cross_dimensional_range)
) +
  geom_histogram(
    binwidth = 0.05,
    boundary = 0,
    closed = "right",
    fill = dark_red,
    colour = "white",
    linewidth = 0.2
  ) +
  geom_vline(
    data = range_medians,
    aes(xintercept = median_range),
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey25"
  ) +
  facet_wrap(~ source, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    x = "Cross-dimensional range",
    y = "Country-years",
    title = "Distribution of cross-dimensional ranges among fragility components",
    subtitle = "Range = maximum component score minus minimum component score within each country-year"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey25"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.major.x = element_line(colour = "grey85")
  )

p_range_distribution

ggsave(
  filename = "fragility_component_cross_dimensional_range_distribution.png",
  plot = p_range_distribution,
  width = 9,
  height = 7,
  dpi = 300
)
