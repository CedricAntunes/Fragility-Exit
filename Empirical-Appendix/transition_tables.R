# Transition Matrices for Binary Fragility Status
# Horizons: 1 to 5 years
# Author: Cedric Antunes (FGV-CEPESP) ------------------------------------------
# Date: May, 2026 --------------------------------------------------------------

# Assumes standardized_data_final is loaded in the environment

# Required packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(knitr)
  library(kableExtra)
  library(scales)
  library(stringr)
})

# ------------------------------------------------------------------------------
# Defining binary fragility outcomes -------------------------------------------
# ------------------------------------------------------------------------------
fragility_outcomes <- c(
  "VDEM_FRAGILE_IDEAL",
  "VDEM_NORL_FRAGILE_IDEAL",
  "WGI_FRAGILE_IDEAL",
  "WGI_NORL_FRAGILE_IDEAL"
)

fragility_labels <- c(
  "VDEM_FRAGILE_IDEAL" = "V-Dem fragility status",
  "VDEM_NORL_FRAGILE_IDEAL" = "V-Dem fragility status, no rule of law",
  "WGI_FRAGILE_IDEAL" = "WGI fragility status",
  "WGI_NORL_FRAGILE_IDEAL" = "WGI fragility status, no rule of law"
)

transition_horizons <- 1:5

# ------------------------------------------------------------------------------
# Data checks ------------------------------------------------------------------
# ------------------------------------------------------------------------------
missing_outcomes <- setdiff(
  fragility_outcomes,
  names(standardized_data_final)
)

fragility_outcomes <- intersect(
  fragility_outcomes,
  names(standardized_data_final)
)

fragility_labels <- fragility_labels[fragility_outcomes]

# ------------------------------------------------------------------------------
# Robust converters (just in case) ---------------------------------------------
# ------------------------------------------------------------------------------
to_year_int <- function(x) {
  if (inherits(x, "Date")) {
    return(as.integer(format(x, "%Y")))
  }
  
  if (inherits(x, "POSIXct") || inherits(x, "POSIXlt")) {
    return(as.integer(format(x, "%Y")))
  }
  
  x_chr <- as.character(x)
  as.integer(str_extract(x_chr, "\\d{4}"))
}

to_binary_int <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  
  if (is.character(x)) {
    x_clean <- str_to_lower(str_trim(x))
    
    return(case_when(
      x_clean %in% c("1", "fragile", "yes", "true") ~ 1L,
      x_clean %in% c("0", "not fragile", "non fragile", "non-fragile",
                     "no", "false") ~ 0L,
      TRUE ~ NA_integer_
    ))
  }
  
  if (is.logical(x)) {
    return(as.integer(x))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    return(as.integer(x))
  }
  
  as.integer(x)
}

# ------------------------------------------------------------------------------
# Preparing the data -----------------------------------------------------------
# ------------------------------------------------------------------------------
transition_base <- standardized_data_final |>
  ungroup() |>
  mutate(
    YEAR_clean = to_year_int(YEAR),
    ISO_CODE_3 = as.character(ISO_CODE_3)
  ) |>
  select(
    ISO_CODE_3,
    COUNTRY_NAME,
    YEAR_clean,
    all_of(fragility_outcomes)
  ) |>
  mutate(
    across(
      all_of(fragility_outcomes),
      ~ to_binary_int(.x)
    )
  )

# Diagnostic -------------------------------------------------------------------
transition_base |>
  summarise(
    rows = n(),
    countries = n_distinct(ISO_CODE_3),
    year_min = min(YEAR_clean, na.rm = TRUE),
    year_max = max(YEAR_clean, na.rm = TRUE),
    across(
      all_of(fragility_outcomes),
      list(
        non_missing = ~ sum(!is.na(.x)),
        values = ~ paste(sort(unique(.x)), collapse = ", ")
      ),
      .names = "{.col}_{.fn}"
    )
  )

# ------------------------------------------------------------------------------
# Function to create transition data using self-join ---------------------------
# ------------------------------------------------------------------------------
make_transition_data_h <- function(data, outcome_var, outcome_label, h) {
  
  end_data <- data %>%
    transmute(
      ISO_CODE_3,
      COUNTRY_NAME,
      YEAR_end = YEAR_clean,
      fragile_end = .data[[outcome_var]]
    )
  
  start_data <- data %>%
    transmute(
      ISO_CODE_3,
      YEAR_start = YEAR_clean,
      fragile_start = .data[[outcome_var]]
    )
  
  end_data %>%
    mutate(
      YEAR_start = YEAR_end - h
    ) %>%
    left_join(
      start_data,
      by = c("ISO_CODE_3", "YEAR_start")
    ) %>%
    filter(
      !is.na(fragile_end),
      !is.na(fragile_start)
    ) %>%
    mutate(
      Outcome = outcome_label,
      Horizon = h
    )
}

# ------------------------------------------------------------------------------
# Building transition data for all outcomes and horizons -----------------------
# ------------------------------------------------------------------------------
transition_data_all <- bind_rows(
  lapply(
    transition_horizons,
    function(h) {
      bind_rows(
        lapply(
          fragility_outcomes,
          function(v) {
            make_transition_data_h(
              data = transition_base,
              outcome_var = v,
              outcome_label = fragility_labels[[v]],
              h = h
            )
          }
        )
      )
    }
  )
)

# sanity check
transition_data_all |>
  count(Horizon, 
        Outcome, 
        fragile_start, 
        fragile_end)

# ------------------------------------------------------------------------------
# Creating transition matrix for all horizons ----------------------------------
# ------------------------------------------------------------------------------
status_levels <- c("Fragile", "Not fragile")

transition_matrix_all <- transition_data_all |>
  mutate(
    `Status at start` = if_else(fragile_start == 1, "Fragile", "Not fragile"),
    `Status at end`   = if_else(fragile_end == 1, "Fragile", "Not fragile"),
    Outcome = factor(Outcome, levels = unname(fragility_labels)),
    `Status at start` = factor(`Status at start`, levels = status_levels),
    `Status at end` = factor(`Status at end`, levels = status_levels)
  ) |>
  count(
    Horizon,
    Outcome,
    `Status at start`,
    `Status at end`,
    name = "n"
  ) |>
  complete(
    Horizon = transition_horizons,
    Outcome = factor(unname(fragility_labels), levels = unname(fragility_labels)),
    `Status at start` = factor(status_levels, levels = status_levels),
    `Status at end` = factor(status_levels, levels = status_levels),
    fill = list(n = 0)
  ) |>
  group_by(
    Horizon,
    Outcome,
    `Status at start`
  ) |>
  mutate(
    denominator = sum(n),
    Probability_raw = if_else(denominator > 0, n / denominator, NA_real_)
  ) |>
  ungroup() |>
  mutate(
    Transitions = comma(n),
    Probability = percent(Probability_raw, accuracy = 0.1)
  ) |>
  select(
    Horizon,
    Outcome,
    `Status at start`,
    `Status at end`,
    Transitions,
    Probability
  ) |>
  arrange(
    Horizon,
    Outcome,
    `Status at start`,
    `Status at end`
  )

transition_matrix_all

# ------------------------------------------------------------------------------
# LaTeX tables -----------------------------------------------------------------
# ------------------------------------------------------------------------------
make_transition_matrix_latex <- function(matrix_data, h) {
  
  start_label <- paste0("Status at $t-", h, "$")
  
  horizon_label <- case_when(
    h == 1 ~ "One-Year",
    h == 2 ~ "Two-Year",
    h == 3 ~ "Three-Year",
    h == 4 ~ "Four-Year",
    h == 5 ~ "Five-Year",
    TRUE ~ paste0(h, "-Year")
  )
  
  caption_text <- paste0(
    horizon_label,
    " Transition Matrix for Binary Fragility Status"
  )
  
  label_text <- paste0("fragility_transition_matrix_", h, "yr")
  
  note_text <- if (h == 1) {
    paste(
      "The table uses observed one-year country pairs with non-missing fragility status at both endpoints.",
      "It is not restricted to countries that change status;",
      "countries that are always fragile or never fragile contribute to the persistence probabilities."
    )
  } else {
    paste(
      "The table uses observed country windows with non-missing fragility status at both endpoints.",
      "It is not restricted to countries that change status;",
      "countries that are always fragile or never fragile contribute to the persistence probabilities.",
      paste0("Probabilities are endpoint transition probabilities over the ", h, "-year window.")
    )
  }
  
  # Prepare data without Outcome column; Outcome will be used as panel labels
  table_data <- matrix_data %>%
    filter(Horizon == h) %>%
    select(
      Outcome,
      `Status at start`,
      `Status at end`,
      Transitions,
      Probability
    ) %>%
    arrange(
      Outcome,
      `Status at start`,
      `Status at end`
    )
  
  panel_index <- table_data %>%
    count(Outcome) %>%
    mutate(Outcome = as.character(Outcome)) %>%
    deframe()
  
  table_data_clean <- table_data %>%
    select(-Outcome)
  
  table_data_clean %>%
    kbl(
      format = "latex",
      booktabs = TRUE,
      caption = caption_text,
      label = label_text,
      col.names = c(
        start_label,
        "Status at $t$",
        "Transitions",
        "Probability"
      ),
      escape = FALSE
    ) %>%
    kable_styling(
      latex_options = c("hold_position"),
      font_size = 8
    ) %>%
    pack_rows(
      index = panel_index,
      italic = TRUE,
      bold = FALSE,
      escape = FALSE
    ) %>%
    footnote(
      general = note_text,
      general_title = "Note: ",
      footnote_as_chunk = TRUE,
      escape = FALSE
    )
}

# ------------------------------------------------------------------------------
# LaTeX tables for horizons 1 to 5 ---------------------------------------------
# ------------------------------------------------------------------------------
transition_matrix_latex <- lapply(
  transition_horizons,
  function(h) {
    make_transition_matrix_latex(
      matrix_data = transition_matrix_all,
      h = h
    )
  }
)

names(transition_matrix_latex) <- paste0(
  "transition_matrix_",
  transition_horizons,
  "yr_latex"
)

# Print individual tables
transition_matrix_latex$transition_matrix_1yr_latex
transition_matrix_latex$transition_matrix_2yr_latex
transition_matrix_latex$transition_matrix_3yr_latex
transition_matrix_latex$transition_matrix_4yr_latex
transition_matrix_latex$transition_matrix_5yr_latex
