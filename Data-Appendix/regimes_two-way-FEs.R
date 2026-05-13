# Two-Way Static FEs Regime Estimates ------------------------------------------
# Author: Cedric Antunes (FGV-CEPESP) ------------------------------------------
# Date: May 13, 2026 -----------------------------------------------------------

# Assumes standardized_data_final exists in the environment

# Required packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(fixest)
  library(forcats)
  library(scales)
  library(glue)
})

# ------------------------------------------------------------------------------
# Colours ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
main_red  <- "#B40000"
dark_red  <- "#7A0000"
light_red <- "#E6B3B3"

# ------------------------------------------------------------------------------
# Regime-measure specifications ------------------------------------------------
# ------------------------------------------------------------------------------
regime_specs <- tribble(
  ~var,                              ~label,                                ~family,      ~scale_type,              ~standardize,
  "INSTITUTIONAL_AUTOCRACY_SCORE",   "Institutional autocracy (Polity-V)",   "Autocracy",  "Continuous scores",      TRUE,
  "INSTITUTIONAL_DEMOCRACY_SCORE",   "Institutional democracy (Polity-V)",   "Democracy",  "Continuous scores",      TRUE,
  "LIBERAL_DEMOCRACY_SCORE",         "Liberal democracy (V-Dem)",           "Democracy",  "Continuous scores",      TRUE,
  "ELECTORAL_DEMOCRACY_SCORE",       "Electoral democracy (V-Dem)",         "Democracy",  "Continuous scores",      TRUE,
  "FH_DEMOCRACY",                    "Democracy (FH)",                      "Democracy",  "Binary indicators",      FALSE,
  "DD_DEMOCRACY",                    "Democracy (DD)",                      "Democracy",  "Binary indicators",      FALSE,
  "BOIX_DEMOCRACY",                  "Democracy (Boix)",                    "Democracy",  "Binary indicators",      FALSE,
  "GWF_DEMOCRACY",                   "Democracy (GWF)",                     "Democracy",  "Binary indicators",      FALSE,
  "GWF_AUTOCRACY",                   "Autocracy (GWF)",                     "Autocracy",  "Binary indicators",      FALSE,
  "GWF_AUTOCRACY_PARTY",             "Party autocracy (GWF)",               "Autocracy",  "Binary indicators",      FALSE,
  "GWF_MILITARY",                    "Military autocracy (GWF)",            "Autocracy",  "Binary indicators",      FALSE,
  "GWF_PERSONALIST",                 "Personalist autocracy (GWF)",         "Autocracy",  "Binary indicators",      FALSE
)

df_plot <- standardized_data_final

# ------------------------------------------------------------------------------
# Outcome specifications -------------------------------------------------------
# ------------------------------------------------------------------------------
plot_specs <- tribble(
  ~outcome_var,                    ~outcome_label,                 ~include_gdp, ~file_stub,
  "VDEM_FRAGILE_IDEAL",            "V-Dem",                       FALSE,        "vdem",
  "VDEM_FRAGILE_IDEAL",            "V-Dem + GDP per capita",      TRUE,         "vdem_gdp",
  "VDEM_NORL_FRAGILE_IDEAL",       "V-Dem (no rule of law)",      FALSE,        "vdem_norl",
  "VDEM_NORL_FRAGILE_IDEAL",       "V-Dem (no rule of law) + GDP per capita", TRUE, "vdem_norl_gdp",
  "WGI_FRAGILE_IDEAL",             "WGI",                         FALSE,        "wgi",
  "WGI_FRAGILE_IDEAL",             "WGI + GDP per capita",        TRUE,         "wgi_gdp",
  "WGI_NORL_FRAGILE_IDEAL",        "WGI (no rule of law)",        FALSE,        "wgi_norl",
  "WGI_NORL_FRAGILE_IDEAL",        "WGI (no rule of law) + GDP per capita", TRUE, "wgi_norl_gdp"
)

# ------------------------------------------------------------------------------
# Function: estimate one model and extract one coefficient ---------------------
# ------------------------------------------------------------------------------
fit_one_regime_model <- function(data, outcome_var, include_gdp,
                                 var, label, family, scale_type, standardize) {

  term <- var
  
  rhs <- if (include_gdp) {
    paste(term, "+ LOG_GDP_PER_CAPITA")
  } else {
    term
  }
  
  fml <- as.formula(
    paste0(outcome_var, " ~ ", rhs, " | ISO_CODE_3 + YEAR")
  )
  
  mod <- feols(
    fml,
    data = data,
    vcov = ~ ISO_CODE_3 + YEAR,
    warn = FALSE
  )
  
  ct <- fixest::coeftable(mod, vcov = ~ ISO_CODE_3 + YEAR)
  
  if (!term %in% rownames(ct)) {
    return(tibble(
      outcome_var = outcome_var,
      include_gdp = include_gdp,
      var = var,
      label = label,
      family = family,
      scale_type = scale_type,
      estimate = NA_real_,
      se = NA_real_,
      p_value = NA_real_,
      conf95_low = NA_real_,
      conf95_high = NA_real_,
      conf90_low = NA_real_,
      conf90_high = NA_real_,
      nobs = nobs(mod)
    ))
  }
  
  est <- unname(ct[term, "Estimate"])
  se  <- unname(ct[term, "Std. Error"])
  p_col <- grep("^Pr", colnames(ct), value = TRUE)[1]
  pval <- unname(ct[term, p_col])
  
  tibble(
    outcome_var = outcome_var,
    include_gdp = include_gdp,
    var = var,
    label = label,
    family = family,
    scale_type = scale_type,
    estimate = est,
    se = se,
    p_value = pval,
    conf95_low = est - qnorm(0.975) * se,
    conf95_high = est + qnorm(0.975) * se,
    conf90_low = est - qnorm(0.95)  * se,
    conf90_high = est + qnorm(0.95) * se,
    nobs = nobs(mod)
  )
}

# ------------------------------------------------------------------------------
# Function: building coefficient dataframe for one plot ------------------------
# ------------------------------------------------------------------------------
build_coef_df <- function(data, outcome_var, include_gdp) {
  
  pmap_dfr(
    regime_specs,
    function(var, label, family, scale_type, standardize) {
      fit_one_regime_model(
        data = data,
        outcome_var = outcome_var,
        include_gdp = include_gdp,
        var = var,
        label = label,
        family = family,
        scale_type = scale_type,
        standardize = standardize
      )
    }
  ) |>
    filter(!is.na(estimate)) |>
    mutate(
      family = factor(family, levels = c("Democracy", "Autocracy")),
      scale_type = factor(scale_type, levels = c("Continuous scores", "Binary indicators"))
    )
}

# ------------------------------------------------------------------------------
# Building all coefficient dataframes first ------------------------------------
# ------------------------------------------------------------------------------

coef_results <- plot_specs |>
  mutate(
    coef_df = pmap(
      list(outcome_var, include_gdp),
      ~ build_coef_df(
        data = df_plot,
        outcome_var = ..1,
        include_gdp = ..2
      )
    )
  )

# ------------------------------------------------------------------------------
# Computing common x-axis limits across all 8 plots ----------------------------
# ------------------------------------------------------------------------------
all_coefs <- bind_rows(coef_results$coef_df)

x_limits <- range(
  c(all_coefs$conf95_low, all_coefs$conf95_high),
  na.rm = TRUE
)

# Adding small padding
x_pad <- 0.05 * diff(x_limits)
x_limits <- c(x_limits[1] - x_pad, x_limits[2] + x_pad)

# ------------------------------------------------------------------------------
# Plotting ---------------------------------------------------------------------
# ------------------------------------------------------------------------------
make_coef_plot <- function(coef_df, outcome_label, include_gdp, x_limits = NULL) {
  
  coef_df <- coef_df |>
    mutate(label = fct_reorder(label, estimate))
  
  subtitle_text <- if (include_gdp) {
    "Two-way FEs; adjusted for logged GDP per capita; SEs clustered by country and year"
  } else {
    "Two-way FEs; SEs clustered by country and year"
  }
  
  caption_text <- paste()
  
  p <- ggplot(coef_df, aes(x = estimate, y = label)) +
    geom_vline(
      xintercept = 0,
      linetype = "twodash",
      color = "grey50",
      linewidth = 0.4
    ) +
    geom_errorbarh(
      aes(xmin = conf95_low, xmax = conf95_high),
      height = 0,
      linewidth = 1.2,
      color = light_red,
      alpha = 0.95
    ) +
    geom_errorbarh(
      aes(xmin = conf90_low, xmax = conf90_high),
      height = 0,
      linewidth = 2.2,
      color = main_red,
      alpha = 0.75
    ) +
    geom_point(
      aes(shape = family),
      size = 2.8,
      color = dark_red
    ) +
    facet_wrap(
      ~ scale_type,
      scales = "free_y",
      ncol = 1
    ) +
    scale_x_continuous(
      labels = label_number(accuracy = 0.01),
      limits = x_limits
    ) +
    labs(
      title = glue("Alternative Regime Measures and Fragility Status: {outcome_label}"),
      subtitle = subtitle_text,
      x = "Coefficient on regime measure",
      y = NULL,
      shape = NULL,
      caption = caption_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      plot.caption = element_text(size = 8, color = "grey35", hjust = 0),
      axis.title.x = element_text(size = rel(0.85)),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(color = "grey30"),
      panel.grid.minor = element_blank()
    )
  
  p
}

# ------------------------------------------------------------------------------
# Plots ------------------------------------------------------------------------
# ------------------------------------------------------------------------------
coef_results <- coef_results |>
  mutate(
    plot = pmap(
      list(coef_df, outcome_label, include_gdp),
      ~ make_coef_plot(
        coef_df = ..1,
        outcome_label = ..2,
        include_gdp = ..3,
        x_limits = x_limits
      )
    )
  )

# ------------------------------------------------------------------------------
# Loop
# ------------------------------------------------------------------------------
coef_results$plot[[1]]

for (i in seq_len(nrow(coef_results))) {
  print(coef_results$plot[[i]])
}

# ------------------------------------------------------------------------------
# Saving the plots -------------------------------------------------------------
# ------------------------------------------------------------------------------
walk2(
  coef_results$plot,
  coef_results$file_stub,
  ~ ggsave(
    filename = paste0("coefplot_", .y, ".png"),
    plot = .x,
    width = 8.5,
    height = 6.2,
    dpi = 300
  )
)
