# Democracy-Fragility Pearson Correlations -------------------------------------
# Author: Cedric Antunes (FGV-CEPESP) ------------------------------------------
# Date: April, 2026 ------------------------------------------------------------

# Required packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ------------------------------------------------------------------------------
# Data preparation -------------------------------------------------------------
# ------------------------------------------------------------------------------
df_plot <- standardized_data_final |>
  ungroup() |>
  transmute(
    country                        = COUNTRY_NAME,
    year                           = YEAR,
    pca_raw                        = NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y,
    pca_raw_nrl                    = NEW_VDEM_NORL_LOADING_FACTOR_1_NORMALIZED,
    LIBERAL_DEMOCRACY_SCORE,
    ELECTORAL_DEMOCRACY_SCORE,
    INSTITUTIONAL_DEMOCRACY_SCORE
  ) |>
  filter(
    !is.na(pca_raw),
    !is.na(pca_raw_nrl)
  ) |>
  mutate(
    state_quality       = rescale(pca_raw, to = c(0, 1)),
    state_quality_nrl   = rescale(pca_raw_nrl, to = c(0, 1)),
    fragility_index     = 1 - state_quality,
    fragility_index_nrl = 1 - state_quality_nrl
  )

# ------------------------------------------------------------------------------
# Reshaping to long format -----------------------------------------------------
# ------------------------------------------------------------------------------
df_long <- df_plot |>
  pivot_longer(
    cols = c(
      LIBERAL_DEMOCRACY_SCORE,
      ELECTORAL_DEMOCRACY_SCORE,
      INSTITUTIONAL_DEMOCRACY_SCORE
    ),
    names_to  = "democracy_measure",
    values_to = "democracy_score"
  ) |>
  mutate(
    democracy_measure = case_when(
      democracy_measure == "LIBERAL_DEMOCRACY_SCORE" ~ "Liberal democracy (V-Dem)",
      democracy_measure == "ELECTORAL_DEMOCRACY_SCORE" ~ "Electoral democracy (V-Dem)",
      democracy_measure == "INSTITUTIONAL_DEMOCRACY_SCORE" ~ "Institutional democracy (Polity-V)"
    ),
    democracy_measure = factor(
      democracy_measure,
      levels = c(
        "Liberal democracy (V-Dem)",
        "Electoral democracy (V-Dem)",
        "Institutional democracy (Polity-V)"
      )
    )
  ) |>
  pivot_longer(
    cols = c(fragility_index, fragility_index_nrl),
    names_to  = "outcome",
    values_to = "fragility"
  ) |>
  mutate(
    outcome = factor(
      outcome,
      levels = c("fragility_index", "fragility_index_nrl"),
      labels = c(
        "Fragility index",
        "Fragility index (excluding rule of law)"
      )
    )
  )

# ------------------------------------------------------------------------------
# Helper: Year FEs -------------------------------------------------------------
# ------------------------------------------------------------------------------
residualize_year <- function(x, year) {
  resid(lm(x ~ factor(year)))
}

# ------------------------------------------------------------------------------
# Computing year-residualized correlations + 95% CIs ---------------------------
# ------------------------------------------------------------------------------
cor_results_yearfe <- df_long |>
  filter(
    !is.na(democracy_score),
    !is.na(fragility),
    !is.na(year)
  ) |>
  group_by(outcome, democracy_measure) |>
  group_modify(~{
    d <- .x |>
      select(year, democracy_score, fragility) |>
      drop_na()
    
    x_resid <- residualize_year(d$democracy_score, d$year)
    y_resid <- residualize_year(d$fragility, d$year)
    
    ct <- cor.test(x_resid, y_resid)
    
    tibble(
      n         = nrow(d),
      estimate  = unname(ct$estimate),
      conf_low  = ct$conf.int[1],
      conf_high = ct$conf.int[2]
    )
  }) |>
  ungroup() |>
  mutate(
    label_y = pmin(conf_high + 0.05, 0.04),
    label   = sprintf("%.2f", estimate)
  )

# ------------------------------------------------------------------------------
# Theme ------------------------------------------------------------------------
# ------------------------------------------------------------------------------
theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title.position   = "plot",
      plot.title            = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle         = element_text(hjust = 0.5, size = 10.5, colour = "grey25"),
      axis.title.x          = element_text(size = 11),
      axis.title.y          = element_blank(),
      axis.text.x           = element_text(size = 10, colour = "black"),
      axis.text.y           = element_text(size = 10.5, colour = "black"),
      strip.text            = element_text(face = "bold", size = 10.5),
      strip.background      = element_rect(fill = "grey95", colour = "grey75", linewidth = 0.5),
      panel.grid.minor      = element_blank(),
      panel.grid.major.y    = element_blank(),
      panel.grid.major.x    = element_line(colour = "grey85", linewidth = 0.35),
      panel.border          = element_rect(fill = NA, colour = "grey70", linewidth = 0.5),
      legend.position       = "none",
      plot.margin           = margin(10, 18, 10, 10)
    )
}

# ------------------------------------------------------------------------------
# Plot -------------------------------------------------------------------------
# ------------------------------------------------------------------------------
green_dark  <- "#6cbf84"
green_light <- "#A1D99B"

p_cor_yearfe <- ggplot(
  cor_results_yearfe,
  aes(x = democracy_measure, y = estimate)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.45
  ) +
  geom_errorbar(
    aes(ymin = conf_low, ymax = conf_high),
    width = 0.12,
    linewidth = 0.9,
    colour = green_dark
  ) +
  geom_point(
    size = 1.5,
    shape = 21,
    stroke = 1.1,
    fill = green_dark,
    colour = green_dark
  ) +
  geom_text(
    aes(y = label_y, label = label),
    hjust = 0,
    size = 3.6,
    colour = "grey20"
  ) +
  coord_flip(clip = "off") +
  facet_wrap(~ outcome, ncol = 2) +
  scale_y_continuous(
    limits = c(-1, 0.06),
    breaks = seq(-1, 0, by = 0.2),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Pearson correlation",
    title = "Democracy and fragility",
    subtitle = "Correlations after removing year fixed effects from both variables; whiskers show 95% confidence intervals"
  ) +
  theme_bw()

# Plot
p_cor_yearfe

# Saving the plot --------------------------------------------------------------
ggsave(
  filename = file.path("democracy_pearson_correlations.png"),
  plot     = p_cor_yearfe,
  width    = 11,
  height   = 7,
  dpi      = 300
)
