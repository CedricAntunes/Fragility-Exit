# ==============================================================================
# State Capacity-Fragility Pearson Correlations
# Two panels:
#   A. Net of year fixed effects
#   B. Net of country and year fixed effects
# ==============================================================================
# Author: Cedric Antunes (FGV-CEPESP) ------------------------------------------
# Date: May, 2026 --------------------------------------------------------------

# Assumes standardized_data_final is loaded

# Required packages ------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(fixest)
})

# ------------------------------------------------------------------------------
# Required variables -----------------------------------------------------------
# ------------------------------------------------------------------------------
required_vars <- c(
  "COUNTRY_NAME",
  "ISO_CODE_3",
  "YEAR",
  "STATE_CAPACITY",
  "NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y",
  "NEW_VDEM_NORL_LOADING_FACTOR_1_NORMALIZED",
  "WGI_LOADING_FACTOR_1_NORMALIZED",
  "WGI_NORL_LOADING_FACTOR_1_NORMALIZED",
  "VDEM_FRAGILE_IDEAL",
  "VDEM_NORL_FRAGILE_IDEAL",
  "WGI_FRAGILE_IDEAL",
  "WGI_NORL_FRAGILE_IDEAL"
)

# ------------------------------------------------------------------------------
# Helpers ----------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Re-oriating fragility: the higher, the more fragile
orient_higher_fragile <- function(x, fragile_binary) {
  cc <- complete.cases(x, fragile_binary)
  
  if (sum(cc) < 10) return(x)
  
  r <- suppressWarnings(
    cor(x[cc], fragile_binary[cc], use = "complete.obs")
  )
  
  if (is.na(r)) return(x)
  
  if (r < 0) -x else x
}

# FEs specifications 
residualize_fe <- function(data, var, fe_formula) {
  fml <- as.formula(paste0(var, " ~ 1 | ", fe_formula))
  
  resid(
    feols(
      fml,
      data = data,
      warn = FALSE,
      notes = FALSE
    )
  )
}

# ------------------------------------------------------------------------------
# Preparing the data -----------------------------------------------------------
# ------------------------------------------------------------------------------
df_plot <- standardized_data_final |>
  ungroup() |>
  transmute(
    country = COUNTRY_NAME,
    iso3    = as.character(ISO_CODE_3),
    year    = as.integer(YEAR),
    
    state_capacity = STATE_CAPACITY,
    
    vdem_raw      = .data[["NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y"]],
    vdem_norl_raw = .data[["NEW_VDEM_NORL_LOADING_FACTOR_1_NORMALIZED"]],
    wgi_raw       = .data[["WGI_LOADING_FACTOR_1_NORMALIZED"]],
    wgi_norl_raw  = .data[["WGI_NORL_LOADING_FACTOR_1_NORMALIZED"]],
    
    vdem_fragile_binary      = VDEM_FRAGILE_IDEAL,
    vdem_norl_fragile_binary = VDEM_NORL_FRAGILE_IDEAL,
    wgi_fragile_binary       = WGI_FRAGILE_IDEAL,
    wgi_norl_fragile_binary  = WGI_NORL_FRAGILE_IDEAL
  ) |>
  mutate(
    # Higher = more fragile
    vdem_fragility =
      orient_higher_fragile(vdem_raw, vdem_fragile_binary),
    
    vdem_fragility_norl =
      orient_higher_fragile(vdem_norl_raw, vdem_norl_fragile_binary),
    
    wgi_fragility =
      orient_higher_fragile(wgi_raw, wgi_fragile_binary),
    
    wgi_fragility_norl =
      orient_higher_fragile(wgi_norl_raw, wgi_norl_fragile_binary)
  )

# ------------------------------------------------------------------------------
# Long format ------------------------------------------------------------------
# ------------------------------------------------------------------------------
df_long <- df_plot |>
  pivot_longer(
    cols = c(
      vdem_fragility,
      vdem_fragility_norl,
      wgi_fragility,
      wgi_fragility_norl
    ),
    names_to = "outcome",
    values_to = "fragility"
  ) |>
  mutate(
    outcome = factor(
      outcome,
      levels = c(
        "vdem_fragility",
        "vdem_fragility_norl",
        "wgi_fragility",
        "wgi_fragility_norl"
      ),
      labels = c(
        "V-Dem fragility index",
        "V-Dem fragility index\n(excluding rule of law)",
        "WGI fragility index",
        "WGI fragility index\n(excluding rule of law)"
      )
    )
  )

# ------------------------------------------------------------------------------
# Computing residualized Pearson correlations ----------------------------------
# ------------------------------------------------------------------------------
compute_resid_cor <- function(data, fe_formula, adjustment_label) {
  
  data |>
    filter(
      !is.na(state_capacity),
      !is.na(fragility),
      !is.na(year),
      !is.na(iso3)
    ) |>
    group_by(outcome) |>
    group_modify(~{
      d <- .x |>
        select(iso3, year, state_capacity, fragility) |>
        drop_na()
      
      d <- d |>
        mutate(
          state_capacity_resid =
            residualize_fe(d, "state_capacity", fe_formula),
          
          fragility_resid =
            residualize_fe(d, "fragility", fe_formula)
        )
      
      ct <- cor.test(
        d$state_capacity_resid,
        d$fragility_resid
      )
      
      tibble(
        adjustment = adjustment_label,
        n = nrow(d),
        estimate = unname(ct$estimate),
        conf_low = ct$conf.int[1],
        conf_high = ct$conf.int[2]
      )
    }) |>
    ungroup()
}

# ------------------------------------------------------------------------------
# Computing both panels --------------------------------------------------------
# ------------------------------------------------------------------------------
cor_yearfe <- compute_resid_cor(
  data = df_long,
  fe_formula = "year",
  adjustment_label = "Net of year fixed effects"
)

cor_country_yearfe <- compute_resid_cor(
  data = df_long,
  fe_formula = "iso3 + year",
  adjustment_label = "Net of country and year fixed effects"
)

cor_results <- bind_rows(
  cor_yearfe,
  cor_country_yearfe
) |>
  mutate(
    adjustment = factor(
      adjustment,
      levels = c(
        "Net of year fixed effects",
        "Net of country and year fixed effects"
      ),
      labels = c(
        "Net of year fixed effects",
        "Net of country and year fixed effects"
      )
    ),
    label = sprintf("%.2f", estimate),
    label_y = if_else(
      estimate < 0,
      pmax(conf_low - 0.04, -0.98),
      pmin(conf_high + 0.04, 0.08)
    ),
    hjust_val = if_else(estimate < 0, 1, 0)
  )

print(cor_results)

write_csv(
  cor_results,
  "state_capacity_fragility_pearson_correlations_two_panel.csv"
)

# ------------------------------------------------------------------------------
# Plotting ---------------------------------------------------------------------
# ------------------------------------------------------------------------------
theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 10.5, colour = "grey25"),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 10, colour = "black"),
      axis.text.y = element_text(size = 10.5, colour = "black"),
      strip.text = element_text(face = "bold", size = 10.5),
      strip.background = element_rect(
        fill = "grey95",
        colour = "grey75",
        linewidth = 0.5
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(
        colour = "grey85",
        linewidth = 0.35
      ),
      panel.border = element_rect(
        fill = NA,
        colour = "grey70",
        linewidth = 0.5
      ),
      legend.position = "none",
      plot.margin = margin(10, 18, 10, 10)
    )
}

dark_red <- "#8B1A1A"

p_state_capacity_cor_two_panel <- ggplot(
  cor_results,
  aes(x = outcome, y = estimate)
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
    colour = dark_red
  ) +
  geom_point(
    size = 1.8,
    shape = 21,
    stroke = 1.1,
    fill = dark_red,
    colour = dark_red
  ) +
  geom_text(
    aes(
      y = label_y,
      label = label,
      hjust = hjust_val
    ),
    size = 3.6,
    colour = "grey20"
  ) +
  coord_flip(clip = "off") +
  facet_wrap(~ adjustment, ncol = 2) +
  scale_y_continuous(
    limits = c(-1, 0.10),
    breaks = seq(-1, 0, by = 0.2),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(
    x = NULL,
    y = "Pearson correlation",
    title = "State capacity and fragility",
    subtitle = "Pearson correlations after residualizing both variables; whiskers show 95% confidence intervals"
  ) +
  theme_pub()

p_state_capacity_cor_two_panel

# Saving the plot --------------------------------------------------------------
ggsave(
  filename = "state_capacity_fragility_pearson_correlations_two_panel.png",
  plot = p_state_capacity_cor_two_panel,
  width = 12,
  height = 6.5,
  dpi = 300
)
