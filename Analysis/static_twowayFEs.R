# ------------------------------------------------------------------------------
# Setting the empirical model --------------------------------------------------
# ------------------------------------------------------------------------------
# Baseline predictors 
baseline_vars <- c(
  "OIL_RENTS",
  # "COMBINED_POLITY_SCORE",
  #"INSTITUTIONAL_AUTOCRACY_SCORE",
  # "INSTITUTIONAL_DEMOCRACY_SOCRE",
  "GDP_GROWTH", 
  # "MAX_CONFLICT_INTENSITY", 
  "N_TOTAL_TROOPS",
  # "TERRITORIAL_FRAGMENTATION",
  # "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",         
  #"LIBERAL_DEMOCRACY_SCORE",
  "BOIX_DEMOCRACY",
  "REGIME_DURABILITY_YEARS",
  "LOG_GDP_PER_CAPITA",
  # "POLITICAL_REGIME",
  # "N_WAR_FRONTS",
  "ODA_RECEIVED_PER_CAPITA",
  #"LOG_ODA_RECEIVED_PER_CAPITA",
  # "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  #"ELECTORAL_DEMOCRACY_SCORE",
  # "POLITICAL_COMPETITION_SCORE",
  "CONFLICT_INTENSITY_YEAR",
  #"EF_INDEX",
  # "AVG_CONFLICT_INTENSITY",
  "GDP_DEFLATOR",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM"
)

# Baseline empirical specification
form_baseline <- as.formula(
  paste("VDEM_NORL_FRAGILE_IDEAL ~", paste(baseline_vars, 
                                      collapse = " + "))
)

# ------------------------------------------------------------------------------
# OLS (cluster-robust at country-year) -----------------------------------------
# ------------------------------------------------------------------------------
linear_model <- feols(form_baseline, 
                      data = standardized_data_final)

clustered_se_linear <- se(linear_model, 
                          cluster = ~ ISO_CODE_3 + YEAR)

# ------------------------------------------------------------------------------
# Logit (cluster-robust at COUNTRY_NAME) ---------------------------------------
# ------------------------------------------------------------------------------
logit_model  <- feglm(form_baseline, 
                      family = binomial(), 
                      data = standardized_data_final)

clustered_se_logit <- se(logit_model, 
                         cluster = ~ ISO_CODE_3 + YEAR)

# ------------------------------------------------------------------------------
# Fixed-effects (within) + Driscoll–Kraay SEs ----------------------------------
# ------------------------------------------------------------------------------
fe_mod <- feols(form_baseline, 
                data = standardized_data_final,
                # Two-way FEs
                fixef = c("ISO_CODE_3",
                          "YEAR"))

# Two-way clustered SEs 
se_tw <- se(fe_mod, cluster = ~ ISO_CODE_3 + YEAR)

# HC1
se_hc1 <- se(fe_mod, vcov = "HC1")

# DK robustness ----------------------------------------------------------------
plm_fe <- plm(form_baseline, 
              data = standardized_data_final, 
              model = "within",
              effect = "twoways", 
              index = c("ISO_CODE_3",
                        "YEAR"))

dk_hc1 <- vcovSCC(plm_fe, 
                  type = "HC1", 
                  maxlag = 3)

se_dk  <- sqrt(diag(dk_hc1))

# ------------------------------------------------------------------------------
# Logit + two-way fixed effects (country & year) -------------------------------
# ------------------------------------------------------------------------------

logit_fe_model <- feglm(
  form_baseline,
  family = binomial(),
  data   = standardized_data_final,
  fixef = c("ISO_CODE_3",
            "YEAR"))

clustered_se_logit_fe <- se(logit_fe_model, cluster = ~ ISO_CODE_3 + YEAR)

# ------------------------------------------------------------------------------
# Fit statistics ---------------------------------------------------------------
# ------------------------------------------------------------------------------

R2_linear <- summary(linear_model)$r.squared

# McFadden PR2 (pooled logit) - your original approach
logit_null <- glm(VDEM_FRAGILE_IDEAL ~ 1,
                  family = binomial,
                  data   = standardized_data_final)
McFadden_R2 <- 1 - (as.numeric(logLik(logit_model)) /
                    as.numeric(logLik(logit_null)))

# McFadden-type PR2 for FE logit from fixest (recommended)
McFadden_R2_FE <- as.numeric(fixest::fitstat(logit_fe_model, "pr2"))

add_stats <- list(
  c("R-squared (OLS)",               sprintf("%.3f", R2_linear),      "",                         ""),
  c("Pseudo R-squared (Logit)",      "",                              sprintf("%.3f", McFadden_R2), ""),
  c("Pseudo R-squared (Logit + FE)", "",                              "",                         sprintf("%.3f", McFadden_R2_FE)),
  c("AIC",                           sprintf("%.2f", AIC(linear_model)), sprintf("%.2f", AIC(logit_model)), sprintf("%.2f", AIC(logit_fe_model))),
  c("BIC",                           sprintf("%.2f", BIC(linear_model)), sprintf("%.2f", BIC(logit_model)), sprintf("%.2f", BIC(logit_fe_model)))
)

# ------------------------------------------------------------------------------
# LaTeX baseline table ---------------------------------------------------------
# ------------------------------------------------------------------------------
fixest::etable(
  list(
    "OLS (LPM)"        = linear_model,
    "Logit"            = logit_model,
    "FE (cty & year)"  = fe_mod,
    "Logit + FE"       = logit_fe_model
  ),
  vcov = list(
    ~ ISO_CODE_3 + YEAR,  # OLS 2-way
    ~ ISO_CODE_3 + YEAR,  # pooled Logit 2-way
    ~ ISO_CODE_3   + YEAR,  # FE-OLS 2-way
    ~ ISO_CODE_3   + YEAR   # FE-Logit 2-way
  ),
  dict = c(
    INSTITUTIONAL_AUTOCRACY_SCORE = "Institutional autocracy score",
    OIL_RENTS = "Oil rents (% of GDP)",
    FH_DEMOCRACY = "Democracy",
    GDP_GROWTH = "GDP Growth",
    N_TOTAL_TROOPS = "N total troops involved",
    LOG_GDP_PER_CAPITA = "log(GDP per capita)",
    REGIME_DURABILITY_YEARS = "Regime durability (in years)",
    BOIX_DEMOCRACY = "Democracy",
    ELECTORAL_DEMOCRACY_SCORE = "Electoral democracy",
    LOG_ODA_RECEIVED_PER_CAPITA = "log(ODA received per capita)",
    CONFLICT_INTENSITY_YEAR = "Conflict intensity (year)",
    GDP_DEFLATOR = "GDP deflator",
    EF_INDEX = "Ethnic fractionalization",
    PARTIAL_DEMOCRACY_WITH_FACTIONALISM = "Partial democracy with factionalism"
  ),
  digits = 2,
  fixef_sizes = TRUE,
  family = TRUE,
  tex = TRUE,
  fontsize = "footnotesize",
  notes = "All models use the same unbalanced country-year panel (1971-2022). Columns (1) and (2) report pooled OLS and pooled logit estimates; column (3) reports a two-way FE linear model with country and year effects absorbed; column (4) reports a two-way FE logit with country and year fixed effects. Standard errors are two-way clustered by country and year. Driscoll-Kraay HAC SEs with a three-year bandwidth yield nearly identical inference and are shown in the appendix. All covariates are standardized; logit coefficients are log-odds. Significance: $^{***} p<0.01$, $^{**} p<0.05$, $^{*} p<0.1$.",
  title = "Linear, Logistic, and Fixed-Effects Estimates of Fragility Status (1971–2022)"
)
