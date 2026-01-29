```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ==========================================================
# Fragility Exit Hazard — WP-style core ID + robustness
# Data object: standardized_data_final
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)      # feglm (GLM) with clustering and FE
  library(brglm2)      # bias-reduced logistic regression
  library(glmmTMB)     # random-intercept logit
  library(margins)     # marginal effects
  library(sandwich)    # clustered vcov
  library(lmtest)      # coeftest with custom vcov
  library(forcats)
})

# ------------------------------------------------------------------------------
# Settings ---------------------------------------------------------------------
# ------------------------------------------------------------------------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
# "Fragile","Transitioning","Non Fragile"
STATUS  <- "VDEM_STATUS_IDEAL"     
# Lag coverage threshold
COVERAGE_MIN <- 0.70               
# sustained-exit year horizons for robustness
K_GRID       <- c(3L, 5L, 7L)      
# 2% leverage drop in influence check
DROP_TOP_LEVERAGE <- 0.02          

# Full predictors menu (we will prune by coverage and then select a CORE)
PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "LOG_GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "COMBINED_POLITY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

# Theory-first CORE spec (kept if coverage permits)
CORE_NAMES <- c(
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "GDP_GROWTH_L1",            
)

# -----------------
# Helpers
# -----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  # lag predictors by 1 year (keep NA for coverage audit)
  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")
  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 1L), .names = "{.col}_L1")) %>%
    ungroup()

  # Clean factor lags (avoid weird numeric-level labels / NA as level)
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")

  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  # helper to print clustered coeftest for glm/brglm2
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

# -----------------
# Load & prep
# -----------------
stopifnot(exists("standardized_data_final"))
df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# regime-like vars as factors (pre-lag)
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors (pre-lag)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# -----------------
# Run the full pipeline for each K in K_GRID
# -----------------
for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")
  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # Coverage audit (top 10)
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # ---- CORE model selection (based on theory & your early signal)
  core <- intersect(CORE_NAMES, keep_lag)
  # add GDP_DEFLATOR_L1 if available (proxy for nominal/ToT)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) {
    # fall back to a larger kept set if coverage excluded too much
    core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  }
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # -----------------
  # POOLED (year FE on by default), clustered by country
  # -----------------
  f_core_yearfe <- as.formula(paste("EXIT_SUSTAINED ~", paste(core, collapse = " + "), "+ i(", YEAR, ")", sep = ""))
  m_pool_core <- feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # Influence diag: drop top 2% leverage and refit pooled without year FE (GLM) for transparency
  f_core <- as.formula(paste("EXIT_SUSTAINED ~", paste(core, collapse = " + ")))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =", length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # -----------------
  # RARE-EVENTS (bias-reduced) + clustered SEs — CORE
  # -----------------
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))
  
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  
  if (!exists("COLLECT")) COLLECT <- tmp else COLLECT <- bind_rows(COLLECT, tmp)


  # -----------------
  # RANDOM-INTERCEPT (country) — CORE (graceful fallback)
  # -----------------
  f_re <- as.formula(paste("EXIT_SUSTAINED ~", paste(core, collapse = " + "), "+ (1 |", COUNTRY, ")"))
  m_re_core <- try(
    glmmTMB(f_re, data = haz, family = binomial(),
            control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))),
    silent = TRUE
  )
  cat("\n=== Random-intercept (country) — CORE ===\n")
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge from default starts (likely sparse categories / quasi-separation).\n")
  } else {
    print(summary(m_re_core))
  }

  # -----------------
  # AMEs for CORE pooled GLM (no year FE; AMEs are about X, not time dummies)
  # -----------------
  ame <- margins(m_glm_core)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  print(summary(ame))

  cat("\n================= End K =", K, "=================\n")
}


```

## Full Predictors 

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ================================================================
# Discrete-time hazard of sustained exit from fragility (WP-ready)
# ================================================================

# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "LOG_GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  #"POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  #"TERRITORIAL_FRAGMENTATION",
  #"INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  #"POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  #"CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  #"MAX_CONFLICT_INTENSITY",
  #"AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "LOG_GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  #"TERRITORIAL_FRAGMENTATION_L1",
  #"INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  #"POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  #"CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  #"MAX_CONFLICT_INTENSITY_L1",
  #"AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df |>
    arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df |>
    group_by(.data[[COUNTRY]]) |>
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) |>
    ungroup() |>
    group_by(.data[[COUNTRY]]) |>
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) |>
    ungroup()

  haz <- df |>
    filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz |>
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 1L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT <- bind_rows(COLLECT, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

```


```{r, hazard, echo = FALSE, message = FALSE, warning = FALSE, results = 'asis'}

stopifnot(exists("COLLECT"))

# Keep only the CORE variables you plotted in the models
core_terms <- c(
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "GDP_GROWTH_L1",
  "GDP_DEFLATOR_L1"
)

# Nice labels for rows
var_labs <- c(
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-1)",
  REGIME_DURABILITY_YEARS_L1   = "Regime durability, years (t-1)",
  CONFLICT_INTENSITY_YEAR_L1   = "Conflict intensity (acute, t-1)",
  AVG_CONFLICT_INTENSITY_L1    = "Conflict intensity (chronic avg., t-1)",
  GDP_GROWTH_L1                = "Real GDP per capita growth (t-1)",
  GDP_DEFLATOR_L1              = "GDP deflator / price tailwind (t-1)"
)

# Compute ORs and 95% clustered CIs
tab <- COLLECT %>%
  filter(term %in% core_terms) %>%
  mutate(
    OR  = exp(beta),
    LCL = exp(beta - 1.96 * se),
    UCL = exp(beta + 1.96 * se),
    cell = sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL),
    term_label = var_labs[term]
  ) %>%
  select(term_label, K, cell)

# If a core variable was dropped for some K (coverage/collinearity), mark as blank
all_rows <- tibble(term_label = unname(var_labs))
all_cols <- tibble(K = sort(unique(tab$K)))
tab_full <- tidyr::complete(tab, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# Wide layout: one column per K
tab_wide <- tab_full %>%
  tidyr::pivot_wider(names_from = K, values_from = cell, names_prefix = "K = ")

# Order rows like CORE
tab_wide <- tab_wide %>%
  mutate(row_order = match(term_label, unname(var_labs))) %>%
  arrange(row_order) %>%
  select(-row_order)

# --- 3) Render LaTeX (booktabs) and write to file ---

library(kableExtra)

latex_tbl <- kable(
  tab_wide,
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (CORE rare-events logit, clustered by country)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c")
) %>%
  kable_styling(latex_options = c("hold_position","striped")) %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3)) %>%
  footnote(
    general = "Entries are odds ratios with 95\\% clustered confidence intervals. Models are bias-reduced (brglm2, AS_mixed). Covariates standardized before lagging. Outcome: first sustained exit (Non-Fragile) at t persisting for K years. Sample includes risk set (Fragile/Transitioning) and the first exit year.",
    threeparttable = TRUE, escape = FALSE
  )

cat(latex_tbl)

```

## Final Table

```{r, echo = FALSE, message = FALSE}
# ============================================================
# OR table (ALL predictors) from COLLECT — vectorized + robust
# ============================================================
stopifnot(exists("COLLECT"))

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-1)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-1)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-1)",
  GDP_DEFLATOR_L1            = "GDP deflator / price tailwind (t-1)",
  POLITICAL_REGIME_L1        = "Political regime (t-1)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-1)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-1)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-1)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-1)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-1)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-1)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-1)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-1)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-1)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-1)",
  N_WAR_FRONTS_L1                 = "# war fronts (t-1)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-1)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-1)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-1)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## Two Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 2L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_two_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_two_year <- bind_rows(COLLECT_two_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_two_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_two_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_two_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_two_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_two_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-2)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-2)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-2)",
  GDP_DEFLATOR_L1            = "GDP deflator / price tailwind (t-2)",
  POLITICAL_REGIME_L1        = "Political regime (t-2)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-2)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-2)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-2)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-2)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-2)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-2)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-2)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-2)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-2)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-2)",
  N_WAR_FRONTS_L1                 = "# war fronts (t-2)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-2)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-2)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-2)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)
```

## Three-Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 3L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_three_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_three_year <- bind_rows(COLLECT_three_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_three_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_three_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_three_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_three_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_three_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-3)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-3)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-3)",
  GDP_DEFLATOR_L1            = "GDP deflator (t-3)",
  POLITICAL_REGIME_L1        = "Political regime (t-3)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-3)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-3)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-3)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-3)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-3)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-3)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-3)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-3)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-3)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-3)",
  N_WAR_FRONTS_L1                 = "N war fronts (t-3)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-3)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-3)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-2)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## Four-Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 4L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_four_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_four_year <- bind_rows(COLLECT_four_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_four_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_four_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_four_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_four_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_four_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-4)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-4)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-4)",
  GDP_DEFLATOR_L1            = "GDP deflator (t-4)",
  POLITICAL_REGIME_L1        = "Political regime (t-4)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-4)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-4)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-4)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-4)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-4)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-4)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-4)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-4)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-4)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-4)",
  N_WAR_FRONTS_L1                 = "N war fronts (t-4)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-4)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-4)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-4)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## Five-Year Lags

```{r, echo = FALSE, message = FALSE, warning = FALSE}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",     # keep literal if that's in your data
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE",
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_INTENSITY_YEAR",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "AVG_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "TERRITORIAL_FRAGMENTATION_L1",
  "INSTITUTIONAL_DEMOCRACY_SOCRE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_INTENSITY_YEAR_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "AVG_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 5L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_five_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_five_year <- bind_rows(COLLECT_five_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_five_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_five_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_five_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_five_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_five_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-4)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-4)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-4)",
  GDP_DEFLATOR_L1            = "GDP deflator (t-4)",
  POLITICAL_REGIME_L1        = "Political regime (t-4)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-4)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-4)",
  TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-4)",
  INSTITUTIONAL_DEMOCRACY_SOCRE_L1 = "Institutional democracy (t-4)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-4)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-4)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-4)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-4)",
  CONFLICT_INTENSITY_YEAR_L1      = "Conflict intensity (acute, t-4)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-4)",
  N_WAR_FRONTS_L1                 = "N war fronts (t-4)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-4)",
  AVG_CONFLICT_INTENSITY_L1       = "Conflict intensity (chronic avg., t-4)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-4)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)

```

## LASSO Predictors Only

```{r}
# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"  # "Fragile","Transitioning","Non Fragile"

COVERAGE_MIN       <- 0.70
K_GRID             <- c(3L, 5L, 7L)
DROP_TOP_LEVERAGE  <- 0.02

PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH", # Drop
  "GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME", # Drop
  "ELECTORAL_DEMOCRACY_SCORE", # Drop
  "LIBERAL_DEMOCRACY_SCORE",
  "INSTITUTIONAL_AUTOCRACY_SCORE",
  "REGIME_DURABILITY_YEARS",
  "POLITICAL_COMPETITION_SCORE", # Drop
  "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS", # Drop
  "N_WAR_FRONTS",
  "MAX_CONFLICT_INTENSITY",
  "N_TOTAL_TROOPS"
)

CORE_NAMES <- c(
  "ODA_RECEIVED_PER_CAPITA_L1",
  "GDP_GROWTH_L1",
  "GDP_PER_CAPITA_L1",
  "GDP_DEFLATOR_L1",
  "ELECTORAL_DEMOCRACY_SCORE_L1",
  "LIBERAL_DEMOCRACY_SCORE_L1",
  "INSTITUTIONAL_AUTOCRACY_SCORE_L1",
  "REGIME_DURABILITY_YEARS_L1",
  "POLITICAL_COMPETITION_SCORE_L1",
  "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1",
  "N_WAR_FRONTS_L1",
  "MAX_CONFLICT_INTENSITY_L1",
  "N_TOTAL_TROOPS_L1"
)

# ---------------- Helpers ----------------
as_exit_dummy <- function(vec, k = 5L) {
  n <- length(vec); exit <- integer(n); nf <- (vec == "Non Fragile")
  for (t in seq_len(n)) if (nf[t] && t + k - 1 <= n && all(nf[t:(t + k - 1)])) { exit[t] <- 1L; break }
  exit
}

build_hazard_raw <- function(df, predictors, k_sustain) {
  df <- df %>% arrange(.data[[COUNTRY]], .data[[YEAR]])
  df <- df %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(EXIT_SUSTAINED = as_exit_dummy(.data[[STATUS]], k = k_sustain)) %>%
    ungroup() %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = dplyr::if_else(
        ever_exited, min(.data[[YEAR]][EXIT_SUSTAINED == 1L], na.rm = TRUE), NA_integer_
      ),
      at_risk = dplyr::if_else(
        .data[[STATUS]] %in% c("Fragile","Transitioning") &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>% ungroup()

  haz <- df %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  predictors <- intersect(predictors, names(haz))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  haz <- haz %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors), ~ dplyr::lag(.x, 1L), .names = "{.col}_L1")) %>%
    ungroup()

  # clean factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz <- refactor_if_exists(haz, "POLITICAL_REGIME_L1")
  haz <- refactor_if_exists(haz, "PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1")
  haz
}

prune_by_coverage <- function(haz_raw, predictors, cov_min) {
  lag_cols <- paste0(intersect(predictors, names(haz_raw)), "_L1")
  lag_cols <- intersect(lag_cols, names(haz_raw))
  if (!length(lag_cols)) stop("After lagging, no candidate predictors remained.")

  coverage <- sapply(lag_cols, function(v) mean(!is.na(haz_raw[[v]])))
  keep_lag <- names(coverage)[coverage >= cov_min]
  if (!length(keep_lag)) stop(paste0("No lagged predictors meet coverage >= ", cov_min))

  # strict complete-case on kept lags
  haz <- haz_raw %>% filter(!if_any(all_of(keep_lag), is.na))

  # drop zero-variance among kept lags
  nzv <- vapply(haz[keep_lag], function(x) { ux <- unique(x[!is.na(x)]); length(ux) > 1 }, logical(1))
  list(haz = haz, keep_lag = keep_lag[nzv], coverage = sort(coverage, decreasing = TRUE))
}

fit_clustered <- function(model) {
  vc <- vcovCL(model, cluster = model$data[[COUNTRY]])
  coeftest(model, vcov = vc)
}

influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

diag_metrics <- function(y, p) {
  out <- list(
    brier = mean((y - p)^2)
  )
  if (have_pROC)  out$auc_roc <- as.numeric(pROC::roc(y, p, quiet = TRUE)$auc)
  if (have_PRROC) out$auc_pr  <- PRROC::pr.curve(scores.class0 = p[y == 1],
                                                 scores.class1 = p[y == 0])$auc.integral
  out
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity: warn if any named predictors missing (pre-lag)
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse=", "))

# factorize regime-like vars pre-lag
fac_candidates <- c("POLITICAL_REGIME","PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
for (v in fac_candidates) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, fac_candidates)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# ---------------- Master loop ----------------
COLLECT_two_year <- NULL
MODELS  <- list()

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard_raw(df, predictors = PREDICTORS, k_sustain = K)
  cat("[Audit] Risk-set + event rows (pre-lag-filter):", nrow(haz_raw), "\n")

  # coverage audit
  lag_cols_all <- paste0(intersect(PREDICTORS, names(haz_raw)), "_L1")
  lag_cols_all <- intersect(lag_cols_all, names(haz_raw))
  lag_cov <- sapply(lag_cols_all, function(v) mean(!is.na(haz_raw[[v]])))
  cat("[Audit] Lagged non-missing (top 10):\n")
  print(sort(round(lag_cov, 3), decreasing = TRUE)[1:min(10, length(lag_cov))])

  pruned <- prune_by_coverage(haz_raw, predictors = PREDICTORS, cov_min = COVERAGE_MIN)
  haz      <- pruned$haz
  keep_lag <- pruned$keep_lag
  cat("[Pruning] Kept lagged predictors (", length(keep_lag), ")\n", sep = "")
  cat(paste(keep_lag, collapse = ", "), "\n")
  cat("[Pruning] Rows after pruning:", nrow(haz), "\n")

  # --------- Add duration dependence (spell clock) ----------
  haz <- haz %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(new_spell = (lag(!!sym(STATUS), default = first(.data[[STATUS]])) == "Non Fragile") &
                       .data[[STATUS]] %in% c("Fragile","Transitioning"),
           spell_id  = cumsum(replace_na(new_spell, FALSE))) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup()

  # --------- End-of-panel trim to avoid immortal-time bias ----------
  max_year <- max(haz[[YEAR]], na.rm = TRUE)
  haz <- haz %>% filter(.data[[YEAR]] <= (max_year - (K - 1L)))
  cat("[Trim] Dropped end-of-panel years to avoid immortal-time bias. Remaining rows:", nrow(haz), "\n")

  # --------- CORE selection ----------
  core <- intersect(CORE_NAMES, keep_lag)
  if ("GDP_DEFLATOR_L1" %in% keep_lag) core <- union(core, "GDP_DEFLATOR_L1")
  if (length(core) < 3L) core <- unique(c(core, head(setdiff(keep_lag, core), 5L)))
  cat("[Core] Variables used:", paste(core, collapse = ", "), "\n")

  # --------- Models ----------
  # (1) Pooled logit w/ Year FE + ns(t_at_risk,3), clustered by country
  f_core_yearfe <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3) + i(", YEAR, ")", sep = ""
  ))
  m_pool_core <- fixest::feglm(f_core_yearfe, data = haz, family = "logit", cluster = COUNTRY)
  cat("\n=== Pooled logit + Year FE + duration spline (clustered) — CORE ===\n")
  print(summary(m_pool_core))

  # (1b) Influence refit (drop top leverage) — GLM without year FE (transparent influence check)
  f_core <- as.formula(paste(
  "EXIT_SUSTAINED ~", paste(core, collapse = " + "),
  "+ splines::ns(t_at_risk, 3)"
  ))
  m_glm_core <- glm(f_core, data = haz, family = binomial("logit"))
  idx_hi <- influence_drop(m_glm_core, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    m_pool_core_rb <- fixest::feglm(f_core_yearfe, data = haz_rb, family = "logit", cluster = COUNTRY)
    cat("\n--- Refit after dropping top", round(100*DROP_TOP_LEVERAGE,1), "% leverage (n drop =",
        length(idx_hi), ") ---\n")
    print(summary(m_pool_core_rb))
  }

  # (2) Rare-events (bias-reduced) + clustered SEs — CORE
  m_rare_core <- glm(f_core, data = haz, family = binomial("logit"),
                     method = "brglmFit", control = brglmControl(type = "AS_mixed"))
  cat("\n=== Rare-events (bias-reduced) — CORE ===\n")
  print(summary(m_rare_core))
  cat("\n--- Rare-events + clustered SEs (country) ---\n")
  print(fit_clustered(m_rare_core))

  # collect coef + clustered SE for rare-events model
  co  <- coef(m_rare_core)
  V   <- vcovCL(m_rare_core, cluster = m_rare_core$data[[COUNTRY]])
  se  <- sqrt(diag(V))
  tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
  COLLECT_two_year <- bind_rows(COLLECT_two_year, tmp)

  # (3) Random-intercept (country frailty) — CORE
  cat("\n=== Random-intercept (country) — CORE ===\n")
  m_re_core <- try(
    glmmTMB(f_core, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re_core, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re_core))
  }

  # (4) Conditional FE logit (country FE absorbed) — robustness
  cat("\n=== Conditional FE logit (absorbing country + year) — CORE ===\n")
  m_condfe <- try(
    fixest::feglm(
      as.formula(paste(
        "EXIT_SUSTAINED ~", paste(core, collapse=" + "),
        "+ ns(t_at_risk,3) |", COUNTRY, "+", YEAR)),
      data = haz, family = "logit", cluster = COUNTRY
    ),
    silent = TRUE
  )
  if (inherits(m_condfe, "try-error")) {
    cat("Conditional FE logit failed (possibly all-country separation in some years).\n")
  } else {
    print(summary(m_condfe))
  }

  # (5) AMEs (from plain GLM without year FE; AMEs are about X, not time dummies)
  cat("\n=== Average Marginal Effects — CORE GLM ===\n")
  ame <- margins::margins(m_glm_core, data = haz)  # <-- add data=
  print(summary(ame))

  # (6) Predictive diagnostics
  phat <- predict(m_glm_core, type = "response")
  y    <- haz$EXIT_SUSTAINED
  dm   <- diag_metrics(y, phat)
  cat("\n=== Diagnostics (GLM w/ duration) ===\n")
  print(dm)

  # store models
  MODELS[[paste0("K", K, "_poolFE")]]  <- m_pool_core
  MODELS[[paste0("K", K, "_rare")]]    <- m_rare_core
  MODELS[[paste0("K", K, "_glm")]]     <- m_glm_core
  MODELS[[paste0("K", K, "_frailty")]] <- m_re_core
  MODELS[[paste0("K", K, "_condFE")]]  <- m_condfe

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefficient collection available in `COLLECT`; models in `MODELS`.\n")

# Example: tidy coef table across K for a few key variables
if (exists("COLLECT_two_year")) {
  key_terms <- c("(Intercept)", CORE_NAMES[CORE_NAMES %in% COLLECT_two_year$term], "ns(t_at_risk, 3)1", "ns(t_at_risk, 3)2", "ns(t_at_risk, 3)3")
  print(
    COLLECT_two_year %>%
      filter(term %in% key_terms) %>%
      mutate(z = beta / se, p = 2*pnorm(abs(z), lower.tail = FALSE)) %>%
      arrange(term, K)
  )
}

# ------------------------------------------------------------------------------
# LATEX Table ------------------------------------------------------------------
# ------------------------------------------------------------------------------

stopifnot(exists("COLLECT_two_year"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(kableExtra); library(stringr)
})

# 1) Terms to include: drop intercept, spline basis, and year FE artifacts
omit_patterns <- c("^\\(Intercept\\)$", "splines::ns\\(", "^ns\\(", "^factor\\(", "^i\\(")
omit_re <- paste(omit_patterns, collapse="|")

terms_all <- COLLECT_two_year %>%
  mutate(term = as.character(term)) %>%
  filter(!str_detect(term, omit_re))

# 2) Pretty labels (override where you want; others auto-format)
var_labs_override <- c(
  ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-1)",
  GDP_GROWTH_L1              = "Real GDP per capita growth (t-1)",
  GDP_PER_CAPITA_L1          = "GDP per capita (t-1)",
  GDP_DEFLATOR_L1            = "GDP deflator / price tailwind (t-1)",
  POLITICAL_REGIME_L1        = "Political regime (t-1)",
  ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-1)",
  LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-1)",
  INSTITUTIONAL_AUTOCRACY_SCORE_L1 = "Institutional autocracy (t-1)",
  REGIME_DURABILITY_YEARS_L1      = "Regime durability, years (t-1)",
  POLITICAL_COMPETITION_SCORE_L1  = "Political competition (t-1)",
  PARTIAL_DEMOCRACY_WITH_FACTIONALISM_L1 = "Partial democracy w/ factionalism (t-1)",
  CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS_L1 = "Cumulative conflict intensity (t-1)",
  N_WAR_FRONTS_L1                 = "# war fronts (t-1)",
  MAX_CONFLICT_INTENSITY_L1       = "Max conflict intensity (t-1)",
  N_TOTAL_TROOPS_L1               = "Total troops (t-1)"
)

auto_label <- function(x) {
  x %>%
    str_replace("_L1$", " (t-1)") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

# Vectorized labeler: named lookup with [v], fill NAs with auto labels
make_label <- function(v) {
  v <- as.character(v)
  lbl <- unname(var_labs_override[v])      # returns named vector with NAs for misses
  miss <- is.na(lbl)
  lbl[miss] <- auto_label(v[miss])
  lbl
}

# 3) Format OR [LCL, UCL], bold if CI excludes 1
fmt_or_ci <- function(b, s) {
  OR  <- exp(b); LCL <- exp(b - 1.96*s); UCL <- exp(b + 1.96*s)
  cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
  if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
  cell
}

# 4) Build long table (K x term)
tab_long <- terms_all %>%
  mutate(
    term_label = make_label(term),
    cell = mapply(fmt_or_ci, beta, se)
  ) %>%
  select(term_label, K, cell)

# Ensure K columns exist in order 3/5/7 and blanks where missing
all_rows <- tibble(term_label = sort(unique(tab_long$term_label)))
all_cols <- tibble(K = c(3L, 5L, 7L))
tab_full <- tidyr::complete(tab_long, term_label = all_rows$term_label, K = all_cols$K, fill = list(cell = ""))

# 5) Optional buckets for section headers
bucket_of <- function(lbl) {
  lbl <- as.character(lbl)
  out <- rep("Other", length(lbl))
  out[stringr::str_detect(lbl, stringr::regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
  out[stringr::str_detect(lbl, stringr::regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
  out[stringr::str_detect(lbl, stringr::regex("democ|autoc|polity|regime|competition|recruitment|fragment", ignore_case = TRUE))] <- "Institutions & politics"
  factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
}

tab_wide <- tab_full %>%
  dplyr::mutate(
    K = factor(K, levels = c(3,5,7), labels = c("K = 3","K = 5","K = 7")),
    bucket = bucket_of(term_label)
  ) %>%
  tidyr::pivot_wider(names_from = K, values_from = cell) %>%
  dplyr::arrange(bucket, term_label)

# 6) Render LaTeX with grouped sections (robust start/end calc)
kb <- kable(
  tab_wide %>% select(term_label, `K = 3`, `K = 5`, `K = 7`),
  format = "latex",
  booktabs = TRUE,
  linesep = "",
  caption = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  col.names = c("Predictors", "K = 3", "K = 5", "K = 7"),
  escape = FALSE,
  align = c("l","c","c","c"),
  label = "tab:all_or"
) %>%
  kable_styling(latex_options = c("hold_position","striped","scale_down"),
                font_size = 8) %>%                                 # smaller font + auto-resize
  column_spec(1, width = "4.5cm") %>%                              # wrap predictor names
  column_spec(2, width = "3.2cm") %>%                              # narrow K columns
  column_spec(3, width = "3.2cm") %>%
  column_spec(4, width = "3.2cm") %>%
  add_header_above(c(" " = 1, "Sustain horizon (years)" = 3))

# Add group headers without rle()
bucket_levels <- c("Economy","Institutions & politics","Conflict & security","Other")
# make sure ordering is locked
tab_wide$bucket <- factor(tab_wide$bucket, levels = bucket_levels)

idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

for (lev in bucket_levels) {
  idx <- idx_by_bucket[[lev]]
  if (!is.null(idx) && length(idx)) {
    kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
  }
}

kb <- footnote(
  kb,
  general = "Entries are odds ratios with 95% clustered confidence intervals. Models are bias-reduced (brglm2, ASmixed). Covariates standardized before lagging. Spline of time-at-risk included; estimates correspond to the rare-events specification collected in COLLECT (no year FE). For the year-FE version, rebuild with COLLECTYFE and rerun. Outcome: first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. We exclude Executive Recruitment and Polity Score from the core specification because their coefficients are not well identified: Executive Recruitment exhibits complete or quasi-separation in the risk set, yielding boundary estimates and uninformative confidence intervals, while Polity Score is highly collinear with our V-Dem democracy measures. Both variables are reported in Appendix Table A.x for completeness.",
  threeparttable = TRUE, escape = FALSE
)

cat(kb)
```
