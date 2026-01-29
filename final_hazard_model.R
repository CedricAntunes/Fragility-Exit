suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(splines)
  library(fixest)
  library(brglm2)
  library(glmmTMB)
  library(margins)
  library(sandwich)
  library(lmtest)
  library(stringr)
})

# ---------------- Settings ----------------
COUNTRY <- "COUNTRY_NAME"
YEAR    <- "YEAR"
STATUS  <- "VDEM_STATUS_IDEAL"   # expected: "Fragile","Transitioning","Non Fragile"

K_GRID            <- c(3L, 5L, 7L)
LAG_YEARS         <- 1L
COVERAGE_MIN      <- 0.70
DROP_TOP_LEVERAGE <- 0.02

# Full menu (used for coverage audits + optional expansion)
PREDICTORS <- c(
  "ODA_RECEIVED_PER_CAPITA",
  "GDP_GROWTH",
  "LOG_GDP_PER_CAPITA",
  "GDP_DEFLATOR",
  "POLITICAL_REGIME",
  "ELECTORAL_DEMOCRACY_SCORE",
  "LIBERAL_DEMOCRACY_SCORE",
  "TERRITORIAL_FRAGMENTATION",
  "INSTITUTIONAL_DEMOCRACY_SOCRE",      # keep literal if that's in your data
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

# Theory-first CORE base names (lag suffix appended automatically)
CORE_VARS     <- c("LIBERAL_DEMOCRACY_SCORE",
                   "REGIME_DURABILITY_YEARS",
                   "CONFLICT_INTENSITY_YEAR",
                   "AVG_CONFLICT_INTENSITY",
                   "GDP_GROWTH")

CORE_OPTIONAL <- c("GDP_DEFLATOR")

FAC_VARS <- c("POLITICAL_REGIME", "PARTIAL_DEMOCRACY_WITH_FACTIONALISM")
FRAG_SET <- c("Fragile", "Transitioning")
NF_LAB   <- "Non Fragile"

# ---------------- Helpers ----------------

# A) Ensure consecutive years are truly consecutive
complete_country_year <- function(df) {
  df %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    tidyr::complete(
      !!rlang::sym(YEAR) := seq(min(.data[[YEAR]], na.rm = TRUE),
                                max(.data[[YEAR]], na.rm = TRUE),
                                by = 1L)
    ) %>%
    ungroup() %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]])
}

# B) 1 at FIRST year t starting a sustained run of k Non-Fragile years
as_sustained_nf_start <- function(status_vec, k) {
  n <- length(status_vec)
  out <- integer(n)
  nf <- (status_vec == NF_LAB)
  nf[is.na(nf)] <- FALSE
  if (n < k) return(out)
  for (t in seq_len(n - k + 1L)) {
    if (nf[t] && all(nf[t:(t + k - 1L)])) {
      out[t] <- 1L
      break
    }
  }
  out
}

# C) Choose spline df that won’t crash on short spells
choose_spline_df <- function(x, max_df = 3L) {
  u <- length(unique(x[!is.na(x)]))
  max(1L, min(max_df, u - 1L))
}

# D) Cluster extraction for glm/brglm2
cluster_from_glm <- function(model, data, cluster_col) {
  ridx <- as.integer(rownames(model.frame(model)))
  data[[cluster_col]][ridx]
}

coeftest_cluster_country <- function(model, data) {
  cl <- cluster_from_glm(model, data, COUNTRY)
  vc <- sandwich::vcovCL(model, cluster = cl)
  lmtest::coeftest(model, vcov. = vc)
}

# E) Influence drop (top leverage)
influence_drop <- function(glm_model, drop_frac = 0.02) {
  h <- hatvalues(glm_model)
  cutoff <- stats::quantile(h, 1 - drop_frac, na.rm = TRUE)
  which(h >= cutoff)
}

# F) Coverage audit
coverage_audit <- function(haz, base_predictors, lag_years) {
  lag_suffix <- paste0("_L", lag_years)
  lag_cols <- paste0(intersect(base_predictors, names(haz)), lag_suffix)
  lag_cols <- intersect(lag_cols, names(haz))
  if (!length(lag_cols)) return(setNames(numeric(0), character(0)))
  sapply(lag_cols, function(v) mean(!is.na(haz[[v]])))
}

# G) Prune ONLY on vars used in the model (prevents killing events)
prune_for_model <- function(haz, model_vars) {
  model_vars <- intersect(model_vars, names(haz))
  haz2 <- haz %>% filter(!if_any(all_of(model_vars), is.na))

  # drop zero-variance vars (after pruning)
  nzv <- vapply(haz2[model_vars], function(x) {
    ux <- unique(x)
    length(ux) > 1
  }, logical(1))

  list(
    haz  = haz2,
    vars = model_vars[nzv],
    dropped_nzv = model_vars[!nzv]
  )
}

# H) BUILD HAZARD (FIXED!)
# Key fix: compute sustained-run start (needs future years) on FULL panel,
# and compute last_year_cty on FULL panel, not on risk-set-truncated data.
build_hazard <- function(df_full, predictors, k_sustain, lag_years = 1L) {

  # last observed year per country (STATUS non-missing)
  end_year <- df_full %>%
    group_by(.data[[COUNTRY]]) %>%
    summarize(
      last_year_cty = max(.data[[YEAR]][!is.na(.data[[STATUS]])], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(last_eligible = last_year_cty - (k_sustain - 1L))

  # compute sustained-run starts on FULL panel (so we can verify K future years)
  df1 <- df_full %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      nf_start    = as_sustained_nf_start(.data[[STATUS]], k = k_sustain),
      prev_status = dplyr::lag(.data[[STATUS]])
    ) %>%
    ungroup() %>%
    # true exit: NF-run start preceded by Fragile/Transitioning
    mutate(
      EXIT_SUSTAINED = if_else(
        nf_start == 1L &
          !is.na(prev_status) &
          prev_status %in% FRAG_SET,
        1L, 0L
      )
    ) %>%
    select(-nf_start, -prev_status) %>%
    left_join(end_year, by = setNames(COUNTRY, COUNTRY))

  # censor last (K-1) years of each country in the ESTIMATION sample
  # (but do NOT censor away the years needed to detect nf_start — already done)
  df1_est <- df1 %>%
    filter(.data[[YEAR]] <= last_eligible) %>%
    select(-last_year_cty, -last_eligible)

  # first exit year (in estimation window)
  df1_est <- df1_est %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      ever_exited = any(EXIT_SUSTAINED == 1L),
      first_exit_year = {
        yy <- .data[[YEAR]][EXIT_SUSTAINED == 1L]
        if (length(yy)) min(yy) else NA_integer_
      },
      # risk set: only fragile/transitioning strictly before first sustained exit year
      at_risk = if_else(
        .data[[STATUS]] %in% FRAG_SET &
          (is.na(first_exit_year) | .data[[YEAR]] < first_exit_year),
        1L, 0L
      )
    ) %>%
    ungroup()

  # keep risk-set rows + event row
  haz0 <- df1_est %>% filter(at_risk == 1L | EXIT_SUSTAINED == 1L)

  # drop countries never at risk and never event (purely always non-fragile)
  haz0 <- haz0 %>%
    group_by(.data[[COUNTRY]]) %>%
    filter(any(at_risk == 1L) | any(EXIT_SUSTAINED == 1L)) %>%
    ungroup()

  # lag predictors
  predictors <- intersect(predictors, names(haz0))
  if (!length(predictors)) stop("None of the listed predictors are in the data.")

  lag_suffix <- paste0("_L", lag_years)

  haz1 <- haz0 %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(across(all_of(predictors),
                  ~ dplyr::lag(.x, lag_years),
                  .names = paste0("{.col}", lag_suffix))) %>%
    ungroup()

  # refactor factor lags
  refactor_if_exists <- function(df, var) {
    if (var %in% names(df)) df[[var]] <- forcats::fct_drop(factor(df[[var]], exclude = NULL))
    df
  }
  haz1 <- refactor_if_exists(haz1, paste0("POLITICAL_REGIME", lag_suffix))
  haz1 <- refactor_if_exists(haz1, paste0("PARTIAL_DEMOCRACY_WITH_FACTIONALISM", lag_suffix))

  # spell clock: attach event row to previous at-risk spell
  haz1 <- haz1 %>%
    arrange(.data[[COUNTRY]], .data[[YEAR]]) %>%
    group_by(.data[[COUNTRY]]) %>%
    mutate(
      lag_at_risk  = dplyr::lag(at_risk, default = 0L),
      spell_start  = (at_risk == 1L & lag_at_risk == 0L),
      spell_id_tmp = cumsum(spell_start),
      spell_id     = if_else(at_risk == 1L, spell_id_tmp, dplyr::lag(spell_id_tmp)),
      spell_id     = replace_na(spell_id, 0L)
    ) %>%
    group_by(.data[[COUNTRY]], spell_id) %>%
    mutate(t_at_risk = row_number() - 1L) %>%
    ungroup() %>%
    select(-lag_at_risk, -spell_start, -spell_id_tmp)

  stopifnot(all(haz1$EXIT_SUSTAINED %in% c(0L, 1L, NA)))
  haz1
}

# ---------------- Load & prep ----------------
stopifnot(exists("standardized_data_final"))

df <- standardized_data_final %>%
  mutate(
    !!YEAR   := suppressWarnings(as.integer(.data[[YEAR]])),
    !!STATUS := as.character(.data[[STATUS]])
  )

# sanity check: status labels
cat("\n[Audit] Unique STATUS labels:\n")
print(sort(unique(df[[STATUS]])))

# warn if predictors missing
missing_preds <- setdiff(PREDICTORS, names(df))
if (length(missing_preds)) warning("Predictors missing in data (pre-lag): ",
                                   paste(missing_preds, collapse = ", "))

# factorize regime-like vars pre-lag
for (v in FAC_VARS) if (v %in% names(df)) df[[v]] <- factor(df[[v]], exclude = NULL)

# standardize numeric predictors pre-lag (global z-scores)
num_vars <- setdiff(PREDICTORS, FAC_VARS)
num_vars <- intersect(num_vars, names(df))
for (v in num_vars) if (is.numeric(df[[v]])) df[[v]] <- as.numeric(scale(df[[v]]))

cat("\n[Audit] Rows in original df:", nrow(df), "\n")

# complete panel
df_full <- complete_country_year(df)
cat("[Audit] Rows after country-year completion:", nrow(df_full), "\n")

# summarize censored always-risk countries (never Non-Fragile)
sum_cty <- df_full %>%
  group_by(.data[[COUNTRY]]) %>%
  summarize(
    ever_nf   = any(.data[[STATUS]] == NF_LAB, na.rm = TRUE),
    ever_risk = any(.data[[STATUS]] %in% FRAG_SET, na.rm = TRUE),
    always_risk = {
      s <- .data[[STATUS]]
      s <- s[!is.na(s)]
      length(s) > 0 && all(s %in% FRAG_SET)
    },
    .groups = "drop"
  )
cat("[Audit] Countries ever at risk:", sum(sum_cty$ever_risk), "\n")
cat("[Audit] Countries always Fragile/Transitioning (never Non-Fragile):", sum(sum_cty$always_risk), "\n")

# ---------------- Master loop ----------------
COLLECT <- NULL
MODELS  <- list()

lag_suffix <- paste0("_L", LAG_YEARS)
CORE_NAMES <- paste0(CORE_VARS, lag_suffix)
CORE_OPTIONAL_NAMES <- paste0(CORE_OPTIONAL, lag_suffix)

for (K in K_GRID) {
  cat("\n================= K =", K, "=================\n")

  haz_raw <- build_hazard(df_full, predictors = PREDICTORS, k_sustain = K, lag_years = LAG_YEARS)

  evt_raw <- sum(haz_raw$EXIT_SUSTAINED == 1L, na.rm = TRUE)
  cat("[Audit] Risk-set + event rows (post-censor, pre-prune):", nrow(haz_raw), "\n")
  cat("[Audit] Event count (pre-prune):", evt_raw, "\n")

  if (evt_raw > 0) {
    cat("[Audit] First few events:\n")
    print(haz_raw %>%
            filter(EXIT_SUSTAINED == 1L) %>%
            select(.data[[COUNTRY]], .data[[YEAR]], .data[[STATUS]], at_risk, t_at_risk) %>%
            head(15))
  }

  # coverage audit (informational)
  cov <- coverage_audit(haz_raw, PREDICTORS, LAG_YEARS)
  if (length(cov)) {
    cat("[Audit] Lagged non-missing (top 10):\n")
    print(sort(round(cov, 3), decreasing = TRUE)[1:min(10, length(cov))])
  }

  # eligible predictors by coverage
  keep_lag <- names(cov)[cov >= COVERAGE_MIN]

  # CORE selection (coverage-aware)
  core <- intersect(CORE_NAMES, keep_lag)
  core <- unique(c(core, intersect(CORE_OPTIONAL_NAMES, keep_lag)))

  # if core too thin, expand with highest-coverage eligible lags
  if (length(core) < 3L) {
    extras <- setdiff(names(sort(cov, decreasing = TRUE)), core)
    extras <- extras[extras %in% keep_lag]
    core <- unique(c(core, head(extras, 5L)))
  }
  cat("[Core] Candidate vars:", paste(core, collapse = ", "), "\n")

  # prune ONLY on variables used in the model (plus t_at_risk is always present)
  pr <- prune_for_model(haz_raw, core)
  haz <- pr$haz
  core_used <- pr$vars

  if (length(pr$dropped_nzv)) {
    cat("[Core] Dropped zero-variance after pruning:", paste(pr$dropped_nzv, collapse = ", "), "\n")
  }

  evt <- sum(haz$EXIT_SUSTAINED == 1L, na.rm = TRUE)
  cat("[Core] Rows after CORE-only pruning:", nrow(haz), "\n")
  cat("[Core] Event count after pruning:", evt, "\n")

  # If DV constant, skip cleanly
  if (evt == 0L || evt == nrow(haz)) {
    cat("Skipping K =", K, "because EXIT_SUSTAINED is constant after pruning.\n")
    next
  }

  # adaptive spline df
  df_spline <- choose_spline_df(haz$t_at_risk, max_df = 3L)
  cat("[Duration] Using spline df =", df_spline, "\n")

  # ---------------- Models ----------------

  # (1) Pooled logit + Year FE + duration spline (clustered by country)
  f_pool <- as.formula(paste0(
    "EXIT_SUSTAINED ~ ", paste(core_used, collapse = " + "),
    " + splines::ns(t_at_risk, ", df_spline, ") + i(", YEAR, ")"
  ))
  m_pool <- try(
    fixest::feglm(f_pool, data = haz, family = "logit", cluster = COUNTRY),
    silent = TRUE
  )
  cat("\n=== (1) Pooled logit + Year FE + duration spline (clustered) ===\n")
  if (inherits(m_pool, "try-error")) {
    cat("Pooled FE model failed (separation / too few events / FE cell issues).\n")
  } else {
    print(summary(m_pool))
  }

  # (1b) Influence check: GLM (no year FE) with duration spline
  f_glm <- as.formula(paste0(
    "EXIT_SUSTAINED ~ ", paste(core_used, collapse = " + "),
    " + splines::ns(t_at_risk, ", df_spline, ")"
  ))
  m_glm <- glm(f_glm, data = haz, family = binomial("logit"))

  idx_hi <- influence_drop(m_glm, drop_frac = DROP_TOP_LEVERAGE)
  if (length(idx_hi)) {
    haz_rb <- haz[-idx_hi, , drop = FALSE]
    evt_rb <- sum(haz_rb$EXIT_SUSTAINED == 1L, na.rm = TRUE)
    cat("\n--- Influence: dropped top ", round(100 * DROP_TOP_LEVERAGE, 1),
        "% leverage (n=", length(idx_hi), "); events now=", evt_rb, " ---\n", sep = "")
    m_pool_rb <- try(
      fixest::feglm(f_pool, data = haz_rb, family = "logit", cluster = COUNTRY),
      silent = TRUE
    )
    if (!inherits(m_pool_rb, "try-error")) print(summary(m_pool_rb))
  }

  # (2) Bias-reduced logit (brglm2) + clustered SEs (no year FE; more stable)
  m_brglm <- try(
    glm(f_glm, data = haz, family = binomial("logit"),
        method = "brglmFit", control = brglmControl(type = "AS_mixed")),
    silent = TRUE
  )
  cat("\n=== (2) Bias-reduced logit (brglm2) + clustered SEs ===\n")
  if (inherits(m_brglm, "try-error")) {
    cat("brglm2 failed (still possible with extreme separation / too few events).\n")
  } else {
    print(summary(m_brglm))
    cat("\n--- Clustered SEs (country) ---\n")
    print(coeftest_cluster_country(m_brglm, haz))

    # collect coef + clustered SE
    co <- coef(m_brglm)
    cl <- cluster_from_glm(m_brglm, haz, COUNTRY)
    V  <- sandwich::vcovCL(m_brglm, cluster = cl)
    se <- sqrt(diag(V))
    tmp <- data.frame(K = K, term = names(co), beta = unname(co), se = unname(se), row.names = NULL)
    COLLECT <- bind_rows(COLLECT, tmp)
  }

  # (3) Random intercept (country) — CORE
  cat("\n=== (3) Random-intercept logit (country) ===\n")
  f_re <- as.formula(paste0(
    "EXIT_SUSTAINED ~ ", paste(core_used, collapse = " + "),
    " + splines::ns(t_at_risk, ", df_spline, ") + (1 | ", COUNTRY, ")"
  ))
  m_re <- try(
    glmmTMB(f_re, data = haz, family = binomial(),
            control = glmmTMBControl(
              optimizer = optim, optArgs = list(method = "BFGS"),
              optCtrl = list(maxit = 1e4)
            )),
    silent = TRUE
  )
  if (inherits(m_re, "try-error")) {
    cat("glmmTMB failed to converge (sparse categories / quasi-separation likely).\n")
  } else {
    print(summary(m_re))
  }

  # (4) Absorbed FE logit (country + year) — robustness
  cat("\n=== (4) Absorbed FE logit (country + year) ===\n")
  f_abs <- as.formula(paste0(
    "EXIT_SUSTAINED ~ ", paste(core_used, collapse = " + "),
    " + splines::ns(t_at_risk, ", df_spline, ") | ", COUNTRY, " + ", YEAR
  ))
  m_abs <- try(
    fixest::feglm(f_abs, data = haz, family = "logit", cluster = COUNTRY),
    silent = TRUE
  )
  if (inherits(m_abs, "try-error")) {
    cat("Absorbed FE logit failed (separation/perfect prediction in FE cells likely).\n")
  } else {
    print(summary(m_abs))
  }

  # (5) AMEs from GLM (no year FE)
  cat("\n=== (5) Average Marginal Effects (GLM) ===\n")
  ame <- try(margins::margins(m_glm, data = haz), silent = TRUE)
  if (inherits(ame, "try-error")) {
    cat("AME computation failed (often due to spline/rare events). You can still report ORs.\n")
  } else {
    print(summary(ame))
  }

  # store models
  MODELS[[paste0("K", K, "_pool")]]  <- m_pool
  MODELS[[paste0("K", K, "_glm")]]   <- m_glm
  MODELS[[paste0("K", K, "_brglm")]] <- m_brglm
  MODELS[[paste0("K", K, "_re")]]    <- m_re
  MODELS[[paste0("K", K, "_abs")]]   <- m_abs

  cat("\n================= End K =", K, "=================\n")
}

cat("\n[Done] Coefs in `COLLECT`; models in `MODELS`.\n")
