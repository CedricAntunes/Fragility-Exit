# ============================================================
# OR table from COLLECT (vectorized) -> LaTeX (kableExtra)
# ============================================================
library(dplyr)
library(tidyr)
library(stringr)
library(knitr)
library(kableExtra)

make_or_ci_table <- function(
  COLLECT,
  K_levels = c(3L, 5L, 7L),
  ci_level = 0.95,
  caption  = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs",
  label    = "tab:all_or",
  file_out = NULL
){

  stopifnot(is.data.frame(COLLECT))
  req <- c("term","beta","se","K")
  miss <- setdiff(req, names(COLLECT))
  if(length(miss)) stop("COLLECT is missing columns: ", paste(miss, collapse=", "))

  # ---- 1) omit nuisance terms (intercept, spline bases, FE artifacts)
  omit_patterns <- c(
    "^\\(Intercept\\)$",
    "splines::ns\\(",
    "^ns\\(",
    "^factor\\(",
    "^i\\(",
    "^YEAR::"      # <--- IMPORTANT (drops year FE terms like YEAR::1975)
  )
  omit_re <- paste(omit_patterns, collapse="|")

  terms_all <- COLLECT %>%
    mutate(
      term = as.character(term),
      K    = as.integer(K)
    ) %>%
    filter(!str_detect(term, omit_re))

  # ---- 2) pretty labels: override + auto fallback
  var_labs_override <- c(
    ODA_RECEIVED_PER_CAPITA_L1 = "ODA per capita (t-1)",
    GDP_GROWTH_L1              = "Real GDP per capita growth (t-1)",
    GDP_PER_CAPITA_L1          = "GDP per capita (t-1)",
    GDP_DEFLATOR_L1            = "GDP deflator / price tailwind (t-1)",
    POLITICAL_REGIME_L1        = "Political regime (t-1)",
    ELECTORAL_DEMOCRACY_SCORE_L1 = "Electoral democracy (t-1)",
    LIBERAL_DEMOCRACY_SCORE_L1   = "Liberal democracy (t-1)",
    TERRITORIAL_FRAGMENTATION_L1 = "Territorial fragmentation (t-1)",
    INSTITUTIONAL_DEMOCRACY_SCORE_L1 = "Institutional democracy (t-1)",  # <- fix typo if needed
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

  make_label <- function(v) {
    v <- as.character(v)
    lbl <- unname(var_labs_override[v])
    miss <- is.na(lbl)
    lbl[miss] <- auto_label(v[miss])
    lbl
  }

  # ---- 3) OR + CI formatter (bold if CI excludes 1)
  z <- qnorm(1 - (1 - ci_level)/2)

  fmt_or_ci <- function(b, s) {
    if (is.na(b) || is.na(s)) return("")
    OR  <- exp(b)
    LCL <- exp(b - z*s)
    UCL <- exp(b + z*s)
    cell <- sprintf("%.2f [%.2f, %.2f]", OR, LCL, UCL)
    if (!is.na(LCL) && !is.na(UCL) && (UCL < 1 || LCL > 1)) cell <- paste0("\\textbf{", cell, "}")
    cell
  }

  tab_long <- terms_all %>%
    mutate(
      term_label = make_label(term),
      cell       = mapply(fmt_or_ci, beta, se)
    ) %>%
    select(term_label, K, cell)

  # ---- 4) complete grid (all terms x all K)
  tab_full <- tidyr::complete(
    tab_long,
    term_label = sort(unique(tab_long$term_label)),
    K          = K_levels,
    fill = list(cell = "")
  )

  # ---- 5) bucket headers
  bucket_of <- function(lbl) {
    lbl <- as.character(lbl)
    out <- rep("Other", length(lbl))
    out[str_detect(lbl, regex("conflict|war|troop", ignore_case = TRUE))] <- "Conflict & security"
    out[str_detect(lbl, regex("gdp|price|oda|per capita|growth|deflator", ignore_case = TRUE))] <- "Economy"
    out[str_detect(lbl, regex("democ|autoc|polity|regime|competition|recruitment|fragment|institution", ignore_case = TRUE))] <- "Institutions & politics"
    factor(out, levels = c("Economy","Institutions & politics","Conflict & security","Other"))
  }

  tab_wide <- tab_full %>%
    mutate(
      K = factor(K, levels = K_levels, labels = paste0("K = ", K_levels)),
      bucket = bucket_of(term_label)
    ) %>%
    pivot_wider(names_from = K, values_from = cell) %>%
    arrange(bucket, term_label)

  # ---- 6) render LaTeX
  kb <- kable(
    tab_wide %>% select(term_label, starts_with("K = ")),
    format    = "latex",
    booktabs  = TRUE,
    linesep   = "",
    caption   = caption,
    col.names = c("Predictors", paste0("K = ", K_levels)),
    escape    = FALSE,
    align     = c("l", rep("c", length(K_levels))),
    label     = label
  ) %>%
    kable_styling(
      latex_options = c("hold_position","striped","scale_down"),
      font_size = 8
    ) %>%
    column_spec(1, width = "4.7cm") %>%
    add_header_above(c(" " = 1, "Sustain horizon (years)" = length(K_levels)))

  # pack rows by bucket (robust indices)
  bucket_levels <- levels(tab_wide$bucket)
  idx_by_bucket <- split(seq_len(nrow(tab_wide)), tab_wide$bucket)

  for (lev in bucket_levels) {
    idx <- idx_by_bucket[[lev]]
    if (!is.null(idx) && length(idx)) {
      kb <- pack_rows(kb, lev, start = min(idx), end = max(idx), indent = FALSE)
    }
  }

  kb <- footnote(
    kb,
    general = paste0(
      "Entries are odds ratios with ", round(ci_level*100), "\\% confidence intervals computed from clustered SEs. ",
      "Outcome is first sustained exit (Non-Fragile) at t persisting for K years; sample is the risk set plus first exit year. ",
      "Models correspond to the estimates stored in COLLECT; spline time-at-risk and fixed effects (if any) are omitted from the table."
    ),
    threeparttable = TRUE,
    escape = FALSE
  )

  if (!is.null(file_out)) writeLines(as.character(kb), file_out)
  return(kb)
}

# ---- Usage (your case)
stopifnot(exists("COLLECT"))
kb <- make_or_ci_table(
  COLLECT   = COLLECT,
  K_levels  = c(3L,5L,7L),
  ci_level  = 0.95,
  caption   = "Sustained Exit from Fragility: Odds Ratios with 95\\% Clustered CIs (ALL predictors, rare-events logit)",
  label     = "tab:all_or",
  file_out  = "table_all_or.tex"
)

cat(kb)
# In Overleaf: \input{table_all_or.tex}
