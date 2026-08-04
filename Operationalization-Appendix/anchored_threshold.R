# =============================================================================
# Appendix Z: regenerate the anchored-threshold robustness tables (Z1, Z2)
# -----------------------------------------------------------------------------
# Reuses the validated classification engine in anchored_robustness.R (the
# add_flags_* and find_exits functions, and the shared constants) and adds a
# two-measure layer that reproduces the numbers in Tables Z1 and Z2 for both
# the rule-of-law and the no-rule-of-law V-Dem indices.
#
# The rule-of-law relative run is checked against the published 36 before any
# table is built. The no-rule-of-law index has no separately published
# universe, so its baseline is the analogous set identified by the same annual
# rule (39 episodes).
#
# Outputs: appendix_Z_summary.csv, appendix_Z_movers.csv, and (optionally)
#          appendix_Z_tables_generated.tex with the two tabular bodies.
#
# Base R only; sourcing anchored_robustness.R does not execute its main().
# =============================================================================

source("anchored_robustness.R")   # engine: add_flags_relative/_fixed, find_exits, constants

RDS_PATH  <- "standardized_data_final.RDS"
OUT_DIR   <- "."
WRITE_TEX <- TRUE

SCORE_RL   <- "NEW_VDEM_LOADING_FACTOR_1_NORMALIZED.y"        # rule of law
SCORE_NORL <- "NEW_VDEM_NORL_LOADING_FACTOR_1_NORMALIZED"     # no rule of law

# --------------------------------------------------------------- load once ---
.e <- new.env(); .nm <- load(RDS_PATH, envir = .e); D <- get(.nm[1], envir = .e)

make_df <- function(score_col) {
  df <- data.frame(
    country = as.character(D[[COUNTRY_COL]]),
    year    = as.integer(as.character(D[[YEAR_COL]])),
    score   = as.numeric(D[[score_col]]),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$score), ]
  df[order(df$country, df$year), ]
}

# ------------------------------------------------- classify one measure ------
classify <- function(score_col) {
  df   <- make_df(score_col)
  rel  <- find_exits(add_flags_relative(df))
  p30p <- quantile(df$score, START_PERCENTILE, names = FALSE, type = 7)
  p35p <- quantile(df$score, EXIT_PERCENTILE,  names = FALSE, type = 7)
  bm   <- df$year >= BASE_PERIOD[1] & df$year <= BASE_PERIOD[2]
  p30b <- quantile(df$score[bm], START_PERCENTILE, names = FALSE, type = 7)
  p35b <- quantile(df$score[bm], EXIT_PERCENTILE,  names = FALSE, type = 7)
  pool <- find_exits(add_flags_fixed(df, p30p, p35p))
  base <- find_exits(add_flags_fixed(df, p30b, p35b))
  list(rel = rel, pool = pool, base = base,
       cut_pool = c(p30p, p35p), cut_base = c(p30b, p35b))
}

# baseline exits dated to the 1990s that survive a given anchor
n90_survive <- function(rel, anchor) {
  in90 <- names(rel)[rel >= 1990 & rel <= 1999]
  c(total = length(in90), survive = sum(in90 %in% names(anchor)))
}

RL   <- classify(SCORE_RL)
NORL <- classify(SCORE_NORL)

# --- reproduction check on the rule-of-law relative run ----------------------
stopifnot(length(RL$rel) == 36)
diffs <- vapply(union(names(RL$rel), names(PUBLISHED_36)), function(k)
  !identical(as.integer(RL$rel[k]), as.integer(PUBLISHED_36[k])), logical(1))
if (any(diffs)) stop("Rule-of-law relative run does not reproduce Table 1.")
cat("Reproduction check (rule of law): PASSED (36 episodes match Table 1)\n\n")

# ------------------------------------------------- Table Z1: summary ---------
summary_row <- function(label, m, anchor_name) {
  anc <- if (anchor_name == "Pooled") m$pool else m$base
  R <- names(m$rel); A <- names(anc)
  n90 <- n90_survive(m$rel, anc)
  data.frame(
    Index    = label,
    Anchor   = anchor_name,
    Baseline = length(R),
    Anchored = length(A),
    Retained = length(intersect(R, A)),
    Dropped  = length(setdiff(R, A)),
    Added    = length(setdiff(A, R)),
    y1990s   = sprintf("%d/%d", n90["survive"], n90["total"]),
    stringsAsFactors = FALSE
  )
}
tabZ1 <- rbind(
  summary_row("Rule of law",    RL,   "Pooled"),
  summary_row("Rule of law",    RL,   "Base-period"),
  summary_row("No rule of law", NORL, "Pooled"),
  summary_row("No rule of law", NORL, "Base-period")
)
cat("== Table Z1: episode counts by index and anchor ==\n")
print(tabZ1, row.names = FALSE)
cat(sprintf("\nCutoffs  RL pooled (%.3f, %.3f)  base (%.3f, %.3f)\n",
            RL$cut_pool[1], RL$cut_pool[2], RL$cut_base[1], RL$cut_base[2]))
cat(sprintf("         NORL pooled (%.3f, %.3f)  base (%.3f, %.3f)\n",
            NORL$cut_pool[1], NORL$cut_pool[2], NORL$cut_base[1], NORL$cut_base[2]))

# ------------------------------------------------- Table Z2: movers ----------
# Under the POOLED anchor: dropped baseline episodes (with baseline year) and
# added country-years (with anchored year).
# Display names so the tables match the appendix prose (data uses WB names).
DISPLAY <- c("Gambia, The" = "The Gambia",
             "Congo, Rep." = "Republic of the Congo",
             "Iran, Islamic Rep." = "Iran")
disp <- function(nm) ifelse(nm %in% names(DISPLAY), DISPLAY[nm], nm)
fmt_list <- function(nm, yr) paste(sprintf("%s (%d)", disp(nm), as.integer(yr)),
                                   collapse = ", ")
movers_row <- function(label, m, direction) {
  R <- names(m$rel); P <- names(m$pool)
  if (direction == "Dropped") {
    ctry <- sort(setdiff(R, P)); yr <- m$rel[ctry]
  } else {
    ctry <- sort(setdiff(P, R)); yr <- m$pool[ctry]
  }
  data.frame(direction = direction, index = label,
             episodes = fmt_list(ctry, yr), stringsAsFactors = FALSE)
}
tabZ2 <- rbind(
  movers_row("Rule of law",    RL,   "Dropped"),
  movers_row("No rule of law", NORL, "Dropped"),
  movers_row("Rule of law",    RL,   "Added"),
  movers_row("No rule of law", NORL, "Added")
)
cat("\n== Table Z2: movers under the pooled anchor ==\n")
for (i in seq_len(nrow(tabZ2)))
  cat(sprintf("[%s] %-14s %s\n", tabZ2$direction[i], tabZ2$index[i], tabZ2$episodes[i]))

# ------------------------------------------------- write outputs -------------
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
write.csv(tabZ1, file.path(OUT_DIR, "appendix_Z_summary.csv"), row.names = FALSE)
write.csv(tabZ2, file.path(OUT_DIR, "appendix_Z_movers.csv"),  row.names = FALSE)

# ------------------------------------------------- emit LaTeX bodies ---------
if (WRITE_TEX) {
  con <- file(file.path(OUT_DIR, "appendix_Z_tables_generated.tex"), open = "w")
  wl <- function(...) writeLines(sprintf(...), con)
  wl("%% Auto-generated by appendix_Z_tables.R -- requires \\usepackage{booktabs}")
  # Z1
  wl("\\begin{tabular}{llcccccc}")
  wl("\\toprule")
  wl("Index & Anchor & Baseline & Anchored & Retained & Dropped & Added & 1990s \\\\")
  wl("\\midrule")
  for (i in seq_len(nrow(tabZ1)))
    wl("%s & %s & %d & %d & %d & %d & %d & %s \\\\",
       tabZ1$Index[i], tabZ1$Anchor[i], tabZ1$Baseline[i], tabZ1$Anchored[i],
       tabZ1$Retained[i], tabZ1$Dropped[i], tabZ1$Added[i], tabZ1$y1990s[i])
  wl("\\bottomrule")
  wl("\\end{tabular}")
  wl("")
  # Z2
  wl("\\begin{tabular}{l p{10.2cm}}")
  wl("\\toprule")
  wl(" & Episodes \\\\")
  wl("\\midrule")
  wl("\\multicolumn{2}{l}{\\textit{Dropped from baseline}} \\\\")
  wl("Rule of law & %s \\\\",    tabZ2$episodes[tabZ2$direction=="Dropped" & tabZ2$index=="Rule of law"])
  wl("No rule of law & %s \\\\", tabZ2$episodes[tabZ2$direction=="Dropped" & tabZ2$index=="No rule of law"])
  wl("\\addlinespace")
  wl("\\multicolumn{2}{l}{\\textit{Added by the anchor}} \\\\")
  wl("Rule of law & %s \\\\",    tabZ2$episodes[tabZ2$direction=="Added" & tabZ2$index=="Rule of law"])
  wl("No rule of law & %s \\\\", tabZ2$episodes[tabZ2$direction=="Added" & tabZ2$index=="No rule of law"])
  wl("\\bottomrule")
  wl("\\end{tabular}")
  close(con)
  cat("\nWrote appendix_Z_summary.csv, appendix_Z_movers.csv, "
      , "appendix_Z_tables_generated.tex\n", sep = "")
}
