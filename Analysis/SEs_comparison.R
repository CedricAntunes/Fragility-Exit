> # ------------------------------------------------------------
> # 0) Force fixest + plm to use EXACT same estimation sample
> # ------------------------------------------------------------
> vars_needed <- c("VDEM_FRAGILE_IDEAL", baseline_vars, "ISO_CODE_3", "YEAR")
> 
> df_fe <- standardized_data_final %>%
+     drop_na(all_of(vars_needed))
> 
> # ------------------------------------------------------------
> # 1) Two-way FE OLS (within)
> # ------------------------------------------------------------
> fe_mod <- feols(
+     form_baseline,
+     data  = df_fe,
+     fixef = c("ISO_CODE_3", "YEAR")
+ )
> 
> # SE variants from fixest
> se_cty <- fixest::se(fe_mod, cluster = ~ ISO_CODE_3)          # country-only cluster
> se_tw  <- fixest::se(fe_mod, cluster = ~ ISO_CODE_3 + YEAR)   # two-way cluster
> se_hc1 <- fixest::se(fe_mod, vcov = "HC1")                    # HC1
> 
> # ------------------------------------------------------------
> # 2) Driscoll–Kraay SEs from plm
> # ------------------------------------------------------------
> plm_fe <- plm(
+     form_baseline,
+     data   = df_fe,
+     model  = "within",
+     effect = "twoways",
+     index  = c("ISO_CODE_3","YEAR")
+ )
> 
> V_dk  <- vcovSCC(plm_fe, type = "HC1", maxlag = 3)
> se_dk <- sqrt(diag(V_dk))
> 
> # ------------------------------------------------------------
> # 3) Build table (same coefficients; different SEs/stars)
> # ------------------------------------------------------------
> b <- coef(fe_mod)
> terms <- names(b)
> 
> # Ensure all SE vectors are named + aligned
> if (is.null(names(se_cty))) names(se_cty) <- terms
> if (is.null(names(se_tw)))  names(se_tw)  <- terms
> if (is.null(names(se_hc1))) names(se_hc1) <- terms
> 
> se_cty <- se_cty[terms]
> se_tw  <- se_tw[terms]
> se_hc1 <- se_hc1[terms]
> se_dk  <- se_dk[terms]
> 
> get_star <- function(est, se) {
+     ok <- is.finite(est) & is.finite(se) & se > 0
+     p  <- rep(NA_real_, length(est))
+     z  <- abs(est[ok] / se[ok])
+     p[ok] <- 2 * pnorm(-z)
+     ifelse(is.na(p), "",
+            ifelse(p < .01, "***",
+                   ifelse(p < .05, "**",
+                          ifelse(p < .10, "*", ""))))
+ }
> 
> pretty_term <- function(x) {
+     x %>%
+         str_replace_all("_", " ") %>%
+         str_to_sentence()
+ }
> 
> tab <- tibble(
+     term     = terms,
+     estimate = unname(b),
+     se_cty   = unname(se_cty),
+     se_tw    = unname(se_tw),
+     se_hc1   = unname(se_hc1),
+     se_dk    = unname(se_dk)
+ ) %>%
+     mutate(
+         Predictor = pretty_term(term),
+         
+         est_cty = sprintf("%.3f%s", estimate, get_star(estimate, se_cty)),
+         est_tw  = sprintf("%.3f%s", estimate, get_star(estimate, se_tw)),
+         est_hc1 = sprintf("%.3f%s", estimate, get_star(estimate, se_hc1)),
+         est_dk  = sprintf("%.3f%s", estimate, get_star(estimate, se_dk)),
+         
+         se_cty = ifelse(is.finite(se_cty), sprintf("(%.3f)", se_cty), ""),
+         se_tw  = ifelse(is.finite(se_tw),  sprintf("(%.3f)", se_tw),  ""),
+         se_hc1 = ifelse(is.finite(se_hc1), sprintf("(%.3f)", se_hc1), ""),
+         se_dk  = ifelse(is.finite(se_dk),  sprintf("(%.3f)", se_dk),  "")
+     ) %>%
+     select(Predictor,
+            est_cty, se_cty,
+            est_tw,  se_tw,
+            est_hc1, se_hc1,
+            est_dk,  se_dk)
> 
> stopifnot(nrow(tab) > 0)
> 
> # ------------------------------------------------------------
> # 4) LaTeX output
> # ------------------------------------------------------------
> kbl_out <- kable(
+     tab,
+     format   = "latex",
+     booktabs = TRUE,
+     escape   = FALSE,
+     align    = c("l", rep("c", 8)),
+     col.names = c(
+         "Predictor",
+         "Country cluster", " ",
+         "Two-way cluster", " ",
+         "HC1", " ",
+         "Driscoll--Kraay", " "
+     ),
+     caption = "Two-way Fixed Effects OLS (within): Alternative Standard Errors",
+     label   = "tab:fe_robust_se"
+ ) %>%
+     add_header_above(c(" " = 1, "SE variant" = 8)) %>%
+     kable_styling(latex_options = c("hold_position","scale_down"), font_size = 8) %>%
+     column_spec(1, width = "6.5cm") %>%
+     footnote(
+         general = paste(
+             "All columns report the same two-way FE-OLS point estimates (country and year fixed effects); only the variance estimator differs.",
+             "Country cluster: SEs clustered by country.",
+             "Two-way cluster: SEs clustered by country and year.",
+             "HC1: heteroskedasticity-robust.",
+             "Driscoll--Kraay: vcovSCC (HC1), maxlag = 3.",
+             "Significance: *** p<0.01, ** p<0.05, * p<0.10."
+         ),
+         threeparttable = TRUE,
+         escape = FALSE
+     )
> 
> writeLines(as.character(kbl_out), "table_fe_se_variants_with_country_cluster.tex")
> cat(as.character(kbl_out))
