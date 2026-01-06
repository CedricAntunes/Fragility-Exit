# ------------------------------------------------------------------------------
# Core fitter (with α-tuning) --------------------------------------------------
# ------------------------------------------------------------------------------
# Setting glmnet panel 
 fit_glmnet_panel <- function(
   data,
   predictors,
   outcome               = "VDEM_STATUS_IDEAL",
   country_col           = "country",
   year_col              = "year",
   # Sticking to a 80/20 time split
   holdout_year          = NULL,
   # Character vector of factor predictors 
   categorical_predictors = NULL,
   # Country and year FEs
   include_fixed_effects  = FALSE,       
   kfold_countries        = 5,
   # For LASSO and Ridge 
   alphas                 = c(1, .75, .5, .25, 0),
   # Re-scaling
   scale_with_train       = TRUE          
 ) {
 
   # Data preparation
   df <- data |>
     select(all_of(c(country_col, 
                     year_col, 
                     predictors, 
                     outcome))) |>
     drop_na()
 
   df[[outcome]] <- factor(df[[outcome]], levels = c("Non Fragile", "Fragile"))
 if (any(is.na(df[[outcome]]))) stop("Unexpected outcome labels in ", outcome)
 if (nlevels(df[[outcome]]) != 2L) stop("Outcome must be binary (2 levels).")
 
 lv  <- levels(df[[outcome]])
 pos <- "Fragile"
 
   if (!is.null(categorical_predictors)) {
     for (v in intersect(categorical_predictors, names(df))) {
       df[[v]] <- factor(df[[v]])
     }
   }
 
   # time split (forward)
   if (is.null(holdout_year)) {
     yrs <- sort(unique(df[[year_col]]))
     holdout_year <- yrs[ceiling(0.8 * length(yrs))] 
   }
   train <- filter(df, .data[[year_col]] <= holdout_year)
   test  <- filter(df, .data[[year_col]] >  holdout_year)
   stopifnot(nrow(train) > 0L, nrow(test) > 0L)
 
   # formula (optionally with FE)
   rhs_vars <- setdiff(names(df), c(outcome, country_col, year_col))
   fe_term <- if (include_fixed_effects) {
     paste0(" + factor(", country_col, ") + factor(", year_col, ")")
   } else ""
   form <- as.formula(paste0(outcome, " ~ ", paste(rhs_vars, collapse = " + "), fe_term, " - 1"))
 
   # optional train-only scaling (prevents leakage even if upstream z-scoring was global)
   num_vars <- rhs_vars[sapply(train[, rhs_vars, drop = FALSE], is.numeric)]
   scaling <- list(enabled = scale_with_train, mu = NULL, sd = NULL, num_vars = num_vars)
 
   if (scale_with_train && length(num_vars) > 0) {
     mu <- sapply(train[, num_vars, drop = FALSE], mean)
     sd <- sapply(train[, num_vars, drop = FALSE], sd); sd[sd == 0] <- 1
     scale_with <- function(d) {
       d2 <- d
       d2[, num_vars] <- sweep(sweep(d2[, num_vars, drop = FALSE], 2, mu, "-"), 2, sd, "/")
       d2
     }
     train <- scale_with(train)
     test  <- scale_with(test)
     scaling$mu <- mu; scaling$sd <- sd
   }
 
   # matrices
   # matrices (build on TRAIN only; then align TEST)
   train_mm <- train
   test_mm  <- test
 
   # Drop unused factor levels in TRAIN so we don't create "future-level" columns
   fact_vars <- rhs_vars[sapply(train_mm[, rhs_vars, drop = FALSE], is.factor)]
   if (include_fixed_effects) fact_vars <- unique(c(fact_vars, country_col, year_col))
   if (length(fact_vars) > 0) {
     train_mm[, fact_vars] <- lapply(train_mm[, fact_vars, drop = FALSE], droplevels)
   }
 
   X_train <- model.matrix(form, data = train_mm)
   X_test  <- model.matrix(form, data = test_mm)
 
   # Align columns: add missing columns as zeros; drop extras
   miss <- setdiff(colnames(X_train), colnames(X_test))
   if (length(miss)) {
     X_test <- cbind(
       X_test,
       matrix(0, nrow(X_test), length(miss), dimnames = list(NULL, miss))
     )
   }
   X_test <- X_test[, colnames(X_train), drop = FALSE]
   X_cols <- colnames(X_train)
 
   y_train <- droplevels(train[[outcome]])
   y_test  <- droplevels(test[[outcome]])
   #lv <- levels(y_train); pos <- lv[2]
 
   # country-wise folds
   set.seed(42)
   K <- kfold_countries
   countries_tr <- unique(train[[country_col]])
   foldmap <- setNames(sample(rep(1:K, length.out = length(countries_tr))), countries_tr)
   foldid  <- as.integer(foldmap[train[[country_col]]])
 
   # class weights (minority upweighted)
   n_pos <- sum(y_train == pos); n_neg <- sum(y_train != pos)
   #w <- ifelse(y_train == pos, n_neg / max(1, n_pos), 1)
   w <- rep(1, length(y_train))
 
   # alpha tuning (AUC CV)
   fits <- lapply(alphas, function(a)
     cv.glmnet(
       x = X_train, y = y_train,
       family = "binomial",
       alpha = a,
       type.measure = "deviance",
       standardize = FALSE,               
       foldid = foldid,
       nlambda = 200,
       lambda.min.ratio = 1e-5,
       weights = w,
       keep = TRUE
     )
   )
   # <- sapply(fits, function(f) max(f$cvm))
   #best_i <- which.max(cv_auc)
   cv_loss <- sapply(fits, function(f) min(f$cvm))
   best_i  <- which.min(cv_loss)
   best_alpha <- alphas[best_i]
   best_fit   <- fits[[best_i]]
   best_alpha <- alphas[best_i]
   best_fit <- fits[[best_i]]
   # s_use <- best_fit$lambda.1se
   s_use <- best_fit$lambda.min
   
   # --- OOF predictions for calibration (robust) ---
   p_train_oof <- rep(NA_real_, length(y_train))
 
   for (k in sort(unique(foldid))) {
   idx_val <- which(foldid == k)
   idx_tr  <- which(foldid != k)
 
   fit_k <- glmnet(
     x = X_train[idx_tr, , drop = FALSE],
     y = y_train[idx_tr],
     family = "binomial",
     alpha = best_alpha,
     lambda = s_use,
     standardize = FALSE,
     weights = w[idx_tr]
   )
 
   p_train_oof[idx_val] <- as.numeric(
     predict(fit_k, newx = X_train[idx_val, , drop = FALSE], type = "response")
   )
 }
 
 # safety
 stopifnot(all(is.finite(p_train_oof)))
 
   # -- TRAIN thresholds (choose here; apply to TEST)
   # In-sample train probs (still useful to keep)
   p_train <- as.numeric(predict(best_fit, s = s_use, newx = X_train, type = "response"))
 
   # Use out-of-fold probs for threshold selection (less optimistic)
   p_thr <- as.numeric(p_train_oof)
 
   roc_tr  <- roc(response = y_train, predictor = p_thr, levels = lv, quiet = TRUE)
   thr_youden_tr <- as.numeric(
     coords(roc_tr, "best", best.method = "youden", transpose = TRUE)["threshold"]
   )
 
   grid <- seq(0.05, 0.95, by = 0.01)
   f1_train <- function(th){
     pr <- factor(ifelse(p_thr >= th, pos, lv[1]), levels = lv)
     tp <- sum(pr==pos & y_train==pos); fp <- sum(pr==pos & y_train!=pos); fn <- sum(pr!=pos & y_train==pos)
     prec <- ifelse(tp+fp==0, 0, tp/(tp+fp)); rec <- ifelse(tp+fn==0, 0, tp/(tp+fn))
     ifelse(prec+rec==0, 0, 2*prec*rec/(prec+rec))
   }
   thr_f1_tr <- grid[which.max(vapply(grid, f1_train, numeric(1)))]
 
 
   # TEST metrics
   p_test <- as.numeric(predict(best_fit, s = s_use, newx = X_test, type = "response"))
   roc_te <- roc(response = y_test, predictor = p_test, levels = lv, quiet = TRUE)
   auc_te <- as.numeric(auc(roc_te))
   ci_te  <- as.numeric(ci.auc(roc_te)) 
   
   # Out-of-fold train probabilities at the chosen lambda (for calibration)
   #idx_lambda <- which.min(abs(best_fit$lambda - s_use))
   #p_train_oof <- as.numeric(best_fit$fit.preval[, idx_lambda])
 
   # thresholds learned on TRAIN, applied to TEST
   pred_y <- factor(ifelse(p_test >= thr_youden_tr, pos, lv[1]), levels = lv)
   pred_f <- factor(ifelse(p_test >= thr_f1_tr,     pos, lv[1]), levels = lv)
   cm_y   <- confusionMatrix(pred_y, y_test, positive = pos)
   cm_f   <- confusionMatrix(pred_f, y_test, positive = pos)
 
   # PR-AUC + Brier
   y_bin <- as.integer(y_test == pos)
   pr    <- pr.curve(scores.class0 = p_test[y_bin==1],
                     scores.class1 = p_test[y_bin==0], curve = FALSE)
   pr_auc <- pr$auc.integral
   brier  <- mean((y_bin - p_test)^2)
   pi_test    <- mean(y_bin)
   brier_null <- pi_test * (1 - pi_test)
   bss        <- 1 - brier / brier_null
 
   # simple calibration (test)
   eps_cal <- 1e-6
   p_test_c <- pmin(pmax(p_test, eps_cal), 1 - eps_cal)
 
   cal_fit <- glm((y_test == pos) ~ qlogis(p_test_c), family = binomial())
   cal_intercept <- unname(coef(cal_fit)[1])
   cal_slope     <- unname(coef(cal_fit)[2])
 
   # coefficients at s_use
   nz <- coef(best_fit, s = s_use)
   nz_df <- data.frame(feature = rownames(nz), beta = as.numeric(nz), row.names = NULL) |>
     filter(beta != 0)
 
   list(
     # model + design metadata
     model = best_fit,
     best_alpha = best_alpha,
     lambda = s_use,
     model_formula = form,
     X_cols = X_cols,
     scaling = scaling,
     levels = lv, pos = pos,
     include_fixed_effects = include_fixed_effects,
 
     # CV + Test metrics
     #cv_auc_by_alpha = data.frame(alpha = alphas, cv_auc = cv_auc),
     cv_loss_by_alpha = data.frame(alpha = alphas, cv_loss = cv_loss),
     auc_test = auc_te,
     auc_test_ci = ci_te,
     pr_auc_test = pr_auc,
     brier_test  = brier,
     brier_null_test   = brier_null,
     brier_skill_test  = bss,
     pi_test           = pi_test,
     calib_intercept = cal_intercept,
     calib_slope     = cal_slope,
     
     p_train      = p_train,
     p_train_oof  = p_train_oof,
     p_test       = p_test,
     y_train      = y_train,
     y_test       = y_test,
     holdout_year = holdout_year,
 
     # thresholds (chosen on TRAIN)
     threshold_youden_train = thr_youden_tr,
     threshold_f1_train     = thr_f1_tr,
 
     # confusion matrices (TEST)
     confusion_youden = cm_y,
     confusion_f1     = cm_f,
 
     # extras
     nonzero_coefs = nz_df,
     X_train_dim = dim(X_train),
     X_test_dim  = dim(X_test)
   )
 }
 
 # ------------------------------------------------------------------------------
 # Helper: predict/evaluate on NEW data block -----------------------------------
 # ------------------------------------------------------------------------------
 eval_on_new <- function(fit, newdata, outcome, country_col, year_col) {
   df <- newdata |>
     select(all_of(c(country_col, year_col, outcome, all.vars(update(fit$model_formula, . ~ .)))))
   df <- drop_na(df)
 
   # applying train scaling if used
   if (isTRUE(fit$scaling$enabled) && length(fit$scaling$num_vars) > 0) {
     nv <- intersect(fit$scaling$num_vars, names(df))
     if (length(nv) > 0) {
       df[, nv] <- sweep(sweep(df[, nv, drop = FALSE], 2, fit$scaling$mu[nv], "-"),
                         2, fit$scaling$sd[nv], "/")
     }
   }
 
   # build design, then align cols to training X
   X_new <- model.matrix(fit$model_formula, data = df)
   miss  <- setdiff(fit$X_cols, colnames(X_new))
   if (length(miss)) X_new <- cbind(X_new, matrix(0, nrow(X_new), length(miss), dimnames = list(NULL, miss)))
   X_new <- X_new[, fit$X_cols, drop = FALSE]
 
   y_new <- factor(df[[outcome]], levels = fit$levels)
   p_new <- as.numeric(predict(fit$model, s = fit$lambda, newx = X_new, type = "response"))
 
   roc_new <- roc(response = y_new, predictor = p_new, levels = fit$levels, quiet = TRUE)
   auc_new <- as.numeric(auc(roc_new))
 
   list(n = length(p_new), auc = auc_new, probs = p_new, y = y_new)
 }
 
# ------------------------------------------------------------------------------
# Plotting helpers -------------------------------------------------------------
# ------------------------------------------------------------------------------
 plot_cv_and_paths <- function(fit, main_prefix = "Best model") {
   op <- par(mfrow = c(1,2), mar = c(5,4,3,1)); on.exit(par(op), add = TRUE)
   plot(fit$model, main = paste0(main_prefix, " — CV (Deviance), alpha=", fit$best_alpha))
   abline(v = log(fit$lambda), lty = 2)
   plot(fit$model$glmnet.fit, xvar = "lambda", main = "Coefficient paths")
   abline(v = log(fit$lambda), lty = 2)
 }
 
 # ------------------------------------------------------------------------------
 # Time-cut sensitivity (forward) -----------------------------------------------
 # ------------------------------------------------------------------------------
 time_cut_sensitivity <- function(data, predictors, outcome,
                                  country_col, year_col,
                                  include_fixed_effects = FALSE,
                                  scale_with_train = TRUE,
                                  cuts = NULL, min_test_n = 100) {
   df_clean <- data |>
     select(all_of(c(country_col, year_col, predictors, outcome))) |>
     drop_na()
   yrs <- sort(unique(df_clean[[year_col]]))
   if (is.null(cuts)) cuts <- tail(yrs[yrs < max(yrs)], 6)
 
   out <- lapply(cuts, function(cut) {
     fit <- try(fit_glmnet_panel(
       data, predictors, outcome, country_col, year_col,
       holdout_year = cut,
       include_fixed_effects = include_fixed_effects,
       scale_with_train = scale_with_train
     ), silent = TRUE)
     if (inherits(fit, "try-error")) return(NULL)
     if (fit$X_test_dim[1] < min_test_n) return(NULL)
     data.frame(cut = cut,
                n_train = fit$X_train_dim[1],
                n_test  = fit$X_test_dim[1],
                auc     = fit$auc_test,
                pr_auc  = fit$pr_auc_test,
                brier   = fit$brier_test)
   })
   do.call(rbind, Filter(Negate(is.null), out))
 }
 
 # ------------------------------------------------------------------------------
 # Region-out (LORO) evaluation -------------------------------------------------
 # ------------------------------------------------------------------------------
 region_out <- function(data, predictors, outcome,
                        country_col, year_col, region_col,
                        holdout_year, include_fixed_effects = FALSE,
                        scale_with_train = TRUE, min_test_n = 150) {
   regs <- sort(unique(data[[region_col]]))
   res <- lapply(regs, function(R) {
     tr <- subset(data, .data[[region_col]] != R)
     te <- subset(data, .data[[region_col]] == R & .data[[year_col]] > holdout_year)
     if (nrow(te) < min_test_n) return(NULL)
 
     fit <- fit_glmnet_panel(
       tr, predictors, outcome, country_col, year_col,
       holdout_year = holdout_year,
       include_fixed_effects = include_fixed_effects,
       scale_with_train = scale_with_train
     )
     ev <- eval_on_new(fit, te, outcome, country_col, year_col)
     data.frame(region = R, n_test = ev$n, auc = ev$auc)
   })
   do.call(rbind, Filter(Negate(is.null), res))
 }
 
# ------------------------------------------------------------------------------
# Running the specification ----------------------------------------------------
# ------------------------------------------------------------------------------
 predictors_baseline <- c(
   "ODA_RECEIVED_PER_CAPITA",
   "GDP_GROWTH",
   "LOG_GDP_PER_CAPITA",
   "GDP_DEFLATOR",
   #"POLITICAL_REGIME",
   #"FH_DEMOCRACY",
   #"FH_AUTOCRACY",
   #"DD_DEMOCRACY",
   "ELECTORAL_DEMOCRACY_SCORE",
   "LIBERAL_DEMOCRACY_SCORE",
   "TERRITORIAL_FRAGMENTATION",
   #"INSTITUTIONAL_DEMOCRACY_SOCRE",
   "INSTITUTIONAL_AUTOCRACY_SCORE",
   #"COMBINED_POLITY_SCORE",
   "REGIME_DURABILITY_YEARS",
   "COMPETITIVENESS_PARTICIPATION",
   "INSTITUTIONAL_EXECUTIVE_RECRUTIMENT",
   "INSTITUTIONAL_PARTICIPATION",
   "EXECUTIVE_CONSTRAINT_SCORE",
   "POLITICAL_COMPETITION_SCORE",
   "PARTIAL_DEMOCRACY_WITH_FACTIONALISM",
   "CONFLICT_INTENSITY_YEAR",
   "CONFLICT_CUMULATIVE_INTENSITY_ACROSS_YEARS",
   "N_WAR_FRONTS",
   #"MAX_CONFLICT_INTENSITY",
   #"AVG_CONFLICT_INTENSITY",
   "N_TOTAL_TROOPS"
 )
 
 set.seed(42)
 
 baseline_fit <- fit_glmnet_panel(
   data = final_clean_percentiles_data_normalized_1,
   predictors = predictors_baseline,
   outcome = "VDEM_STATUS_IDEAL",
   country_col = "COUNTRY_NAME",
   year_col = "YEAR",
   holdout_year = 2016,          # or NULL to auto 80/20
   include_fixed_effects = FALSE,
   kfold_countries = 5,
   scale_with_train = TRUE       # prevents any scaling leakage
 )
 
 # Inspecting & plottting
 # baseline_fit$cv_auc_by_alpha
 baseline_fit$cv_loss_by_alpha
 baseline_fit$best_alpha
 baseline_fit$auc_test; baseline_fit$auc_test_ci
 baseline_fit$pr_auc_test; baseline_fit$brier_test
 baseline_fit$brier_null_test
 baseline_fit$brier_skill_test
 baseline_fit$pi_test
 baseline_fit$calib_intercept; baseline_fit$calib_slope
 baseline_fit$confusion_youden    
                                          
 head(baseline_fit$nonzero_coefs)
 plot_cv_and_paths(baseline_fit, "Baseline")
 
 # Time-cut stability
 time_cut_summary <- time_cut_sensitivity(
   data = final_clean_percentiles_data_normalized_1,
   predictors = predictors_baseline,
   outcome = "VDEM_STATUS_IDEAL",
   country_col = "COUNTRY_NAME",
   year_col = "YEAR",
   include_fixed_effects = FALSE,
   scale_with_train = TRUE
 )
 time_cut_summary
 
 
 # Getting all coefficients (including zero) at chosen lambda
 all_coefs <- as.matrix(coef(baseline_fit$model, s = baseline_fit$lambda))
 
 coef_table <- data.frame(
   feature = rownames(all_coefs),
   beta    = as.numeric(all_coefs)
 ) |>
   # Dropping intercepet
   filter(feature != "(Intercept)") |>
   mutate(shrunk_to_zero = beta == 0)
 
 # Inspect
 print(head(coef_table, 20))  # first 20 rows
 table(coef_table$shrunk_to_zero)

 # ------------------------------------------------------------------------------
 # Cross-validated AUC by alpha table -------------------------------------------
 # ------------------------------------------------------------------------------
 baseline_fit$cv_loss_by_alpha |>
   mutate(
     alpha   = round(alpha, 2),
     cv_loss = round(cv_loss, 3)
   )
 
 # ------------------------------------------------------------------------------
 # Time-cut summary table -------------------------------------------------------
 # ------------------------------------------------------------------------------
 time_cut_summary |>
   mutate(
     cut    = as.integer(cut),
     auc    = round(auc, 3),
     pr_auc = round(pr_auc, 3),
     brier  = round(brier, 3)
   )
 
 # ------------------------------------------------------------------------------
 # Confusion Matrix output ------------------------------------------------------
 # ------------------------------------------------------------------------------
 extract_cm_row <- function(cm, label) {
   # cm is a caret::confusionMatrix object
   tibble(
     Thresholding     = label,
     Accuracy         = unname(cm$overall["Accuracy"]),
     Sensitivity      = unname(cm$byClass["Sensitivity"]),
     Specificity      = unname(cm$byClass["Specificity"]),
     `Pos Pred Value` = unname(cm$byClass["Pos Pred Value"]),
     `Neg Pred Value` = unname(cm$byClass["Neg Pred Value"])
   )
 }
 
 tab_cm <- bind_rows(
   extract_cm_row(baseline_fit$confusion_youden, "Youden's J"),
   extract_cm_row(baseline_fit$confusion_f1,     "F1-optimal")
 ) |>
   mutate(across(-Thresholding, ~ round(.x, 3)))
 
 tab_cm
 
 # ------------------------------------------------------------------------------
 # Platt scalling ---------------------------------------------------------------
 # ------------------------------------------------------------------------------
 eps <- 1e-6
 clip <- function(p) pmin(pmax(p, eps), 1 - eps)
 
 pos <- baseline_fit$pos
 
 # Use OOF train probs for calibration (best). Fall back to in-sample if needed.
 p_tr <- if (!is.null(baseline_fit$p_train_oof)) baseline_fit$p_train_oof else baseline_fit$p_train
 
 y_tr <- as.integer(baseline_fit$y_train == pos)
 y_te <- as.integer(baseline_fit$y_test  == pos)
 
 p_tr <- clip(as.numeric(p_tr))
 p_te <- clip(as.numeric(baseline_fit$p_test))
 
 # Platt scaler fit on TRAIN (OOF probs)
 cal <- glm(y_tr ~ qlogis(p_tr), family = binomial())
 
 # Apply to TEST
 p_te_cal <- plogis(coef(cal)[1] + coef(cal)[2] * qlogis(p_te))
 
 # Brier + Skill (vs prevalence baseline) on TEST
 brier_raw <- mean((y_te - p_te)^2)
 brier_cal <- mean((y_te - p_te_cal)^2)
 
 pi_te <- mean(y_te)
 brier_null <- pi_te * (1 - pi_te)
 
 bss_raw <- 1 - brier_raw / brier_null
 bss_cal <- 1 - brier_cal / brier_null
 
 c(pi_test = pi_te,
   brier_null = brier_null,
   brier_raw = brier_raw, bss_raw = bss_raw,
   brier_cal = brier_cal, bss_cal = bss_cal)
 
 # ------------------------------------------------------------------------------
 # Isotonic regression ----------------------------------------------------------
 #-------------------------------------------------------------------------------
 eps <- 1e-6
 clip <- function(p) pmin(pmax(p, eps), 1 - eps)
 
 pos <- baseline_fit$pos
 p_tr <- clip(as.numeric(baseline_fit$p_train_oof))
 y_tr <- as.integer(baseline_fit$y_train == pos)
 
 p_te <- clip(as.numeric(baseline_fit$p_test))
 y_te <- as.integer(baseline_fit$y_test == pos)
 
 # Fit isotonic calibration map on (p_tr, y_tr)
 iso <- isoreg(p_tr, y_tr)
 
 # Predict calibrated probs by interpolating the step function
 p_te_iso <- pmin(pmax(approx(iso$x, iso$yf, xout = p_te, rule = 2)$y, 0), 1)

 brier_raw <- mean((y_te - p_te)^2)
 brier_iso <- mean((y_te - p_te_iso)^2)
 
 c(brier_raw = brier_raw, brier_iso = brier_iso) 
