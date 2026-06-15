# ==============================================================================
# functions.R
# Structured utilities and estimators for transfer learning with censored survival
# outcomes and fixed-time logistic risk models.
# ==============================================================================

local_r_lib <- file.path(getwd(), "_r_lib")
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Please install it first.", pkg))
  }
}

# ==============================================================================
# 1. Basic helpers
# ==============================================================================

choose_nfolds_binomial <- function(y, max_folds = 5) {
  y <- as.integer(y)
  n1 <- sum(y == 1, na.rm = TRUE)
  n0 <- sum(y == 0, na.rm = TRUE)
  if (n1 < 2 || n0 < 2) return(2)
  max(3, min(max_folds, n1, n0))
}

align_beta_to_X <- function(beta, X_colnames) {
  beta_names <- names(beta)
  beta <- as.numeric(beta)
  if (!is.null(beta_names)) names(beta) <- beta_names
  beta[!is.finite(beta)] <- 0
  if (is.null(beta_names)) {
    if (length(beta) != length(X_colnames)) {
      stop("beta has no names and its length does not match design matrix columns.")
    }
    names(beta) <- X_colnames
    return(beta)
  }
  out <- rep(0, length(X_colnames)); names(out) <- X_colnames
  common <- intersect(beta_names, X_colnames)
  out[common] <- beta[common]
  if ("(Intercept)" %in% beta_names && "Intercept" %in% X_colnames) {
    out["Intercept"] <- beta["(Intercept)"]
  }
  out[!is.finite(out)] <- 0
  out
}

make_design_matrix <- function(data, covariates, intercept_name = "Intercept") {
  X <- model.matrix(as.formula(paste("~", paste(covariates, collapse = " + "))), data = data)
  colnames(X)[colnames(X) == "(Intercept)"] <- intercept_name
  X
}

linear_lp_from_beta <- function(beta, newdata, covariates = paste0("X", 1:7)) {
  X <- make_design_matrix(newdata, covariates)
  beta_aligned <- align_beta_to_X(beta, colnames(X))
  as.vector(X %*% beta_aligned)
}

make_scaler <- function(train_data, cont_covariates = c("X3", "X4", "X5", "X6", "X7")) {
  scaler <- lapply(cont_covariates, function(v) {
    s <- stats::sd(train_data[[v]], na.rm = TRUE)
    if (is.na(s) || s < 1e-8) s <- 1
    list(mean = mean(train_data[[v]], na.rm = TRUE), sd = s)
  })
  names(scaler) <- cont_covariates
  scaler
}

apply_scaler <- function(data, scaler, suffix = "_sc") {
  out <- data
  for (v in names(scaler)) {
    out[[paste0(v, suffix)]] <- (out[[v]] - scaler[[v]]$mean) / scaler[[v]]$sd
  }
  out
}

clip_probs <- function(p, eps = 1e-3) pmin(pmax(p, eps), 1 - eps)

# ==============================================================================
# 2. IPCW and fixed-time binary outcome
# ==============================================================================

estimate_weights <- function(data, t, Y_col = "Y", delta_col = "delta", eps_G = 1e-3, w_cap = Inf) {
  Y_vals <- data[[Y_col]]
  delta_vals <- data[[delta_col]]
  km_fit <- survival::survfit(survival::Surv(Y_vals, 1 - delta_vals) ~ 1, data = data)
  time_point <- pmin(Y_vals, t)
  Ghat <- vapply(time_point, function(x) {
    s <- summary(km_fit, times = x)$surv
    if (length(s) == 0 || is.na(s)) eps_G else max(s, eps_G)
  }, numeric(1))
  indicator <- ifelse((Y_vals > t) | ((Y_vals <= t) & (delta_vals == 1)), 1, 0)
  weights <- indicator / Ghat
  if (is.finite(w_cap)) weights <- pmin(weights, w_cap)
  weights
}

make_Dt_weights <- function(data, t, Y_col = "Y", delta_col = "delta", eps_G = 1e-3, w_cap = 20) {
  list(
    D_t = as.integer(data[[Y_col]] <= t & data[[delta_col]] == 1),
    weights = estimate_weights(data, t, Y_col, delta_col, eps_G, w_cap)
  )
}

# ==============================================================================
# 3. Source models and covariate-shift weighting
# ==============================================================================

estimate_betaS <- function(data, t, Y_col = "Y", delta_col = "delta", covariates = NULL, eps_G = 1e-3, w_cap = Inf) {
  aux <- make_Dt_weights(data, t, Y_col, delta_col, eps_G, w_cap)
  df <- data
  df$D_t <- aux$D_t
  df$w_ipcw <- aux$weights
  if (is.null(covariates)) covariates <- setdiff(names(df), c(Y_col, delta_col, "D_t", "w_ipcw"))
  fit <- glm(
    as.formula(paste("D_t ~", paste(covariates, collapse = " + "))),
    data = df, family = quasibinomial(), weights = w_ipcw,
    control = glm.control(maxit = 100)
  )
  beta <- stats::coef(fit)
  names(beta)[names(beta) == "(Intercept)"] <- "Intercept"
  beta[!is.finite(beta)] <- 0
  beta
}

estimate_wcs <- function(source_data, target_data, covariates, lambda = 1,
                         method = c("tilt", "logistic"), stabilize = TRUE,
                         eps_p = 1e-3, clip_eta = 30, nfolds_dom = 5) {
  method <- match.arg(method)
  Z_src <- as.matrix(cbind(Intercept = 1, source_data[, covariates, drop = FALSE]))
  Z_tgt <- as.matrix(cbind(Intercept = 1, target_data[, covariates, drop = FALSE]))
  nS <- nrow(Z_src); nT <- nrow(Z_tgt)

  if (method == "tilt") {
    mu_tgt <- colMeans(Z_tgt)
    clip <- function(x, lo = -clip_eta, hi = clip_eta) pmin(pmax(x, lo), hi)
    objective <- function(theta) {
      theta <- as.numeric(theta)
      eta <- clip(as.vector(Z_src %*% theta))
      -(sum(mu_tgt * theta) - mean(exp(eta)) - lambda * sum(theta^2))
    }
    fit <- nlm(f = objective, p = rep(0, ncol(Z_src)), iterlim = 500)
    theta_hat <- fit$estimate
    w_src <- exp(clip(as.vector(Z_src %*% theta_hat)))
    if (stabilize) w_src <- w_src / mean(w_src)
    return(list(theta = theta_hat, w_source = w_src, w_target = rep(1, nT),
                converged = fit$code == 1, nlm_code = fit$code, method = "tilt"))
  }

  X_all <- rbind(source_data[, covariates, drop = FALSE], target_data[, covariates, drop = FALSE])
  X_dm <- model.matrix(~ ., data = X_all)[, -1, drop = FALSE]
  R <- c(rep(0, nS), rep(1, nT))
  nfolds_use <- choose_nfolds_binomial(R, max_folds = nfolds_dom)
  cv_fit <- glmnet::cv.glmnet(X_dm, R, family = "binomial", alpha = 0,
                              nfolds = nfolds_use, type.measure = "deviance", standardize = TRUE)
  fit <- glmnet::glmnet(X_dm, R, family = "binomial", alpha = 0,
                        lambda = cv_fit$lambda.min, standardize = TRUE)
  p_hat <- clip_probs(as.vector(stats::predict(fit, X_dm, type = "response")), eps_p)
  pi_hat <- mean(R == 1)
  w_src <- ((p_hat / (1 - p_hat)) * ((1 - pi_hat) / pi_hat))[1:nS]
  if (stabilize) w_src <- w_src / mean(w_src)
  list(theta = as.numeric(stats::coef(fit)), w_source = w_src, w_target = rep(1, nT),
       converged = TRUE, method = "logistic", lambda_dom = cv_fit$lambda.min,
       nfolds_dom = nfolds_use)
}

estimate_wcs_tilt_safe <- function(source_data, target_data, covariates, lambda = 1, clip_eta = 30) {
  estimate_wcs(source_data, target_data, covariates, lambda = lambda, method = "tilt", clip_eta = clip_eta)
}

.fit_weighted_source_logistic <- function(source_data, t, covariates, w_tot, Y_col = "Y", delta_col = "delta") {
  D_t <- as.integer(source_data[[Y_col]] <= t & source_data[[delta_col]] == 1)
  X <- make_design_matrix(source_data, covariates)
  dfX <- as.data.frame(X)
  dfX$D_t <- D_t
  dfX$w_tot <- w_tot
  fit <- glm(
    as.formula(paste("D_t ~ 0 +", paste(colnames(X), collapse = " + "))),
    data = dfX, family = quasibinomial(), weights = w_tot,
    control = glm.control(maxit = 100)
  )
  b <- stats::coef(fit)
  out_names <- c("Intercept", covariates)
  beta_out <- rep(0, length(out_names)); names(beta_out) <- out_names
  beta_out[intersect(names(b), out_names)] <- b[intersect(names(b), out_names)]
  beta_out[!is.finite(beta_out)] <- 0
  beta_out
}

estimate_betaS_tcs_call <- function(source_data, target_data, t, Y_col = "Y", delta_col = "delta",
                                    covariates, lambda_w = 1, w_cap = 50, clip_eta = 30,
                                    eps_G = 1e-3) {
  w_obj <- estimate_wcs_tilt_safe(source_data, target_data, covariates, lambda_w, clip_eta)
  w_cs <- as.numeric(w_obj$w_source)
  w_ipcw <- estimate_weights(source_data, t, Y_col, delta_col, eps_G, w_cap)
  w_tot <- w_cs * w_ipcw
  if (is.finite(w_cap)) w_tot <- pmin(w_tot, w_cap)
  beta_out <- .fit_weighted_source_logistic(source_data, t, covariates, w_tot, Y_col, delta_col)
  list(beta = beta_out, theta_tilt = w_obj$theta, w_cs = w_cs, w_ipcw = w_ipcw,
       w_tot = w_tot, converged = w_obj$converged, nlm_code = w_obj$nlm_code)
}

estimate_betaS_cs <- function(source_data, target_data, t, Y_col = "Y", delta_col = "delta",
                              covariates, eps_G = 1e-3, w_cap = 50, max_folds = 5) {
  w_obj <- estimate_wcs(source_data, target_data, covariates, method = "logistic", nfolds_dom = max_folds)
  w_cs <- as.numeric(w_obj$w_source)
  w_ipcw <- estimate_weights(source_data, t, Y_col, delta_col, eps_G, w_cap)
  w_tot <- w_cs * w_ipcw
  if (is.finite(w_cap)) w_tot <- pmin(w_tot, w_cap)
  beta_out <- .fit_weighted_source_logistic(source_data, t, covariates, w_tot, Y_col, delta_col)
  list(beta = beta_out, w_cs = w_cs, w_ipcw = w_ipcw, w_tot = w_tot,
       lambda_dom = w_obj$lambda_dom, nfolds_dom = w_obj$nfolds_dom)
}

# ==============================================================================
# 4. Linear target adaptation
# ==============================================================================

estimate_eta_beta <- function(betaS, target_data, t, Y_col = "Y", delta_col = "delta",
                              covariates, eps_G = 1e-3, w_cap = Inf, max_folds = 5,
                              alpha_eta = 1, lambda_choice = c("lambda.1se", "lambda.min"),
                              standardize_eta = TRUE) {
  lambda_choice <- match.arg(lambda_choice)
  aux <- make_Dt_weights(target_data, t, Y_col, delta_col, eps_G, w_cap)
  X_target <- make_design_matrix(target_data, covariates)
  betaS_aligned <- align_beta_to_X(betaS, colnames(X_target))
  offset_target <- as.vector(X_target %*% betaS_aligned)
  nfolds_use <- choose_nfolds_binomial(aux$D_t, max_folds)
  cv_fit <- glmnet::cv.glmnet(
    x = X_target, y = aux$D_t, family = "binomial", weights = aux$weights,
    offset = offset_target, alpha = alpha_eta, nfolds = nfolds_use,
    type.measure = "deviance", intercept = FALSE, standardize = standardize_eta
  )
  lambda_use <- if (lambda_choice == "lambda.1se") cv_fit$lambda.1se else cv_fit$lambda.min
  fit_final <- glmnet::glmnet(
    x = X_target, y = aux$D_t, family = "binomial", weights = aux$weights,
    offset = offset_target, alpha = alpha_eta, lambda = lambda_use,
    intercept = FALSE, standardize = standardize_eta
  )
  coef_mat <- as.matrix(stats::coef(fit_final)); rn <- rownames(coef_mat)
  eta_est <- rep(0, ncol(X_target)); names(eta_est) <- colnames(X_target)
  for (nm in names(eta_est)) {
    if (nm %in% rn) eta_est[nm] <- coef_mat[nm, 1]
    else if (nm == "Intercept" && "(Intercept)" %in% rn) eta_est[nm] <- coef_mat["(Intercept)", 1]
  }
  list(eta = eta_est, beta = betaS_aligned + eta_est, lambda = lambda_use,
       lambda_choice = lambda_choice, alpha_eta = alpha_eta, nfolds = nfolds_use,
       weights = aux$weights, offset = offset_target,
       n_event = sum(aux$D_t == 1), n_nonevent = sum(aux$D_t == 0))
}

estimate_eta_beta_df <- function(betaS, target_data, t, Y_col = "Y", delta_col = "delta",
                                 common_cov, target_only_cov = character(0), eps_G = 1e-3,
                                 w_cap = Inf, max_folds = 5, alpha_eta = 1,
                                 lambda_choice = c("lambda.1se", "lambda.min"),
                                 standardize_eta = TRUE) {
  lambda_choice <- match.arg(lambda_choice)
  X_common <- as.matrix(target_data[, common_cov, drop = FALSE])
  X_only <- as.matrix(target_data[, target_only_cov, drop = FALSE])
  X_target <- cbind(Intercept = 1, X_common, X_only)
  aux <- make_Dt_weights(target_data, t, Y_col, delta_col, eps_G, w_cap)
  betaS_common <- align_beta_to_X(betaS, c("Intercept", common_cov))
  betaS_padded <- c(betaS_common, rep(0, ncol(X_only))); names(betaS_padded) <- colnames(X_target)
  k0 <- 1 + length(common_cov)
  offset <- as.vector(X_target[, 1:k0, drop = FALSE] %*% betaS_padded[1:k0])
  penalty_factor <- c(0, rep(1, length(common_cov)), rep(0, ncol(X_only)))
  nfolds_use <- choose_nfolds_binomial(aux$D_t, max_folds)
  cv_fit <- glmnet::cv.glmnet(
    x = X_target, y = aux$D_t, family = "binomial", weights = aux$weights,
    offset = offset, penalty.factor = penalty_factor, alpha = alpha_eta,
    nfolds = nfolds_use, type.measure = "deviance", intercept = FALSE,
    standardize = standardize_eta
  )
  lambda_use <- if (lambda_choice == "lambda.1se") cv_fit$lambda.1se else cv_fit$lambda.min
  fit <- glmnet::glmnet(
    x = X_target, y = aux$D_t, family = "binomial", weights = aux$weights,
    offset = offset, penalty.factor = penalty_factor, alpha = alpha_eta,
    lambda = lambda_use, intercept = FALSE, standardize = standardize_eta
  )
  coef_mat <- as.matrix(stats::coef(fit)); rn <- rownames(coef_mat)
  eta_vec <- rep(0, ncol(X_target)); names(eta_vec) <- colnames(X_target)
  for (nm in names(eta_vec)) {
    if (nm %in% rn) eta_vec[nm] <- coef_mat[nm, 1]
    else if (nm == "Intercept" && "(Intercept)" %in% rn) eta_vec[nm] <- coef_mat["(Intercept)", 1]
  }
  list(beta = betaS_padded + eta_vec, eta = eta_vec, lambda = lambda_use,
       lambda_choice = lambda_choice, alpha_eta = alpha_eta, nfolds = nfolds_use,
       weights = aux$weights, offset = offset,
       n_event = sum(aux$D_t == 1), n_nonevent = sum(aux$D_t == 0))
}

# ==============================================================================
# 5. Transfer GAM family, S3-style
# ============================================================================== 

make_source_gam_formula <- function(use_interaction = FALSE, include_binary = TRUE,
                                    smooth_bs = c("ts", "tp")) {
  smooth_bs <- match.arg(smooth_bs)
  rhs <- character(0)
  if (include_binary) rhs <- c(rhs, "X1", "X2")
  rhs <- c(rhs, sprintf("s(%s_sc, k = 5, bs = '%s')", paste0("X", 3:7), smooth_bs))
  if (use_interaction) rhs <- c(rhs, sprintf("s(X3_sc, by = X1, k = 5, bs = '%s')", smooth_bs))
  as.formula(paste("D_t ~", paste(rhs, collapse = " + ")))
}

make_target_gam_formula <- function(calibrated_source = FALSE, use_interaction = FALSE,
                                    include_binary = TRUE, source_term = NULL,
                                    smooth_bs = c("ts", "tp")) {
  smooth_bs <- match.arg(smooth_bs)
  rhs <- if (!is.null(source_term)) {
    source_term
  } else if (calibrated_source) {
    "source_lp"
  } else {
    "offset(source_lp)"
  }
  if (include_binary) rhs <- c(rhs, "X1", "X2")
  rhs <- c(rhs, sprintf("s(%s_sc, k = 5, bs = '%s')", paste0("X", 3:7), smooth_bs))
  if (use_interaction) rhs <- c(rhs, sprintf("s(X3_sc, by = X1, k = 5, bs = '%s')", smooth_bs))
  as.formula(paste("D_t ~", paste(rhs, collapse = " + ")))
}

make_source_recal_formula <- function(recalibration = c("intercept", "slope")) {
  recalibration <- match.arg(recalibration)
  if (recalibration == "intercept") return(stats::as.formula("D_t ~ 1 + offset(source_lp)"))
  stats::as.formula("D_t ~ source_lp")
}

extract_source_alpha <- function(object) {
  out <- NA_real_
  valid_alpha <- function(x) {
    length(x) > 0 && is.finite(as.numeric(x)[1])
  }

  if (!is.null(object$source_alpha) && valid_alpha(object$source_alpha)) {
    return(as.numeric(object$source_alpha)[1])
  }

  if (!is.null(object$target_obj) &&
      !is.null(object$target_obj$source_alpha) &&
      valid_alpha(object$target_obj$source_alpha)) {
    return(as.numeric(object$target_obj$source_alpha)[1])
  }

  if (!is.null(object$target_obj) && !is.null(object$target_obj$fit)) {
    cc <- stats::coef(object$target_obj$fit)
    if ("source_lp" %in% names(cc)) out <- as.numeric(cc["source_lp"])
    if ("source_lp_delta" %in% names(cc)) out <- 1 + as.numeric(cc["source_lp_delta"])
  }

  if (!is.null(object$fit)) {
    cc <- stats::coef(object$fit)
    if ("source_lp" %in% names(cc)) out <- as.numeric(cc["source_lp"])
    if ("source_lp_delta" %in% names(cc)) out <- 1 + as.numeric(cc["source_lp_delta"])
  }

  out
}

source_weight_summary <- function(w) {
  w <- as.numeric(w)
  if (length(w) == 0 || all(is.na(w))) {
    return(list(
      min = NA_real_, q25 = NA_real_, median = NA_real_, q75 = NA_real_,
      q95 = NA_real_, max = NA_real_, ess = NA_real_
    ))
  }
  list(
    min = unname(stats::quantile(w, 0.00, na.rm = TRUE)),
    q25 = unname(stats::quantile(w, 0.25, na.rm = TRUE)),
    median = unname(stats::quantile(w, 0.50, na.rm = TRUE)),
    q75 = unname(stats::quantile(w, 0.75, na.rm = TRUE)),
    q95 = unname(stats::quantile(w, 0.95, na.rm = TRUE)),
    max = unname(stats::quantile(w, 1.00, na.rm = TRUE)),
    ess = (sum(w, na.rm = TRUE)^2) / sum(w^2, na.rm = TRUE)
  )
}

fit_source_linear_score <- function(source_data, target_train, t, covariates = paste0("X", 1:7),
                                    use_tcs_weight = TRUE, lambda_w = 1,
                                    eps_G = 1e-3, w_cap = 20) {
  if (use_tcs_weight) {
    fit <- estimate_betaS_tcs_call(source_data, target_train, t, covariates = covariates,
                                   lambda_w = lambda_w, w_cap = w_cap, eps_G = eps_G)
    betaS <- fit$beta
  } else {
    betaS <- estimate_betaS(source_data, t, covariates = covariates, eps_G = eps_G, w_cap = w_cap)
    fit <- list(beta = betaS)
  }
  list(beta = betaS, fit = fit, source_type = "linear", use_tcs_weight = use_tcs_weight)
}

fit_source_gam_score <- function(source_data, target_train = NULL, t,
                                 covariates = paste0("X", 1:7),
                                 cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                 use_interaction = FALSE, include_binary = TRUE,
                                 use_tcs_weight = TRUE, lambda_w = 1, clip_eta = 30,
                                 eps_G = 1e-3, w_cap = 20,
                                 gamma = 1.4, select = FALSE,
                                 smooth_bs = c("ts", "tp")) {
  .require_pkg("mgcv")
  smooth_bs <- match.arg(smooth_bs)

  aux <- make_Dt_weights(source_data, t, eps_G = eps_G, w_cap = w_cap)
  w_ipcw <- aux$weights
  w_cs <- rep(1, nrow(source_data))
  w_obj <- NULL

  if (isTRUE(use_tcs_weight)) {
    if (is.null(target_train)) {
      stop("target_train is required when use_tcs_weight = TRUE for source GAM fitting.")
    }
    w_obj <- estimate_wcs_tilt_safe(
      source_data = source_data,
      target_data = target_train,
      covariates = covariates,
      lambda = lambda_w,
      clip_eta = clip_eta
    )
    w_cs <- as.numeric(w_obj$w_source)
  }

  w_tot <- w_ipcw * w_cs
  if (is.finite(w_cap)) w_tot <- pmin(w_tot, w_cap)

  df <- source_data
  df$D_t <- aux$D_t
  df$w_ipcw <- w_ipcw
  df$w_cs <- w_cs
  df$w_tot <- w_tot

  scaler <- make_scaler(df, cont_covariates)
  df <- apply_scaler(df, scaler)

  fit <- mgcv::gam(
    formula = make_source_gam_formula(use_interaction, include_binary, smooth_bs = smooth_bs),
    data = df,
    family = quasibinomial(),
    weights = w_tot,
    method = "REML",
    gamma = gamma,
    select = select
  )

  list(
    fit = fit,
    scaler = scaler,
    source_type = "gam",
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    smooth_bs = smooth_bs,
    use_tcs_weight = use_tcs_weight,
    w_obj = w_obj,
    w_ipcw = w_ipcw,
    w_cs = w_cs,
    w_tot = w_tot,
    weight_summary = source_weight_summary(w_tot),
    n_event = sum(aux$D_t == 1),
    n_nonevent = sum(aux$D_t == 0)
  )
}

predict_source_gam_lp <- function(source_gam_obj, newdata) {
  newdata_sc <- apply_scaler(newdata, source_gam_obj$scaler)
  as.vector(stats::predict(source_gam_obj$fit, newdata = newdata_sc, type = "link"))
}

fit_target_gam_correction <- function(target_train, t, source_lp_train, calibrated_source = FALSE,
                                      cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                      use_interaction = FALSE, include_binary = TRUE,
                                      eps_G = 1e-3, w_cap = 20,
                                      gamma = 1.4, select = FALSE,
                                      smooth_bs = c("ts", "tp")) {
  .require_pkg("mgcv")
  smooth_bs <- match.arg(smooth_bs)
  aux <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  df <- target_train
  df$D_t <- aux$D_t
  df$w_ipcw <- aux$weights
  df$source_lp <- as.numeric(source_lp_train)

  scaler <- make_scaler(df, cont_covariates)
  df <- apply_scaler(df, scaler)

  fit <- mgcv::gam(
    formula = make_target_gam_formula(calibrated_source, use_interaction, include_binary,
                                      smooth_bs = smooth_bs),
    data = df,
    family = quasibinomial(),
    weights = w_ipcw,
    method = "REML",
    gamma = gamma,
    select = select
  )

  obj <- list(
    fit = fit,
    scaler = scaler,
    calibrated_source = calibrated_source,
    source_alpha = NA_real_,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    smooth_bs = smooth_bs,
    n_event = sum(aux$D_t == 1),
    n_nonevent = sum(aux$D_t == 0)
  )
  if (isTRUE(calibrated_source)) obj$source_alpha <- extract_source_alpha(obj)
  obj
}

fit_target_gam_penalized_slope_correction <- function(target_train, t, source_lp_train,
                                                       cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                                       use_interaction = FALSE, include_binary = TRUE,
                                                       eps_G = 1e-3, w_cap = 20,
                                                       gamma = 1.4, select = FALSE,
                                                       alpha_penalty = getOption("translogistic.alpha_penalty", 10)) {
  .require_pkg("mgcv")
  aux <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  df <- target_train
  df$D_t <- aux$D_t
  df$w_ipcw <- aux$weights
  df$source_lp <- as.numeric(source_lp_train)
  df$source_lp_delta <- as.numeric(source_lp_train)

  scaler <- make_scaler(df, cont_covariates)
  df <- apply_scaler(df, scaler)

  alpha_penalty <- as.numeric(alpha_penalty)[1]
  if (!is.finite(alpha_penalty)) alpha_penalty <- -1
  sp_alpha <- if (alpha_penalty > 0) alpha_penalty else -1

  fit <- mgcv::gam(
    formula = make_target_gam_formula(
      calibrated_source = FALSE,
      use_interaction = use_interaction,
      include_binary = include_binary,
      source_term = c("offset(source_lp)", "source_lp_delta")
    ),
    data = df,
    family = quasibinomial(),
    weights = w_ipcw,
    method = "REML",
    gamma = gamma,
    select = select,
    paraPen = list(source_lp_delta = list(matrix(1, 1, 1), sp = sp_alpha))
  )

  alpha_penalty_estimated <- NA_real_
  if (!is.null(fit$sp) && length(fit$sp) > 0) {
    sp_names <- names(fit$sp)
    idx <- grep("source_lp_delta", sp_names, fixed = TRUE)
    if (length(idx) > 0) {
      alpha_penalty_estimated <- as.numeric(fit$sp[idx[1]])
    } else if (sp_alpha < 0) {
      alpha_penalty_estimated <- as.numeric(fit$sp[length(fit$sp)])
    }
  }

  obj <- list(
    fit = fit,
    scaler = scaler,
    calibrated_source = TRUE,
    source_calibration = "penalized_slope",
    source_alpha = extract_source_alpha(list(fit = fit)),
    alpha_penalty = alpha_penalty,
    alpha_penalty_estimated = alpha_penalty_estimated,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    n_event = sum(aux$D_t == 1),
    n_nonevent = sum(aux$D_t == 0)
  )
  obj
}

predict_target_gam_lp <- function(target_gam_obj, newdata, source_lp_new) {
  pred_data <- newdata
  pred_data$source_lp <- as.numeric(source_lp_new)
  if (identical(target_gam_obj$source_calibration, "penalized_slope")) {
    pred_data$source_lp_delta <- as.numeric(source_lp_new)
  }
  pred_data <- apply_scaler(pred_data, target_gam_obj$scaler)
  as.vector(stats::predict(target_gam_obj$fit, newdata = pred_data, type = "link"))
}

estimate_linear_correction_from_source_lp <- function(source_lp_train, target_data, t,
                                                      covariates = paste0("X", 1:7),
                                                      calibrated_source = FALSE,
                                                      Y_col = "Y", delta_col = "delta",
                                                      eps_G = 1e-3, w_cap = 20,
                                                      max_folds = 5, alpha_eta = 1,
                                                      lambda_choice = c("lambda.1se", "lambda.min"),
                                                      standardize_eta = TRUE) {
  lambda_choice <- match.arg(lambda_choice)
  aux <- make_Dt_weights(target_data, t, Y_col, delta_col, eps_G, w_cap)
  X_target <- make_design_matrix(target_data, covariates)
  nfolds_use <- choose_nfolds_binomial(aux$D_t, max_folds)

  if (isTRUE(calibrated_source)) {
    x <- cbind(source_lp = as.numeric(source_lp_train), X_target)
    offset <- rep(0, nrow(X_target))
    penalty_factor <- c(0, 0, rep(1, ncol(X_target) - 1))
  } else {
    x <- X_target
    offset <- as.numeric(source_lp_train)
    penalty_factor <- c(0, rep(1, ncol(X_target) - 1))
  }

  cv_fit <- glmnet::cv.glmnet(
    x = x,
    y = aux$D_t,
    family = "binomial",
    weights = aux$weights,
    offset = offset,
    penalty.factor = penalty_factor,
    alpha = alpha_eta,
    nfolds = nfolds_use,
    type.measure = "deviance",
    intercept = FALSE,
    standardize = standardize_eta
  )

  lambda_use <- if (lambda_choice == "lambda.1se") cv_fit$lambda.1se else cv_fit$lambda.min

  fit <- glmnet::glmnet(
    x = x,
    y = aux$D_t,
    family = "binomial",
    weights = aux$weights,
    offset = offset,
    penalty.factor = penalty_factor,
    alpha = alpha_eta,
    lambda = lambda_use,
    intercept = FALSE,
    standardize = standardize_eta
  )

  coef_mat <- as.matrix(stats::coef(fit))
  rn <- rownames(coef_mat)

  source_alpha <- if (isTRUE(calibrated_source) && "source_lp" %in% rn) {
    as.numeric(coef_mat["source_lp", 1])
  } else if (isTRUE(calibrated_source)) {
    0
  } else {
    1
  }

  eta_names <- colnames(X_target)
  eta <- rep(0, length(eta_names)); names(eta) <- eta_names
  for (nm in eta_names) {
    if (nm %in% rn) eta[nm] <- coef_mat[nm, 1]
    else if (nm == "Intercept" && "(Intercept)" %in% rn) eta[nm] <- coef_mat["(Intercept)", 1]
  }

  obj <- list(
    fit = fit,
    cv_fit = cv_fit,
    x_names = colnames(x),
    covariates = covariates,
    calibrated_source = calibrated_source,
    source_alpha = source_alpha,
    eta = eta,
    lambda = lambda_use,
    lambda_choice = lambda_choice,
    alpha_eta = alpha_eta,
    nfolds = nfolds_use,
    weights = aux$weights,
    offset_train = offset,
    n_event = sum(aux$D_t == 1),
    n_nonevent = sum(aux$D_t == 0)
  )
  class(obj) <- "source_lp_linear_correction"
  obj
}

predict.source_lp_linear_correction <- function(object, newdata, source_lp_new,
                                                type = c("link", "response"), ...) {
  type <- match.arg(type)
  X_new <- make_design_matrix(newdata, object$covariates)

  if (isTRUE(object$calibrated_source)) {
    x_new <- cbind(source_lp = as.numeric(source_lp_new), X_new)
    offset_new <- rep(0, nrow(X_new))
  } else {
    x_new <- X_new
    offset_new <- as.numeric(source_lp_new)
  }

  x_new <- x_new[, object$x_names, drop = FALSE]
  lp <- as.vector(stats::predict(object$fit, newx = x_new, s = object$lambda,
                                 type = "link", newoffset = offset_new))
  if (type == "link") return(lp)
  stats::plogis(lp)
}

fit_source_gam_recalibration <- function(source_data, target_train, t,
                                         recalibration = c("intercept", "slope"),
                                         covariates = paste0("X", 1:7),
                                         cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                         use_interaction = FALSE, include_binary = TRUE,
                                         use_tcs_weight = TRUE, lambda_w = 1,
                                         eps_G = 1e-3, w_cap = 20,
                                         gamma_source = 1.4, select_source = FALSE) {
  recalibration <- match.arg(recalibration)

  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )

  source_lp_train <- predict_source_gam_lp(source_obj, target_train)
  aux <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)

  df <- target_train
  df$D_t <- aux$D_t
  df$w_ipcw <- aux$weights
  df$source_lp <- as.numeric(source_lp_train)

  fit <- stats::glm(
    formula = make_source_recal_formula(recalibration),
    data = df,
    family = quasibinomial(),
    weights = w_ipcw,
    control = stats::glm.control(maxit = 100)
  )

  obj <- list(
    method = paste0("source_gam_recal_", recalibration),
    source_obj = source_obj,
    fit = fit,
    recalibration = recalibration,
    source_alpha = if (recalibration == "slope") as.numeric(stats::coef(fit)["source_lp"]) else 1,
    intercept = as.numeric(stats::coef(fit)["(Intercept)"]),
    covariates = covariates,
    cont_covariates = cont_covariates,
    n_event = sum(aux$D_t == 1),
    n_nonevent = sum(aux$D_t == 0)
  )
  class(obj) <- "source_gam_recalibration_model"
  obj
}

predict.source_gam_recalibration_model <- function(object, newdata,
                                                   type = c("link", "response"), ...) {
  type <- match.arg(type)
  source_lp <- predict_source_gam_lp(object$source_obj, newdata)
  pred_data <- data.frame(source_lp = as.numeric(source_lp))
  lp <- as.vector(stats::predict(object$fit, newdata = pred_data, type = "link"))
  if (type == "link") return(lp)
  stats::plogis(lp)
}

predict_source_gam_recalibration_risk_t <- function(object, newdata) {
  lp <- predict(object, newdata = newdata, type = "link")
  list(lp = lp, p_event_t = stats::plogis(lp))
}

fit_trans_linear_linear <- function(source_data, target_train, t, covariates = paste0("X", 1:7),
                                    lambda_w = 1, eps_G = 1e-3, w_cap = 20,
                                    max_folds = 5, alpha_eta = 1,
                                    lambda_choice = "lambda.1se") {
  source_obj <- fit_source_linear_score(source_data, target_train, t, covariates,
                                        use_tcs_weight = TRUE, lambda_w = lambda_w,
                                        eps_G = eps_G, w_cap = w_cap)
  fit_eta <- estimate_eta_beta(source_obj$beta, target_train, t, covariates = covariates,
                               eps_G = eps_G, w_cap = w_cap, max_folds = max_folds,
                               alpha_eta = alpha_eta, lambda_choice = lambda_choice,
                               standardize_eta = TRUE)
  obj <- list(method = "trans_linear_linear", source_obj = source_obj, target_obj = fit_eta,
              beta = fit_eta$beta, covariates = covariates, source_alpha = 1)
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_linear_gam <- function(source_data, target_train, t, covariates = paste0("X", 1:7),
                                 cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                 use_interaction = FALSE, include_binary = TRUE,
                                 calibrated_source = FALSE,
                                 use_tcs_weight = TRUE, lambda_w = 1,
                                 eps_G = 1e-3, w_cap = 20,
                                 gamma = 1.4, select = FALSE,
                                 method_name = NULL) {
  source_obj <- fit_source_linear_score(source_data, target_train, t, covariates,
                                        use_tcs_weight = use_tcs_weight, lambda_w = lambda_w,
                                        eps_G = eps_G, w_cap = w_cap)
  source_lp_train <- linear_lp_from_beta(source_obj$beta, target_train, covariates)
  target_obj <- fit_target_gam_correction(
    target_train = target_train,
    t = t,
    source_lp_train = source_lp_train,
    calibrated_source = calibrated_source,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma,
    select = select
  )
  if (is.null(method_name)) {
    method_name <- if (calibrated_source) "trans_linear_gam_calibrated" else "trans_linear_gam"
  }
  obj <- list(method = method_name, source_obj = source_obj, target_obj = target_obj,
              covariates = covariates, cont_covariates = cont_covariates,
              source_alpha = extract_source_alpha(target_obj))
  if (!calibrated_source) obj$source_alpha <- 1
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_gam_linear <- function(source_data, target_train, t, covariates = paste0("X", 1:7),
                                 cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                 use_interaction = FALSE, include_binary = TRUE,
                                 calibrated_source = FALSE,
                                 use_tcs_weight = TRUE, lambda_w = 1,
                                 eps_G = 1e-3, w_cap = 20,
                                 max_folds = 5, alpha_eta = 1,
                                 lambda_choice = "lambda.1se",
                                 gamma_source = 1.4, select_source = FALSE,
                                 method_name = NULL) {
  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)
  target_obj <- estimate_linear_correction_from_source_lp(
    source_lp_train = source_lp_train,
    target_data = target_train,
    t = t,
    covariates = covariates,
    calibrated_source = calibrated_source,
    eps_G = eps_G,
    w_cap = w_cap,
    max_folds = max_folds,
    alpha_eta = alpha_eta,
    lambda_choice = lambda_choice,
    standardize_eta = TRUE
  )
  if (is.null(method_name)) {
    method_name <- if (calibrated_source) "trans_gam_linear_calibrated" else "trans_gam_linear"
  }
  obj <- list(method = method_name, source_obj = source_obj, target_obj = target_obj,
              covariates = covariates, cont_covariates = cont_covariates,
              source_alpha = target_obj$source_alpha)
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_gam_gam <- function(source_data, target_train, t, covariates = paste0("X", 1:7),
                              cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                              use_interaction = FALSE, include_binary = TRUE,
                              calibrated_source = FALSE,
                              use_tcs_weight = TRUE, lambda_w = 1,
                              eps_G = 1e-3, w_cap = 20,
                              gamma_source = 1.4, gamma_target = 1.4,
                              select_source = FALSE, select_target = FALSE,
                              source_smooth_bs = c("ts", "tp"),
                              target_smooth_bs = c("ts", "tp"),
                              method_name = NULL) {
  source_smooth_bs <- match.arg(source_smooth_bs)
  target_smooth_bs <- match.arg(target_smooth_bs)
  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source,
    smooth_bs = source_smooth_bs
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)
  target_obj <- fit_target_gam_correction(
    target_train = target_train,
    t = t,
    source_lp_train = source_lp_train,
    calibrated_source = calibrated_source,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_target,
    select = select_target,
    smooth_bs = target_smooth_bs
  )
  if (is.null(method_name)) {
    method_name <- if (calibrated_source) "trans_gam_gam_calibrated" else "trans_gam_gam"
    if (!isTRUE(use_tcs_weight)) method_name <- paste0(method_name, "_no_tilt")
  }
  obj <- list(method = method_name, source_obj = source_obj, target_obj = target_obj,
              covariates = covariates, cont_covariates = cont_covariates,
              source_smooth_bs = source_smooth_bs,
              target_smooth_bs = target_smooth_bs,
              source_alpha = extract_source_alpha(target_obj))
  if (!calibrated_source) obj$source_alpha <- 1
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_gam_gam_calibrated <- function(source_data, target_train, t, covariates = paste0("X", 1:7),
                                         cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                         use_interaction = FALSE, include_binary = TRUE,
                                         use_tcs_weight = TRUE, lambda_w = 1,
                                         eps_G = 1e-3, w_cap = 20,
                                         gamma_source = 1.4, gamma_target = 1.4,
                                         select_source = FALSE, select_target = FALSE) {
  fit_trans_gam_gam(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    calibrated_source = TRUE,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma_source = gamma_source,
    gamma_target = gamma_target,
    select_source = select_source,
    select_target = select_target,
    method_name = if (use_tcs_weight) "trans_gam_gam_calibrated" else "trans_gam_gam_calibrated_no_tilt"
  )
}

cv_select_partial_gam_calibration <- function(source_lp_train,
                                              target_train,
                                              t,
                                              fold_id,
                                              partial_grid = c(0, 0.25, 0.5, 0.75, 1),
                                              cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                              use_interaction = FALSE,
                                              include_binary = TRUE,
                                              eps_G = 1e-3,
                                              w_cap = 20,
                                              gamma_target = 1.4,
                                              select_target = FALSE,
                                              auto_metric = c("logloss", "brier", "auc")) {
  auto_metric <- match.arg(auto_metric)
  higher_is_better <- identical(auto_metric, "auc")
  failed_score <- if (higher_is_better) -Inf else Inf
  partial_grid <- sort(unique(as.numeric(partial_grid)))
  partial_grid <- partial_grid[is.finite(partial_grid) & partial_grid >= 0 & partial_grid <= 1]
  if (length(partial_grid) == 0) partial_grid <- c(0, 0.25, 0.5, 0.75, 1)

  aux_full <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  folds <- sort(unique(fold_id[is.finite(fold_id)]))
  score_mat <- matrix(NA_real_, nrow = length(folds), ncol = length(partial_grid))
  colnames(score_mat) <- paste0("w_", format(partial_grid, trim = TRUE, scientific = FALSE))

  for (i in seq_along(folds)) {
    fold <- folds[i]
    valid_idx <- which(fold_id == fold)
    train_idx <- which(fold_id != fold)

    score_mat[i, ] <- tryCatch({
      target_no_cal <- fit_target_gam_correction(
        target_train = target_train[train_idx, , drop = FALSE],
        t = t,
        source_lp_train = source_lp_train[train_idx],
        calibrated_source = FALSE,
        cont_covariates = cont_covariates,
        use_interaction = use_interaction,
        include_binary = include_binary,
        eps_G = eps_G,
        w_cap = w_cap,
        gamma = gamma_target,
        select = select_target
      )
      target_cal <- fit_target_gam_correction(
        target_train = target_train[train_idx, , drop = FALSE],
        t = t,
        source_lp_train = source_lp_train[train_idx],
        calibrated_source = TRUE,
        cont_covariates = cont_covariates,
        use_interaction = use_interaction,
        include_binary = include_binary,
        eps_G = eps_G,
        w_cap = w_cap,
        gamma = gamma_target,
        select = select_target
      )

      lp_no <- predict_target_gam_lp(
        target_gam_obj = target_no_cal,
        newdata = target_train[valid_idx, , drop = FALSE],
        source_lp_new = source_lp_train[valid_idx]
      )
      lp_cal <- predict_target_gam_lp(
        target_gam_obj = target_cal,
        newdata = target_train[valid_idx, , drop = FALSE],
        source_lp_new = source_lp_train[valid_idx]
      )

      vapply(partial_grid, function(w) {
        lp_w <- (1 - w) * lp_no + w * lp_cal
        score_binary_predictions(
          y = aux_full$D_t[valid_idx],
          p = stats::plogis(lp_w),
          weights = aux_full$weights[valid_idx],
          metric = auto_metric
        )
      }, numeric(1))
    }, error = function(e) rep(failed_score, length(partial_grid)))
  }

  score_mean <- apply(score_mat, 2, mean_score)
  score_se <- apply(score_mat, 2, se_score)
  best_idx <- if (higher_is_better) {
    which.max(score_mean)
  } else {
    which.min(score_mean)
  }
  if (length(best_idx) == 0 || !is.finite(score_mean[best_idx])) {
    best_idx <- which(partial_grid == 0)[1]
    if (is.na(best_idx)) best_idx <- 1
  }

  list(
    selected_weight = partial_grid[best_idx],
    score_table = data.frame(
      partial_weight = partial_grid,
      score_mean = as.numeric(score_mean),
      score_se = as.numeric(score_se),
      auto_metric = auto_metric
    ),
    fold_scores = score_mat
  )
}

fit_trans_gam_gam_partial_calibration <- function(source_data, target_train, t,
                                                  covariates = paste0("X", 1:7),
                                                  cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                                  use_interaction = FALSE, include_binary = TRUE,
                                                  use_tcs_weight = TRUE, lambda_w = 1,
                                                  eps_G = 1e-3, w_cap = 20,
                                                  gamma_source = 1.4, gamma_target = 1.4,
                                                  select_source = FALSE, select_target = FALSE,
                                                  auto_folds = 3,
                                                  auto_metric = c("logloss", "brier", "auc"),
                                                  partial_grid = getOption("translogistic.partial_grid", c(0, 0.25, 0.5, 0.75, 1)),
                                                  auto_seed = 202605,
                                                  method_name = NULL) {
  auto_metric <- match.arg(auto_metric)
  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)

  aux_full <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  n_event <- sum(aux_full$D_t == 1, na.rm = TRUE)
  n_nonevent <- sum(aux_full$D_t == 0, na.rm = TRUE)
  nfolds <- min(auto_folds, n_event, n_nonevent)
  partial_choice <- list(
    selected_weight = 0,
    score_table = data.frame(
      partial_weight = sort(unique(partial_grid)),
      score_mean = NA_real_,
      score_se = NA_real_,
      auto_metric = auto_metric
    )
  )

  if (is.finite(nfolds) && nfolds >= 2) {
    fold_id <- make_stratified_fold_id(aux_full$D_t, nfolds = nfolds, seed = auto_seed)
    partial_choice <- cv_select_partial_gam_calibration(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      partial_grid = partial_grid,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
  }

  target_no_cal_obj <- fit_target_gam_correction(
    target_train = target_train,
    t = t,
    source_lp_train = source_lp_train,
    calibrated_source = FALSE,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_target,
    select = select_target
  )
  target_cal_obj <- fit_target_gam_correction(
    target_train = target_train,
    t = t,
    source_lp_train = source_lp_train,
    calibrated_source = TRUE,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_target,
    select = select_target
  )

  w <- partial_choice$selected_weight
  alpha_cal <- extract_source_alpha(target_cal_obj)
  source_alpha <- if (is.finite(alpha_cal)) (1 - w) * 1 + w * alpha_cal else NA_real_
  score_table <- partial_choice$score_table
  score_no_cal <- score_table$score_mean[which.min(abs(score_table$partial_weight - 0))]
  score_cal <- score_table$score_mean[which.min(abs(score_table$partial_weight - 1))]
  if (is.null(method_name)) {
    method_name <- if (use_tcs_weight) "trans_gam_gam_partial_cal" else "trans_gam_gam_partial_cal_no_tilt"
  }

  obj <- list(
    method = method_name,
    source_obj = source_obj,
    target_no_cal_obj = target_no_cal_obj,
    target_cal_obj = target_cal_obj,
    target_obj = target_cal_obj,
    partial_calibration = TRUE,
    partial_weight = w,
    partial_score_table = score_table,
    selected_target = "gam",
    selected_calibration = "partial_cal",
    auto_metric = auto_metric,
    auto_score_no_cal = if (length(score_no_cal) > 0) score_no_cal[1] else NA_real_,
    auto_score_cal = if (length(score_cal) > 0) score_cal[1] else NA_real_,
    auto_cal_loss_diff = if (length(score_no_cal) > 0 && length(score_cal) > 0) score_cal[1] - score_no_cal[1] else NA_real_,
    auto_cal_loss_diff_se = NA_real_,
    covariates = covariates,
    cont_covariates = cont_covariates,
    source_alpha = source_alpha
  )
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_gam_gam_penalized_calibration <- function(source_data, target_train, t,
                                                    covariates = paste0("X", 1:7),
                                                    cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                                    use_interaction = FALSE, include_binary = TRUE,
                                                    use_tcs_weight = TRUE, lambda_w = 1,
                                                    eps_G = 1e-3, w_cap = 20,
                                                    gamma_source = 1.4, gamma_target = 1.4,
                                                    select_source = FALSE, select_target = FALSE,
                                                    alpha_penalty = getOption("translogistic.alpha_penalty", 10),
                                                    method_name = NULL) {
  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)
  target_obj <- fit_target_gam_penalized_slope_correction(
    target_train = target_train,
    t = t,
    source_lp_train = source_lp_train,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_target,
    select = select_target,
    alpha_penalty = alpha_penalty
  )
  if (is.null(method_name)) {
    method_name <- if (use_tcs_weight) "trans_gam_gam_penalized_cal" else "trans_gam_gam_penalized_cal_no_tilt"
  }
  obj <- list(
    method = method_name,
    source_obj = source_obj,
    target_obj = target_obj,
    selected_target = "gam",
    selected_calibration = "penalized_cal",
    alpha_penalty = alpha_penalty,
    alpha_penalty_estimated = target_obj$alpha_penalty_estimated,
    covariates = covariates,
    cont_covariates = cont_covariates,
    source_alpha = extract_source_alpha(target_obj)
  )
  class(obj) <- "trans_transfer_model"
  obj
}

make_stratified_fold_id <- function(y, nfolds = 3, seed = NULL) {
  y <- as.integer(y)
  n <- length(y)
  fold_id <- rep(NA_integer_, n)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  for (yy in sort(unique(y))) {
    idx <- which(y == yy)
    idx <- sample(idx, length(idx))
    fold_id[idx] <- rep(seq_len(nfolds), length.out = length(idx))
  }

  fold_id
}

weighted_auc_binary <- function(y, score, weights = NULL) {
  y <- as.integer(y)
  score <- as.numeric(score)
  if (is.null(weights)) weights <- rep(1, length(y))
  weights <- as.numeric(weights)
  keep <- y %in% c(0L, 1L) & is.finite(score) & is.finite(weights) & weights > 0
  y <- y[keep]
  score <- score[keep]
  weights <- weights[keep]

  case_idx <- which(y == 1L)
  ctrl_idx <- which(y == 0L)
  if (length(case_idx) == 0 || length(ctrl_idx) == 0) return(NA_real_)

  case_score <- score[case_idx]
  ctrl_score <- score[ctrl_idx]
  case_weight <- weights[case_idx]
  ctrl_weight <- weights[ctrl_idx]

  cmp <- outer(case_score, ctrl_score, FUN = function(a, b) {
    as.numeric(a > b) + 0.5 * as.numeric(a == b)
  })
  w_pair <- outer(case_weight, ctrl_weight, FUN = "*")
  denom <- sum(case_weight) * sum(ctrl_weight)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(w_pair * cmp, na.rm = TRUE) / denom
}

score_binary_predictions <- function(y, p, weights = NULL,
                                     metric = c("brier", "logloss", "auc"),
                                     eps = 1e-5) {
  metric <- match.arg(metric)
  y <- as.integer(y)
  p <- as.numeric(p)
  if (is.null(weights)) weights <- rep(1, length(y))
  weights <- as.numeric(weights)
  weights[!is.finite(weights)] <- 0

  if (sum(weights) <= 0) return(if (metric == "auc") NA_real_ else Inf)
  if (metric == "auc") return(weighted_auc_binary(y = y, score = p, weights = weights))

  p <- clip_probs(p, eps = eps)

  loss <- if (metric == "brier") {
    (y - p)^2
  } else {
    -(y * log(p) + (1 - y) * log(1 - p))
  }

  sum(weights * loss, na.rm = TRUE) / sum(weights, na.rm = TRUE)
}

fit_target_correction_from_source_lp <- function(target_model,
                                                 source_lp_train,
                                                 target_train,
                                                 t,
                                                 covariates = paste0("X", 1:7),
                                                 cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                                 calibrated_source = FALSE,
                                                 use_interaction = FALSE,
                                                 include_binary = TRUE,
                                                 eps_G = 1e-3,
                                                 w_cap = 20,
                                                 max_folds = 5,
                                                 alpha_eta = 1,
                                                 lambda_choice = "lambda.1se",
                                                 gamma_target = 1.4,
                                                 select_target = FALSE) {
  target_model <- match.arg(target_model, c("linear", "gam"))
  if (target_model == "gam") {
    fit_target_gam_correction(
      target_train = target_train,
      t = t,
      source_lp_train = source_lp_train,
      calibrated_source = calibrated_source,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma = gamma_target,
      select = select_target
    )
  } else {
    estimate_linear_correction_from_source_lp(
      source_lp_train = source_lp_train,
      target_data = target_train,
      t = t,
      covariates = covariates,
      calibrated_source = calibrated_source,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      standardize_eta = TRUE
    )
  }
}

predict_target_correction_lp <- function(target_obj,
                                         target_model,
                                         newdata,
                                         source_lp_new) {
  target_model <- match.arg(target_model, c("linear", "gam"))
  if (target_model == "gam") {
    predict_target_gam_lp(
      target_gam_obj = target_obj,
      newdata = newdata,
      source_lp_new = source_lp_new
    )
  } else {
    predict(
      target_obj,
      newdata = newdata,
      source_lp_new = source_lp_new,
      type = "link"
    )
  }
}

cv_score_target_correction <- function(source_lp_train,
                                       target_train,
                                       t,
                                       fold_id,
                                       target_model = c("linear", "gam"),
                                       calibrated_source = FALSE,
                                       covariates = paste0("X", 1:7),
                                       cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                       use_interaction = FALSE,
                                       include_binary = TRUE,
                                       eps_G = 1e-3,
                                       w_cap = 20,
                                       max_folds = 5,
                                       alpha_eta = 1,
                                       lambda_choice = "lambda.1se",
                                       gamma_target = 1.4,
                                       select_target = FALSE,
                                       auto_metric = c("logloss", "brier", "auc")) {
  target_model <- match.arg(target_model)
  auto_metric <- match.arg(auto_metric)
  higher_is_better <- identical(auto_metric, "auc")
  failed_score <- if (higher_is_better) -Inf else Inf
  aux_full <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  folds <- sort(unique(fold_id[is.finite(fold_id)]))
  fold_scores <- rep(NA_real_, length(folds))

  for (i in seq_along(folds)) {
    fold <- folds[i]
    valid_idx <- which(fold_id == fold)
    train_idx <- which(fold_id != fold)

    fold_scores[i] <- tryCatch({
      target_obj <- fit_target_correction_from_source_lp(
        target_model = target_model,
        source_lp_train = source_lp_train[train_idx],
        target_train = target_train[train_idx, , drop = FALSE],
        t = t,
        covariates = covariates,
        cont_covariates = cont_covariates,
        calibrated_source = calibrated_source,
        use_interaction = use_interaction,
        include_binary = include_binary,
        eps_G = eps_G,
        w_cap = w_cap,
        max_folds = min(3, max_folds),
        alpha_eta = alpha_eta,
        lambda_choice = lambda_choice,
        gamma_target = gamma_target,
        select_target = select_target
      )
      lp_valid <- predict_target_correction_lp(
        target_obj = target_obj,
        target_model = target_model,
        newdata = target_train[valid_idx, , drop = FALSE],
        source_lp_new = source_lp_train[valid_idx]
      )
      score_binary_predictions(
        y = aux_full$D_t[valid_idx],
        p = stats::plogis(lp_valid),
        weights = aux_full$weights[valid_idx],
        metric = auto_metric
      )
    }, error = function(e) failed_score)
  }

  fold_scores
}

mean_score <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

se_score <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

choose_gam_by_one_se <- function(linear_scores,
                                 gam_scores,
                                 auto_metric = c("logloss", "brier", "auc"),
                                 tolerance = 1e-3,
                                 se_multiplier = 1) {
  auto_metric <- match.arg(auto_metric)
  higher_is_better <- identical(auto_metric, "auc")
  if (higher_is_better) {
    diff <- gam_scores - linear_scores
  } else {
    diff <- linear_scores - gam_scores
  }
  gain <- mean_score(diff)
  gain_se <- se_score(diff)
  threshold <- max(tolerance, ifelse(is.finite(gain_se), se_multiplier * gain_se, 0))
  selected <- if (is.finite(gain) && gain > threshold) "gam" else "linear"
  list(selected_target = selected, gain = gain, gain_se = gain_se, threshold = threshold)
}

choose_calibration_by_one_se <- function(no_cal_scores,
                                         cal_scores,
                                         auto_metric = c("logloss", "brier", "auc"),
                                         tolerance = 1e-3,
                                         se_multiplier = 1) {
  auto_metric <- match.arg(auto_metric)
  higher_is_better <- identical(auto_metric, "auc")
  if (higher_is_better) {
    diff <- cal_scores - no_cal_scores
    selected <- "cal"
  } else {
    diff <- cal_scores - no_cal_scores
    selected <- "cal"
  }
  diff_mean <- mean_score(diff)
  diff_se <- se_score(diff)
  threshold <- max(tolerance, ifelse(is.finite(diff_se), se_multiplier * diff_se, 0))

  use_cal <- if (higher_is_better) {
    is.finite(diff_mean) && diff_mean > threshold
  } else {
    is.finite(diff_mean) && diff_mean < -threshold
  }

  list(
    selected_calibration = if (use_cal) selected else "no_cal",
    diff = diff_mean,
    diff_se = diff_se,
    threshold = threshold
  )
}

make_auto_scores_table <- function(source_score = NA_real_,
                                   linear_score = NA_real_,
                                   gam_score = NA_real_,
                                   source_se = NA_real_,
                                   linear_se = NA_real_,
                                   gam_se = NA_real_,
                                   nfolds = NA_integer_) {
  data.frame(
    target_model = c("source", "linear", "gam"),
    score_mean = c(source_score, linear_score, gam_score),
    score_sd = c(NA_real_, NA_real_, NA_real_),
    score_se = c(source_se, linear_se, gam_se),
    nfolds = nfolds,
    stringsAsFactors = FALSE
  )
}

fit_trans_gam_auto_calibrated <- function(source_data, target_train, t,
                                          covariates = paste0("X", 1:7),
                                          cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                          use_interaction = FALSE, include_binary = TRUE,
                                          use_tcs_weight = FALSE, lambda_w = 1,
                                          calibrated_source = TRUE,
                                          eps_G = 1e-3, w_cap = 20,
                                          max_folds = 5, alpha_eta = 1,
                                          lambda_choice = "lambda.1se",
                                          gamma_source = 1.4, gamma_target = 1.4,
                                          select_source = FALSE, select_target = FALSE,
                                          auto_folds = 3,
                                          auto_metric = c("logloss", "brier", "auc"),
                                          auto_tolerance = 3e-3,
                                          auto_shift_tolerance = 3e-3,
                                          auto_seed = 202605) {
  auto_metric <- match.arg(auto_metric)
  higher_is_better <- identical(auto_metric, "auc")
  failed_score <- if (higher_is_better) -Inf else Inf

  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)

  aux_full <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  n_event <- sum(aux_full$D_t == 1, na.rm = TRUE)
  n_nonevent <- sum(aux_full$D_t == 0, na.rm = TRUE)
  nfolds <- min(auto_folds, max_folds, n_event, n_nonevent)

  cv_scores <- data.frame(
    target_model = c("source", "linear", "gam"),
    score_mean = NA_real_,
    score_sd = NA_real_,
    nfolds = nfolds,
    stringsAsFactors = FALSE
  )

  if (is.finite(nfolds) && nfolds >= 2) {
    fold_id <- make_stratified_fold_id(aux_full$D_t, nfolds = nfolds, seed = auto_seed)
    fold_scores_source <- rep(NA_real_, nfolds)
    fold_scores_linear <- rep(NA_real_, nfolds)
    fold_scores_gam <- rep(NA_real_, nfolds)

    for (fold in seq_len(nfolds)) {
      valid_idx <- which(fold_id == fold)
      train_idx <- which(fold_id != fold)

      fold_scores_source[fold] <- tryCatch({
        aux_train <- make_Dt_weights(
          target_train[train_idx, , drop = FALSE],
          t,
          eps_G = eps_G,
          w_cap = w_cap
        )
        df_train <- target_train[train_idx, , drop = FALSE]
        df_train$D_t <- aux_train$D_t
        df_train$w_ipcw <- aux_train$weights
        df_train$source_lp <- source_lp_train[train_idx]

        df_valid <- target_train[valid_idx, , drop = FALSE]
        df_valid$source_lp <- source_lp_train[valid_idx]
        p_valid <- if (isTRUE(calibrated_source)) {
          fit_source_cal <- stats::glm(
            D_t ~ source_lp,
            data = df_train,
            family = quasibinomial(),
            weights = w_ipcw,
            control = stats::glm.control(maxit = 100)
          )
          stats::plogis(stats::predict(fit_source_cal, newdata = df_valid, type = "link"))
        } else {
          stats::plogis(df_valid$source_lp)
        }

        score_binary_predictions(
          y = aux_full$D_t[valid_idx],
          p = p_valid,
          weights = aux_full$weights[valid_idx],
          metric = auto_metric
        )
      }, error = function(e) failed_score)

      fold_scores_linear[fold] <- tryCatch({
        fit_linear <- estimate_linear_correction_from_source_lp(
          source_lp_train = source_lp_train[train_idx],
          target_data = target_train[train_idx, , drop = FALSE],
          t = t,
          covariates = covariates,
          calibrated_source = calibrated_source,
          eps_G = eps_G,
          w_cap = w_cap,
          max_folds = min(3, max_folds),
          alpha_eta = alpha_eta,
          lambda_choice = lambda_choice,
          standardize_eta = TRUE
        )
        p_valid <- predict(
          fit_linear,
          newdata = target_train[valid_idx, , drop = FALSE],
          source_lp_new = source_lp_train[valid_idx],
          type = "response"
        )
        score_binary_predictions(
          y = aux_full$D_t[valid_idx],
          p = p_valid,
          weights = aux_full$weights[valid_idx],
          metric = auto_metric
        )
      }, error = function(e) failed_score)

      fold_scores_gam[fold] <- tryCatch({
        fit_gam <- fit_target_gam_correction(
          target_train = target_train[train_idx, , drop = FALSE],
          t = t,
          source_lp_train = source_lp_train[train_idx],
          calibrated_source = calibrated_source,
          cont_covariates = cont_covariates,
          use_interaction = use_interaction,
          include_binary = include_binary,
          eps_G = eps_G,
          w_cap = w_cap,
          gamma = gamma_target,
          select = select_target
        )
        lp_valid <- predict_target_gam_lp(
          target_gam_obj = fit_gam,
          newdata = target_train[valid_idx, , drop = FALSE],
          source_lp_new = source_lp_train[valid_idx]
        )
        score_binary_predictions(
          y = aux_full$D_t[valid_idx],
          p = stats::plogis(lp_valid),
          weights = aux_full$weights[valid_idx],
          metric = auto_metric
        )
      }, error = function(e) failed_score)
    }

    cv_scores$score_mean <- c(mean(fold_scores_source, na.rm = TRUE),
                              mean(fold_scores_linear, na.rm = TRUE),
                              mean(fold_scores_gam, na.rm = TRUE))
    cv_scores$score_sd <- c(stats::sd(fold_scores_source, na.rm = TRUE),
                            stats::sd(fold_scores_linear, na.rm = TRUE),
                            stats::sd(fold_scores_gam, na.rm = TRUE))
  }

  source_score <- cv_scores$score_mean[cv_scores$target_model == "source"]
  linear_score <- cv_scores$score_mean[cv_scores$target_model == "linear"]
  gam_score <- cv_scores$score_mean[cv_scores$target_model == "gam"]
  correction_score <- if (higher_is_better) {
    max(c(linear_score, gam_score), na.rm = TRUE)
  } else {
    min(c(linear_score, gam_score), na.rm = TRUE)
  }
  if (!is.finite(correction_score)) correction_score <- NA_real_
  correction_gain <- if (higher_is_better) {
    correction_score - source_score
  } else {
    source_score - correction_score
  }
  gam_gain <- if (higher_is_better) {
    gam_score - linear_score
  } else {
    linear_score - gam_score
  }

  is_large_shift <- if (is.finite(correction_gain)) {
    correction_gain > auto_shift_tolerance
  } else {
    FALSE
  }

  selected_shift <- if (isTRUE(is_large_shift)) "large" else "small"

  selected_target <- if (!isTRUE(is_large_shift)) {
    "linear"
  } else if (!higher_is_better &&
             is.finite(gam_score) &&
             (!is.finite(linear_score) ||
              gam_score + auto_tolerance < linear_score)) {
    "gam"
  } else if (higher_is_better &&
             is.finite(gam_score) &&
             (!is.finite(linear_score) ||
              gam_score > linear_score + auto_tolerance)) {
    "gam"
  } else {
    "linear"
  }

  target_obj <- if (selected_target == "gam") {
    fit_target_gam_correction(
      target_train = target_train,
      t = t,
      source_lp_train = source_lp_train,
      calibrated_source = calibrated_source,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma = gamma_target,
      select = select_target
    )
  } else {
    estimate_linear_correction_from_source_lp(
      source_lp_train = source_lp_train,
      target_data = target_train,
      t = t,
      covariates = covariates,
      calibrated_source = calibrated_source,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      standardize_eta = TRUE
    )
  }

  obj <- list(
    method = "trans_gam_auto_calibrated",
    source_obj = source_obj,
    target_obj = target_obj,
    selected_target = selected_target,
    selected_shift = selected_shift,
    selected_calibration = if (calibrated_source) "cal" else "no_cal",
    auto_scores = cv_scores,
    auto_metric = auto_metric,
    auto_calibrated_source = calibrated_source,
    auto_tolerance = auto_tolerance,
    auto_shift_tolerance = auto_shift_tolerance,
    auto_correction_gain = correction_gain,
    auto_gam_gain = gam_gain,
    covariates = covariates,
    cont_covariates = cont_covariates,
    source_alpha = extract_source_alpha(target_obj)
  )
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_gam_residual_1se <- function(source_data, target_train, t,
                                       covariates = paste0("X", 1:7),
                                       cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                       use_interaction = FALSE, include_binary = TRUE,
                                       use_tcs_weight = FALSE, lambda_w = 1,
                                       calibrated_source = FALSE,
                                       eps_G = 1e-3, w_cap = 20,
                                       max_folds = 5, alpha_eta = 1,
                                       lambda_choice = "lambda.1se",
                                       gamma_source = 1.4, gamma_target = 1.4,
                                       select_source = FALSE, select_target = FALSE,
                                       auto_folds = 3,
                                       auto_metric = c("logloss", "brier", "auc"),
                                       auto_tolerance = 1e-3,
                                       auto_se_multiplier = 1,
                                       auto_seed = 202605,
                                       method_name = NULL) {
  auto_metric <- match.arg(auto_metric)
  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)

  aux_full <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  n_event <- sum(aux_full$D_t == 1, na.rm = TRUE)
  n_nonevent <- sum(aux_full$D_t == 0, na.rm = TRUE)
  nfolds <- min(auto_folds, max_folds, n_event, n_nonevent)

  linear_scores <- gam_scores <- rep(NA_real_, nfolds)
  if (is.finite(nfolds) && nfolds >= 2) {
    fold_id <- make_stratified_fold_id(aux_full$D_t, nfolds = nfolds, seed = auto_seed)
    linear_scores <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = "linear",
      calibrated_source = calibrated_source,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
    gam_scores <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = "gam",
      calibrated_source = calibrated_source,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
  }

  target_choice <- choose_gam_by_one_se(
    linear_scores = linear_scores,
    gam_scores = gam_scores,
    auto_metric = auto_metric,
    tolerance = auto_tolerance,
    se_multiplier = auto_se_multiplier
  )
  selected_target <- target_choice$selected_target

  target_obj <- fit_target_correction_from_source_lp(
    target_model = selected_target,
    source_lp_train = source_lp_train,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    calibrated_source = calibrated_source,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    max_folds = max_folds,
    alpha_eta = alpha_eta,
    lambda_choice = lambda_choice,
    gamma_target = gamma_target,
    select_target = select_target
  )

  if (is.null(method_name)) {
    method_name <- if (calibrated_source) "trans_gam_auto_1se_calibrated" else "trans_gam_auto_1se_no_cal"
  }

  obj <- list(
    method = method_name,
    source_obj = source_obj,
    target_obj = target_obj,
    selected_target = selected_target,
    selected_shift = "correction",
    selected_calibration = if (calibrated_source) "cal" else "no_cal",
    auto_scores = make_auto_scores_table(
      linear_score = mean_score(linear_scores),
      gam_score = mean_score(gam_scores),
      linear_se = se_score(linear_scores),
      gam_se = se_score(gam_scores),
      nfolds = nfolds
    ),
    auto_metric = auto_metric,
    auto_calibrated_source = calibrated_source,
    auto_tolerance = auto_tolerance,
    auto_se_multiplier = auto_se_multiplier,
    auto_gam_gain = target_choice$gain,
    auto_gam_gain_se = target_choice$gain_se,
    auto_gam_gain_threshold = target_choice$threshold,
    covariates = covariates,
    cont_covariates = cont_covariates,
    source_alpha = extract_source_alpha(target_obj)
  )
  if (!calibrated_source) obj$source_alpha <- 1
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_gam_auto_calibration_fixed_target <- function(source_data, target_train, t,
                                                       target_model = c("linear", "gam"),
                                                       covariates = paste0("X", 1:7),
                                                       cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                                       use_interaction = FALSE, include_binary = TRUE,
                                                       use_tcs_weight = FALSE, lambda_w = 1,
                                                       eps_G = 1e-3, w_cap = 20,
                                                       max_folds = 5, alpha_eta = 1,
                                                       lambda_choice = "lambda.1se",
                                                       gamma_source = 1.4, gamma_target = 1.4,
                                                       select_source = FALSE, select_target = FALSE,
                                                       auto_folds = 3,
                                                       auto_metric = c("logloss", "brier", "auc"),
                                                       auto_tolerance = 1e-3,
                                                       auto_se_multiplier = 1,
                                                       auto_seed = 202605,
                                                       method_name = NULL) {
  target_model <- match.arg(target_model)
  auto_metric <- match.arg(auto_metric)
  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)

  aux_full <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  n_event <- sum(aux_full$D_t == 1, na.rm = TRUE)
  n_nonevent <- sum(aux_full$D_t == 0, na.rm = TRUE)
  nfolds <- min(auto_folds, max_folds, n_event, n_nonevent)

  no_cal_scores <- cal_scores <- rep(NA_real_, nfolds)
  if (is.finite(nfolds) && nfolds >= 2) {
    fold_id <- make_stratified_fold_id(aux_full$D_t, nfolds = nfolds, seed = auto_seed)
    no_cal_scores <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = target_model,
      calibrated_source = FALSE,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
    cal_scores <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = target_model,
      calibrated_source = TRUE,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
  }

  cal_choice <- choose_calibration_by_one_se(
    no_cal_scores = no_cal_scores,
    cal_scores = cal_scores,
    auto_metric = auto_metric,
    tolerance = auto_tolerance,
    se_multiplier = auto_se_multiplier
  )
  calibrated_source <- identical(cal_choice$selected_calibration, "cal")
  target_obj <- fit_target_correction_from_source_lp(
    target_model = target_model,
    source_lp_train = source_lp_train,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    calibrated_source = calibrated_source,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    max_folds = max_folds,
    alpha_eta = alpha_eta,
    lambda_choice = lambda_choice,
    gamma_target = gamma_target,
    select_target = select_target
  )

  if (is.null(method_name)) {
    method_name <- paste0("trans_gam_", target_model, "_auto_cal")
  }

  obj <- list(
    method = method_name,
    source_obj = source_obj,
    target_obj = target_obj,
    selected_target = target_model,
    selected_shift = "correction",
    selected_calibration = cal_choice$selected_calibration,
    auto_scores = make_auto_scores_table(
      linear_score = if (target_model == "linear") mean_score(no_cal_scores) else NA_real_,
      gam_score = if (target_model == "gam") mean_score(no_cal_scores) else NA_real_,
      linear_se = if (target_model == "linear") se_score(no_cal_scores) else NA_real_,
      gam_se = if (target_model == "gam") se_score(no_cal_scores) else NA_real_,
      nfolds = nfolds
    ),
    auto_score_no_cal = mean_score(no_cal_scores),
    auto_score_cal = mean_score(cal_scores),
    auto_cal_loss_diff = cal_choice$diff,
    auto_cal_loss_diff_se = cal_choice$diff_se,
    auto_cal_threshold = cal_choice$threshold,
    auto_metric = auto_metric,
    auto_calibrated_source = calibrated_source,
    auto_tolerance = auto_tolerance,
    auto_se_multiplier = auto_se_multiplier,
    covariates = covariates,
    cont_covariates = cont_covariates,
    source_alpha = extract_source_alpha(target_obj)
  )
  if (!calibrated_source) obj$source_alpha <- 1
  class(obj) <- "trans_transfer_model"
  obj
}

fit_trans_gam_auto_two_layer <- function(source_data, target_train, t,
                                         covariates = paste0("X", 1:7),
                                         cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                         use_interaction = FALSE, include_binary = TRUE,
                                         use_tcs_weight = FALSE, lambda_w = 1,
                                         eps_G = 1e-3, w_cap = 20,
                                         max_folds = 5, alpha_eta = 1,
                                         lambda_choice = "lambda.1se",
                                         gamma_source = 1.4, gamma_target = 1.4,
                                         select_source = FALSE, select_target = FALSE,
                                         auto_folds = 3,
                                         auto_metric = c("logloss", "brier", "auc"),
                                         auto_tolerance = 1e-3,
                                         auto_se_multiplier = 1,
                                         auto_seed = 202605,
                                         method_name = "trans_gam_auto_two_layer") {
  auto_metric <- match.arg(auto_metric)
  source_obj <- fit_source_gam_score(
    source_data = source_data,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    use_tcs_weight = use_tcs_weight,
    lambda_w = lambda_w,
    eps_G = eps_G,
    w_cap = w_cap,
    gamma = gamma_source,
    select = select_source
  )
  source_lp_train <- predict_source_gam_lp(source_obj, target_train)

  aux_full <- make_Dt_weights(target_train, t, eps_G = eps_G, w_cap = w_cap)
  n_event <- sum(aux_full$D_t == 1, na.rm = TRUE)
  n_nonevent <- sum(aux_full$D_t == 0, na.rm = TRUE)
  nfolds <- min(auto_folds, max_folds, n_event, n_nonevent)

  linear_no_cal <- gam_no_cal <- no_cal_scores <- cal_scores <- rep(NA_real_, nfolds)
  if (is.finite(nfolds) && nfolds >= 2) {
    fold_id <- make_stratified_fold_id(aux_full$D_t, nfolds = nfolds, seed = auto_seed)
    linear_no_cal <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = "linear",
      calibrated_source = FALSE,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
    gam_no_cal <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = "gam",
      calibrated_source = FALSE,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
  }

  target_choice <- choose_gam_by_one_se(
    linear_scores = linear_no_cal,
    gam_scores = gam_no_cal,
    auto_metric = auto_metric,
    tolerance = auto_tolerance,
    se_multiplier = auto_se_multiplier
  )
  selected_target <- target_choice$selected_target

  if (is.finite(nfolds) && nfolds >= 2) {
    fold_id <- make_stratified_fold_id(aux_full$D_t, nfolds = nfolds, seed = auto_seed + 17)
    no_cal_scores <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = selected_target,
      calibrated_source = FALSE,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
    cal_scores <- cv_score_target_correction(
      source_lp_train = source_lp_train,
      target_train = target_train,
      t = t,
      fold_id = fold_id,
      target_model = selected_target,
      calibrated_source = TRUE,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_target = gamma_target,
      select_target = select_target,
      auto_metric = auto_metric
    )
  }

  cal_choice <- choose_calibration_by_one_se(
    no_cal_scores = no_cal_scores,
    cal_scores = cal_scores,
    auto_metric = auto_metric,
    tolerance = auto_tolerance,
    se_multiplier = auto_se_multiplier
  )
  calibrated_source <- identical(cal_choice$selected_calibration, "cal")

  target_obj <- fit_target_correction_from_source_lp(
    target_model = selected_target,
    source_lp_train = source_lp_train,
    target_train = target_train,
    t = t,
    covariates = covariates,
    cont_covariates = cont_covariates,
    calibrated_source = calibrated_source,
    use_interaction = use_interaction,
    include_binary = include_binary,
    eps_G = eps_G,
    w_cap = w_cap,
    max_folds = max_folds,
    alpha_eta = alpha_eta,
    lambda_choice = lambda_choice,
    gamma_target = gamma_target,
    select_target = select_target
  )

  obj <- list(
    method = method_name,
    source_obj = source_obj,
    target_obj = target_obj,
    selected_target = selected_target,
    selected_shift = "correction",
    selected_calibration = cal_choice$selected_calibration,
    auto_scores = make_auto_scores_table(
      linear_score = mean_score(linear_no_cal),
      gam_score = mean_score(gam_no_cal),
      linear_se = se_score(linear_no_cal),
      gam_se = se_score(gam_no_cal),
      nfolds = nfolds
    ),
    auto_score_no_cal = mean_score(no_cal_scores),
    auto_score_cal = mean_score(cal_scores),
    auto_cal_loss_diff = cal_choice$diff,
    auto_cal_loss_diff_se = cal_choice$diff_se,
    auto_cal_threshold = cal_choice$threshold,
    auto_metric = auto_metric,
    auto_calibrated_source = calibrated_source,
    auto_tolerance = auto_tolerance,
    auto_se_multiplier = auto_se_multiplier,
    auto_gam_gain = target_choice$gain,
    auto_gam_gain_se = target_choice$gain_se,
    auto_gam_gain_threshold = target_choice$threshold,
    covariates = covariates,
    cont_covariates = cont_covariates,
    source_alpha = extract_source_alpha(target_obj)
  )
  if (!calibrated_source) obj$source_alpha <- 1
  class(obj) <- "trans_transfer_model"
  obj
}

predict.trans_transfer_model <- function(object, newdata, type = c("link", "response"), ...) {
  type <- match.arg(type)
  method <- object$method

  if (method == "trans_linear_linear") {
    lp <- linear_lp_from_beta(object$beta, newdata, object$covariates)
  } else if (method %in% c("trans_linear_gam", "trans_linear_gam_calibrated")) {
    source_lp <- linear_lp_from_beta(object$source_obj$beta, newdata, object$covariates)
    lp <- predict_target_gam_lp(object$target_obj, newdata, source_lp)
  } else if (method %in% c("trans_gam_linear", "trans_gam_linear_calibrated")) {
    source_lp <- predict_source_gam_lp(object$source_obj, newdata)
    lp <- predict(object$target_obj, newdata = newdata, source_lp_new = source_lp, type = "link")
  } else if (method %in% c("trans_gam_gam", "trans_gam_gam_tp",
                           "trans_gam_gam_source_tp", "trans_gam_gam_all_tp",
                           "trans_gam_gam_select",
                           "trans_gam_gam_calibrated",
                           "trans_gam_gam_no_tilt", "trans_gam_gam_tp_no_tilt",
                           "trans_gam_gam_source_tp_no_tilt", "trans_gam_gam_all_tp_no_tilt",
                           "trans_gam_gam_select_no_tilt",
                           "trans_gam_gam_calibrated_no_tilt",
                           "trans_gam_gam_penalized_cal", "trans_gam_gam_penalized_cal_no_tilt")) {
    source_lp <- predict_source_gam_lp(object$source_obj, newdata)
    lp <- predict_target_gam_lp(object$target_obj, newdata, source_lp)
  } else if (method %in% c("trans_gam_gam_partial_cal", "trans_gam_gam_partial_cal_no_tilt")) {
    source_lp <- predict_source_gam_lp(object$source_obj, newdata)
    lp_no <- predict_target_gam_lp(object$target_no_cal_obj, newdata, source_lp)
    lp_cal <- predict_target_gam_lp(object$target_cal_obj, newdata, source_lp)
    w <- if (is.null(object$partial_weight)) 0.5 else object$partial_weight
    lp <- (1 - w) * lp_no + w * lp_cal
  } else if (method %in% c("trans_gam_auto_calibrated",
                           "trans_gam_auto_calibrated_tilt",
                           "trans_gam_auto_auc",
                           "trans_gam_auto_auc_tilt",
                           "trans_gam_auto_no_cal",
                           "trans_gam_auto_no_cal_tilt",
                           "trans_gam_auto_auc_no_cal",
                           "trans_gam_auto_auc_no_cal_tilt",
                           "trans_gam_auto_1se_no_cal",
                           "trans_gam_auto_1se_calibrated",
                           "trans_gam_auto_two_layer",
                           "trans_gam_linear_auto_cal",
                           "trans_gam_gam_auto_cal")) {
    source_lp <- predict_source_gam_lp(object$source_obj, newdata)
    if (identical(object$selected_target, "gam")) {
      lp <- predict_target_gam_lp(object$target_obj, newdata, source_lp)
    } else {
      lp <- predict(object$target_obj, newdata = newdata, source_lp_new = source_lp, type = "link")
    }
  } else {
    stop("Unknown trans_transfer_model method: ", method)
  }

  if (type == "link") return(lp)
  stats::plogis(lp)
}

predict_transfer_risk_t <- function(object, newdata) {
  lp <- predict(object, newdata = newdata, type = "link")
  list(lp = lp, p_event_t = stats::plogis(lp))
}

fit_rsf_survival_model <- function(train_data,
                                   covariates = paste0("X", 1:7),
                                   num.trees = as.integer(getOption("translogistic.rsf_trees", 300)),
                                   min.node.size = as.integer(getOption("translogistic.rsf_min_node_size", 20)),
                                   num.threads = as.integer(getOption("translogistic.rsf_threads", 1)),
                                   seed = 202605) {
  .require_pkg("ranger")

  df <- train_data[, c("Y", "delta", covariates), drop = FALSE]
  fml <- stats::as.formula(
    paste("Surv(Y, delta) ~", paste(covariates, collapse = " + "))
  )
  mtry_use <- max(1, floor(sqrt(length(covariates))))

  fit <- ranger::ranger(
    formula = fml,
    data = df,
    num.trees = num.trees,
    mtry = mtry_use,
    min.node.size = min.node.size,
    splitrule = "logrank",
    write.forest = TRUE,
    num.threads = max(1L, num.threads),
    seed = seed
  )

  list(
    fit = fit,
    covariates = covariates,
    model_type = "rsf",
    num.trees = num.trees,
    min.node.size = min.node.size,
    num.threads = max(1L, num.threads),
    seed = seed
  )
}

predict_rsf_risk_t <- function(object, newdata, t) {
  pred <- stats::predict(
    object$fit,
    data = newdata[, object$covariates, drop = FALSE]
  )
  event_times <- as.numeric(pred$unique.death.times)
  surv_mat <- as.matrix(pred$survival)

  if (length(event_times) == 0 || ncol(surv_mat) == 0) {
    p_event_t <- rep(0, nrow(newdata))
  } else {
    time_idx <- findInterval(t, event_times)
    surv_t <- if (time_idx == 0) {
      rep(1, nrow(newdata))
    } else {
      surv_mat[, time_idx]
    }
    p_event_t <- 1 - as.numeric(surv_t)
  }

  p_event_t <- clip_probs(p_event_t, eps = 1e-5)
  list(lp = stats::qlogis(p_event_t), p_event_t = p_event_t)
}

fit_ipcw_xgb_model <- function(train_data,
                               t,
                               covariates = paste0("X", 1:7),
                               eps_G = 1e-3,
                               w_cap = 20,
                               nrounds = 150,
                               max_depth = 2,
                               eta = 0.05,
                               subsample = 0.8,
                               colsample_bytree = 0.8,
                               min_child_weight = 5,
                               nthread = as.integer(getOption("translogistic.xgb_threads", 1)),
                               seed = 202605) {
  .require_pkg("xgboost")

  aux <- make_Dt_weights(
    data = train_data,
    t = t,
    eps_G = eps_G,
    w_cap = w_cap
  )

  x_train <- as.matrix(train_data[, covariates, drop = FALSE])
  storage.mode(x_train) <- "double"

  dtrain <- xgboost::xgb.DMatrix(
    data = x_train,
    label = aux$D_t,
    weight = aux$weights
  )

  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = max_depth,
    eta = eta,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    min_child_weight = min_child_weight,
    nthread = max(1L, nthread)
  )

  set.seed(seed)
  fit <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )

  list(
    fit = fit,
    covariates = covariates,
    model_type = "ipcw_xgb",
    t = t,
    nrounds = nrounds,
    params = params,
    n_event = sum(aux$D_t == 1, na.rm = TRUE),
    n_nonevent = sum(aux$D_t == 0, na.rm = TRUE)
  )
}

predict_ipcw_xgb_risk_t <- function(object, newdata) {
  x_new <- as.matrix(newdata[, object$covariates, drop = FALSE])
  storage.mode(x_new) <- "double"
  p_event_t <- as.numeric(stats::predict(object$fit, newdata = x_new))
  p_event_t <- clip_probs(p_event_t, eps = 1e-5)
  list(lp = stats::qlogis(p_event_t), p_event_t = p_event_t)
}

load_transcox_sources <- local({
  loaded <- FALSE
  function() {
    if (isTRUE(loaded) && exists("runTransCox_one", mode = "function", inherits = TRUE)) {
      return(invisible(TRUE))
    }

    roots <- unique(c(
      file.path(getwd(), "TransCox-master"),
      file.path(dirname(getwd()), "TransCox-master"),
      getOption("translogistic.transcox_dir", "")
    ))
    roots <- roots[nzchar(roots) & dir.exists(roots)]
    if (length(roots) == 0) {
      stop("Cannot find TransCox-master. Set options(translogistic.transcox_dir=...).")
    }

    r_dir <- file.path(roots[[1]], "R")
    files <- c(
      "zzz_utils.R",
      "dQtocumQ.R",
      "deltaQ.R",
      "GetLogLike.R",
      "GetBIC.R",
      "GetAuxSurv.R",
      "GetPrimaryParam.R",
      "runTransCox_one.R",
      "SelParam_By_BIC.R",
      "SelLR_By_BIC.R",
      "ChekLR_BIC.R"
    )
    for (f in files) {
      path <- file.path(r_dir, f)
      if (file.exists(path)) source(path)
    }
    loaded <<- TRUE
    invisible(TRUE)
  }
})

make_transcox_data <- function(data, covariates, out_covariates = paste0("X", seq_along(covariates))) {
  out <- data.frame(
    time = as.numeric(data$Y),
    status = ifelse(as.integer(data$delta) == 1L, 2L, 1L),
    stringsAsFactors = FALSE
  )
  x <- as.data.frame(data[, covariates, drop = FALSE])
  names(x) <- out_covariates
  cbind(out, x)
}

standardize_transcox_pair <- function(source_tc, target_tc, covariates) {
  pooled_x <- rbind(source_tc[, covariates, drop = FALSE], target_tc[, covariates, drop = FALSE])
  center <- vapply(pooled_x, mean, numeric(1), na.rm = TRUE)
  scale <- vapply(pooled_x, stats::sd, numeric(1), na.rm = TRUE)
  scale[!is.finite(scale) | scale <= 0] <- 1

  apply_one <- function(df) {
    for (v in covariates) df[[v]] <- (df[[v]] - center[[v]]) / scale[[v]]
    df
  }

  list(
    source = apply_one(source_tc),
    target = apply_one(target_tc),
    center = center,
    scale = scale,
    covariates = covariates
  )
}

fit_transcox_model <- function(source_data,
                               target_train,
                               source_covariates = paste0("X", 1:7),
                               target_covariates = paste0("X", 1:7),
                               select_bic = TRUE,
                               lambda1 = getOption("translogistic.transcox_lambda1", 0.1),
                               lambda2 = getOption("translogistic.transcox_lambda2", 0.1),
                               lambda1_vec = getOption("translogistic.transcox_lambda1_vec", c(0.1, 0.5, 1, 2)),
                               lambda2_vec = getOption("translogistic.transcox_lambda2_vec", c(0.1, 0.5, 1, 2)),
                               learning_rate = getOption("translogistic.transcox_learning_rate", 0.004),
                               nsteps = getOption("translogistic.transcox_nsteps", 100),
                               backend = getOption("translogistic.transcox_backend", "r"),
                               standardize = TRUE) {
  load_transcox_sources()

  cov_tc <- paste0("X", seq_along(target_covariates))
  source_tc <- make_transcox_data(source_data, source_covariates, cov_tc)
  target_tc <- make_transcox_data(target_train, target_covariates, cov_tc)

  if (isTRUE(standardize)) {
    sc <- standardize_transcox_pair(source_tc, target_tc, cov_tc)
    source_tc <- sc$source
    target_tc <- sc$target
  } else {
    sc <- list(center = setNames(rep(0, length(cov_tc)), cov_tc),
               scale = setNames(rep(1, length(cov_tc)), cov_tc),
               covariates = cov_tc)
  }

  if (sum(source_tc$status == 2) < 1 || sum(target_tc$status == 2) < 1) {
    stop("TransCox requires at least one event in both source and target training data.")
  }

  Cout <- GetAuxSurv(source_tc, cov = cov_tc)
  Pout <- GetPrimaryParam(target_tc, q = Cout$q, estR = Cout$estR)

  bic <- NULL
  if (isTRUE(select_bic)) {
    bic <- SelParam_By_BIC(
      primData = target_tc,
      auxData = source_tc,
      cov = cov_tc,
      lambda1_vec = lambda1_vec,
      lambda2_vec = lambda2_vec,
      learning_rate = learning_rate,
      nsteps = nsteps,
      backend = backend,
      verbose = FALSE
    )
    lambda1 <- bic$best_la1
    lambda2 <- bic$best_la2
  }

  fit <- runTransCox_one(
    Pout = Pout,
    l1 = lambda1,
    l2 = lambda2,
    learning_rate = learning_rate,
    nsteps = nsteps,
    backend = backend,
    cov = cov_tc
  )

  obj <- list(
    method = if (isTRUE(select_bic)) "transcox_bic" else "transcox_fixed",
    fit = fit,
    Pout = Pout,
    Cout = Cout,
    bic = bic,
    lambda1 = lambda1,
    lambda2 = lambda2,
    learning_rate = learning_rate,
    nsteps = nsteps,
    backend = backend,
    scaler = sc,
    source_covariates = source_covariates,
    target_covariates = target_covariates,
    covariates = cov_tc
  )
  class(obj) <- "transcox_model"
  obj
}

predict_transcox_risk_t <- function(object, newdata, t) {
  x <- as.data.frame(newdata[, object$target_covariates, drop = FALSE])
  names(x) <- object$covariates
  for (v in object$covariates) {
    x[[v]] <- (x[[v]] - object$scaler$center[[v]]) / object$scaler$scale[[v]]
  }

  lp <- as.vector(as.matrix(x[, object$covariates, drop = FALSE]) %*% as.numeric(object$fit$new_beta))
  event_times <- as.numeric(object$fit$time)
  increments <- pmax(as.numeric(object$fit$new_IntH), 0)
  H0_t <- sum(increments[event_times <= t], na.rm = TRUE)
  if (!is.finite(H0_t) || H0_t < 0) H0_t <- 0
  p_event_t <- 1 - exp(-H0_t * exp(lp))
  p_event_t <- clip_probs(p_event_t, eps = 1e-5)
  list(lp = lp, p_event_t = p_event_t)
}

fit_transfer_method <- function(method, source_data, target_train, t,
                                covariates = paste0("X", 1:7),
                                cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                use_interaction = FALSE, include_binary = TRUE,
                                lambda_w = 1, eps_G = 1e-3, w_cap = 20,
                                max_folds = 5, alpha_eta = 1,
                                lambda_choice = "lambda.1se",
                                gamma_source = 1.4, gamma_target = 1.4,
                                select_source = FALSE, select_target = FALSE) {
  if (method == "trans_linear_linear") {
    return(fit_trans_linear_linear(source_data, target_train, t, covariates, lambda_w,
                                   eps_G, w_cap, max_folds, alpha_eta, lambda_choice))
  }

  if (method %in% c("trans_linear_gam", "trans_linear_gam_calibrated")) {
    return(fit_trans_linear_gam(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      calibrated_source = method == "trans_linear_gam_calibrated",
      use_tcs_weight = TRUE,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma = gamma_target,
      select = select_target,
      method_name = method
    ))
  }

  if (method %in% c("trans_gam_linear", "trans_gam_linear_calibrated")) {
    return(fit_trans_gam_linear(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      calibrated_source = method == "trans_gam_linear_calibrated",
      use_tcs_weight = TRUE,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_source = gamma_source,
      select_source = select_source,
      method_name = method
    ))
  }

  if (method %in% c("trans_gam_auto_calibrated",
                    "trans_gam_auto_calibrated_tilt",
                    "trans_gam_auto_auc",
                    "trans_gam_auto_auc_tilt",
                    "trans_gam_auto_no_cal",
                    "trans_gam_auto_no_cal_tilt",
                    "trans_gam_auto_auc_no_cal",
                    "trans_gam_auto_auc_no_cal_tilt")) {
    fit_obj <- fit_trans_gam_auto_calibrated(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      use_tcs_weight = method %in% c("trans_gam_auto_calibrated_tilt",
                                     "trans_gam_auto_auc_tilt",
                                     "trans_gam_auto_no_cal_tilt",
                                     "trans_gam_auto_auc_no_cal_tilt"),
      calibrated_source = !(method %in% c("trans_gam_auto_no_cal",
                                          "trans_gam_auto_no_cal_tilt",
                                          "trans_gam_auto_auc_no_cal",
                                          "trans_gam_auto_auc_no_cal_tilt")),
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target,
      auto_metric = if (method %in% c("trans_gam_auto_auc",
                                      "trans_gam_auto_auc_tilt",
                                      "trans_gam_auto_auc_no_cal",
                                      "trans_gam_auto_auc_no_cal_tilt")) "auc" else "logloss"
    )
    fit_obj$method <- method
    return(fit_obj)
  }

  if (method %in% c("trans_gam_auto_1se_no_cal",
                    "trans_gam_auto_1se_calibrated")) {
    return(fit_trans_gam_residual_1se(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      use_tcs_weight = FALSE,
      calibrated_source = method == "trans_gam_auto_1se_calibrated",
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target,
      auto_metric = "logloss",
      method_name = method
    ))
  }

  if (method == "trans_gam_auto_two_layer") {
    return(fit_trans_gam_auto_two_layer(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      use_tcs_weight = FALSE,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target,
      auto_metric = "logloss",
      method_name = method
    ))
  }

  if (method %in% c("trans_gam_linear_auto_cal", "trans_gam_gam_auto_cal")) {
    return(fit_trans_gam_auto_calibration_fixed_target(
      source_data = source_data,
      target_train = target_train,
      t = t,
      target_model = if (method == "trans_gam_gam_auto_cal") "gam" else "linear",
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      use_tcs_weight = FALSE,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = alpha_eta,
      lambda_choice = lambda_choice,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target,
      auto_metric = "logloss",
      method_name = method
    ))
  }

  if (method %in% c("trans_gam_gam", "trans_gam_gam_tp",
                    "trans_gam_gam_source_tp", "trans_gam_gam_all_tp",
                    "trans_gam_gam_select",
                    "trans_gam_gam_calibrated",
                    "trans_gam_gam_no_tilt", "trans_gam_gam_tp_no_tilt",
                    "trans_gam_gam_source_tp_no_tilt", "trans_gam_gam_all_tp_no_tilt",
                    "trans_gam_gam_select_no_tilt",
                    "trans_gam_gam_calibrated_no_tilt")) {
    use_tcs <- !(method %in% c("trans_gam_gam_no_tilt",
                               "trans_gam_gam_tp_no_tilt",
                               "trans_gam_gam_source_tp_no_tilt",
                               "trans_gam_gam_all_tp_no_tilt",
                               "trans_gam_gam_select_no_tilt",
                               "trans_gam_gam_calibrated_no_tilt"))
    calibrated <- method %in% c("trans_gam_gam_calibrated", "trans_gam_gam_calibrated_no_tilt")
    force_target_select <- method %in% c("trans_gam_gam_select", "trans_gam_gam_select_no_tilt")
    source_bs <- if (method %in% c("trans_gam_gam_source_tp", "trans_gam_gam_source_tp_no_tilt",
                                   "trans_gam_gam_all_tp", "trans_gam_gam_all_tp_no_tilt")) "tp" else "ts"
    target_bs <- if (method %in% c("trans_gam_gam_tp", "trans_gam_gam_tp_no_tilt",
                                   "trans_gam_gam_all_tp", "trans_gam_gam_all_tp_no_tilt")) "tp" else "ts"
    return(fit_trans_gam_gam(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      calibrated_source = calibrated,
      use_tcs_weight = use_tcs,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target || force_target_select,
      source_smooth_bs = source_bs,
      target_smooth_bs = target_bs,
      method_name = method
    ))
  }

  if (method %in% c("trans_gam_gam_partial_cal", "trans_gam_gam_partial_cal_no_tilt")) {
    return(fit_trans_gam_gam_partial_calibration(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      use_tcs_weight = method == "trans_gam_gam_partial_cal",
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target,
      auto_metric = "logloss",
      method_name = method
    ))
  }

  if (method %in% c("trans_gam_gam_penalized_cal", "trans_gam_gam_penalized_cal_no_tilt")) {
    return(fit_trans_gam_gam_penalized_calibration(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      use_tcs_weight = method == "trans_gam_gam_penalized_cal",
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target,
      method_name = method
    ))
  }

  if (method %in% c("source_gam_recal_intercept", "source_gam_recal_slope")) {
    recal <- if (method == "source_gam_recal_intercept") "intercept" else "slope"
    return(fit_source_gam_recalibration(
      source_data = source_data,
      target_train = target_train,
      t = t,
      recalibration = recal,
      covariates = covariates,
      cont_covariates = cont_covariates,
      use_interaction = use_interaction,
      include_binary = include_binary,
      use_tcs_weight = TRUE,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma_source = gamma_source,
      select_source = select_source
    ))
  }

  stop("Unknown transfer method: ", method)
}


# ============================================================================== 
# 6. Target-only GAM model
# ============================================================================== 

fit_target_only_gam_model <- function(target_train,
                                      t,
                                      cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
                                      use_interaction = FALSE,
                                      include_binary = TRUE,
                                      eps_G = 1e-3,
                                      w_cap = 20,
                                      gamma = 1.4,
                                      select = FALSE) {
  .require_pkg("mgcv")

  aux <- make_Dt_weights(
    data = target_train,
    t = t,
    eps_G = eps_G,
    w_cap = w_cap
  )

  df <- target_train
  df$D_t <- aux$D_t
  df$w_ipcw <- aux$weights

  scaler <- make_scaler(df, cont_covariates)
  df <- apply_scaler(df, scaler)

  fit <- mgcv::gam(
    formula = make_source_gam_formula(
      use_interaction = use_interaction,
      include_binary = include_binary
    ),
    data = df,
    family = quasibinomial(),
    weights = w_ipcw,
    method = "REML",
    gamma = gamma,
    select = select
  )

  obj <- list(
    method = "target_only_gam",
    fit = fit,
    scaler = scaler,
    cont_covariates = cont_covariates,
    use_interaction = use_interaction,
    include_binary = include_binary,
    n_event = sum(aux$D_t == 1),
    n_nonevent = sum(aux$D_t == 0)
  )
  class(obj) <- "target_only_gam_model"
  obj
}

predict.target_only_gam_model <- function(object,
                                          newdata,
                                          type = c("link", "response"),
                                          ...) {
  type <- match.arg(type)

  newdata_sc <- apply_scaler(newdata, object$scaler)
  lp <- as.vector(stats::predict(object$fit, newdata = newdata_sc, type = "link"))

  if (type == "link") return(lp)
  stats::plogis(lp)
}

predict_target_only_gam_risk_t <- function(object,
                                           newdata) {
  lp <- predict(object, newdata = newdata, type = "link")
  list(lp = lp, p_event_t = stats::plogis(lp))
}
