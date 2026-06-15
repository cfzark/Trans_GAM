local_r_lib <- file.path(getwd(), "_r_lib")
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

library(MASS)
library(survival)
library(glmnet)

# Use the structured functions file generated in this conversation.
# If you rename it to functions.R, change this line accordingly.
source(file.path(getOption("transgam.r_dir", "."), "transfer_gam_functions.R"))
################################################################################
# sim_setup_aft_nonlinear.R
#
# Formal simulation with:
#   source_data + target_train + target_test
#
# Supports:
#   - covariate shift: moderate / strong
#   - linear or nonlinear truth
#   - AFT-logistic event time as the primary non-Weibull DGP
#   - Weibull PH event time remains supported by the generic generator
#   - log-normal independent censoring calibrated to target censoring rate
#   - non-transfer baselines:
#       source/target/pooled x logistic/GAM/RSF/IPCW-XGBoost/Cox
#   - trans_linear_linear
#   - trans_linear_gam
#   - trans_gam_gam
#   - trans_gam_gam_calibrated
#   - old alias trans_logistic -> trans_linear_linear
#   - oracle
################################################################################


################################################################################
# 1. Default simulation parameters
################################################################################

default_cov_shift7 <- function(shift_level = c("moderate", "strong")) {
  shift_level <- match.arg(shift_level)

  if (shift_level == "moderate") {
    return(list(
      p_bin_source = c(0.35, 0.45),
      p_bin_target = c(0.60, 0.75),

      rho_source   = 0.25,
      rho_target   = 0.65,

      mu34_source  = c(-0.25,  0.20),
      mu34_target  = c( 0.35, -0.25),

      x5_source    = c(mean = 0.0, sd = 1.0),
      x5_target    = c(mean = 0.7, sd = 1.0),

      x6_source    = list(dist = "unif", min = 0, max = 1),
      x6_target    = list(dist = "beta", a = 2.5, b = 2.0),

      x7_source    = c(mean = -0.3, sd = 0.8),
      x7_target    = c(mean = -0.1, sd = 1.2)
    ))
  }

  list(
    p_bin_source = c(0.25, 0.35),
    p_bin_target = c(0.70, 0.85),

    rho_source   = 0.15,
    rho_target   = 0.80,

    mu34_source  = c(-0.45,  0.30),
    mu34_target  = c( 0.55, -0.45),

    x5_source    = c(mean = -0.2, sd = 1.0),
    x5_target    = c(mean =  1.0, sd = 1.1),

    x6_source    = list(dist = "unif", min = 0, max = 1),
    x6_target    = list(dist = "beta", a = 3.0, b = 1.8),

    x7_source    = c(mean = -0.3, sd = 0.7),
    x7_target    = c(mean =  0.0, sd = 1.5)
  )
}


default_truth7 <- function(use_nonlinear_in_truth = TRUE,
                           intercept = 0,
                           beta_linear = c(0.9, -0.8, 0.7, -0.6, 0.5, 0.0, 0.0),
                           beta_nonlinear = c(I_X4_2 = 0.55,
                                              I_X1X3 = -0.45,
                                              I_sinX3 = 0.50)) {
  beta_linear <- as.numeric(beta_linear)

  if (length(beta_linear) != 7) {
    stop("beta_linear must have length 7.")
  }

  if (!use_nonlinear_in_truth) {
    beta <- c(
      Intercept = intercept,
      setNames(beta_linear, paste0("X", 1:7))
    )

    return(list(
      use_nonlinear_in_truth = FALSE,
      beta_true = beta,
      beta_linear = beta_linear,
      beta_nonlinear = numeric(0)
    ))
  }

  beta <- c(
    Intercept = intercept,
    setNames(beta_linear, paste0("X", 1:7)),
    beta_nonlinear
  )

  list(
    use_nonlinear_in_truth = TRUE,
    beta_true = beta,
    beta_linear = beta_linear,
    beta_nonlinear = beta_nonlinear
  )
}


default_time7 <- function(time_model = c("weibull", "aft_logistic"),
                          weib_shape = 1.5,
                          weib_scale = 1.0,
                          t0 = 10,
                          aft_logistic_scale = 0.7) {
  time_model <- match.arg(time_model)

  list(
    time_model = time_model,
    weib_shape = weib_shape,
    weib_scale = weib_scale,
    t0 = t0,
    aft_logistic_scale = aft_logistic_scale
  )
}


################################################################################
# 2. Censoring
################################################################################

simulate_censor7 <- function(T,
                             target_censor_rate,
                             sigma_c = 0.6,
                             tol = 0.005,
                             max_iter = 100,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- length(T)
  Z <- rnorm(n)

  f_mu <- function(mu) {
    C <- exp(mu + sigma_c * Z)
    mean(T >= C)
  }

  L <- -10
  U <- 10

  while (f_mu(L) < target_censor_rate) L <- L - 2
  while (f_mu(U) > target_censor_rate) U <- U + 2

  iter <- 0
  mu_c <- (L + U) / 2
  f_mid <- f_mu(mu_c)

  while (abs(f_mid - target_censor_rate) > tol && iter < max_iter) {
    if (f_mid > target_censor_rate) {
      L <- mu_c
    } else {
      U <- mu_c
    }

    mu_c <- (L + U) / 2
    f_mid <- f_mu(mu_c)
    iter <- iter + 1
  }

  C <- exp(mu_c + sigma_c * Z)

  list(
    C = C,
    mu_c = mu_c,
    sigma_c = sigma_c,
    achieved_censor = mean(T >= C),
    iter = iter,
    target_censor_rate = target_censor_rate
  )
}


################################################################################
# 3. Covariate and truth feature generation
################################################################################

make_design_features7 <- function(df, use_nonlinear = TRUE) {
  X1 <- df$X1
  X2 <- df$X2
  X3 <- df$X3
  X4 <- df$X4
  X5 <- df$X5
  X6 <- df$X6
  X7 <- df$X7

  if (!use_nonlinear) {
    X <- as.matrix(df[, paste0("X", 1:7), drop = FALSE])
    colnames(X) <- paste0("X", 1:7)
    return(X)
  }

  cbind(
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X4 = X4,
    X5 = X5,
    X6 = X6,
    X7 = X7,
    I_X4_2  = X4^2,
    I_X1X3  = X1 * X3,
    I_sinX3 = sin(X3)
  )
}


gen_covariates7 <- function(n,
                            dataset = c("target", "source"),
                            cov_shift = default_cov_shift7("moderate"),
                            seed = NULL) {
  dataset <- match.arg(dataset)

  if (!is.null(seed)) set.seed(seed)

  if (dataset == "target") {
    p_bin <- cov_shift$p_bin_target
    rho   <- cov_shift$rho_target
    mu34  <- cov_shift$mu34_target
    x5par <- cov_shift$x5_target
    x6par <- cov_shift$x6_target
    x7par <- cov_shift$x7_target
  } else {
    p_bin <- cov_shift$p_bin_source
    rho   <- cov_shift$rho_source
    mu34  <- cov_shift$mu34_source
    x5par <- cov_shift$x5_source
    x6par <- cov_shift$x6_source
    x7par <- cov_shift$x7_source
  }

  X1 <- rbinom(n, 1, p_bin[1])
  X2 <- rbinom(n, 1, p_bin[2])

  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  L <- chol(Sigma)
  Z <- matrix(rnorm(n * 2), n, 2) %*% L

  X3 <- Z[, 1] + mu34[1]
  X4 <- Z[, 2] + mu34[2]

  X5 <- rnorm(n, mean = x5par["mean"], sd = x5par["sd"])

  if (x6par$dist == "beta") {
    X6 <- rbeta(n, x6par$a, x6par$b)
  } else if (x6par$dist == "unif") {
    X6 <- runif(n, x6par$min, x6par$max)
  } else {
    stop("Unsupported x6 distribution.")
  }

  X7 <- rnorm(n, mean = x7par["mean"], sd = x7par["sd"])

  data.frame(
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X4 = X4,
    X5 = X5,
    X6 = X6,
    X7 = X7
  )
}


################################################################################
# 4. Event time generation
################################################################################

gen_event_time7 <- function(dfX, truth_par, time_par) {
  X_truth <- make_design_features7(
    df = dfX,
    use_nonlinear = truth_par$use_nonlinear_in_truth
  )
  X_truth <- cbind(Intercept = 1, X_truth)

  beta_true <- truth_par$beta_true

  if (length(beta_true) != ncol(X_truth)) {
    stop(sprintf(
      "beta_true length = %d but truth design has %d columns.",
      length(beta_true),
      ncol(X_truth)
    ))
  }

  eta <- as.vector(X_truth %*% beta_true)

  if (time_par$time_model == "weibull") {
    U <- runif(nrow(dfX))
    T_event <- time_par$weib_scale *
      ((-log(U)) / exp(eta))^(1 / time_par$weib_shape)
  } else if (time_par$time_model == "aft_logistic") {
    U <- runif(nrow(dfX))
    eps <- log(U / (1 - U))
    logT <- log(time_par$t0) - eta + time_par$aft_logistic_scale * eps
    T_event <- exp(logT)
  } else {
    stop("Unsupported time model.")
  }

  list(
    T_event = T_event,
    eta_true = eta,
    X_truth = X_truth
  )
}


################################################################################
# 5. Dataset generation
################################################################################

simulate_surv_data7 <- function(n,
                                dataset = c("target", "source"),
                                truth_par = default_truth7(TRUE),
                                time_par = default_time7("weibull"),
                                cov_shift = default_cov_shift7("moderate"),
                                target_censor_rate = 0.35,
                                censor_sigma = 0.6,
                                seed = NULL) {
  dataset <- match.arg(dataset)

  if (!is.null(seed)) set.seed(seed)

  dfX <- gen_covariates7(
    n = n,
    dataset = dataset,
    cov_shift = cov_shift
  )

  event_obj <- gen_event_time7(
    dfX = dfX,
    truth_par = truth_par,
    time_par = time_par
  )

  if (target_censor_rate > 0) {
    cens_obj <- simulate_censor7(
      T = event_obj$T_event,
      target_censor_rate = target_censor_rate,
      sigma_c = censor_sigma
    )
    C <- cens_obj$C
  } else {
    C <- rep(max(event_obj$T_event) + 1, n)

    cens_obj <- list(
      C = C,
      mu_c = NA_real_,
      sigma_c = NA_real_,
      achieved_censor = 0,
      iter = 0,
      target_censor_rate = 0
    )
  }

  Y <- pmin(event_obj$T_event, C)
  delta <- as.integer(event_obj$T_event <= C)

  data <- data.frame(
    T_event = event_obj$T_event,
    C = C,
    Y = Y,
    delta = delta,
    dfX
  )

  list(
    data = data,
    censor_info = cens_obj,
    eta_true = event_obj$eta_true,
    truth_par = truth_par,
    time_par = time_par,
    cov_shift = cov_shift,
    dataset = dataset
  )
}


generate_train_test_surv_data7 <- function(n_source = 1000,
                                           n_target_train = 200,
                                           n_target_test = 1000,
                                           truth_source = default_truth7(TRUE),
                                           truth_target = truth_source,
                                           time_source = default_time7("weibull"),
                                           time_target = time_source,
                                           cov_shift = default_cov_shift7("moderate"),
                                           censor_source = 0.35,
                                           censor_target_train = 0.35,
                                           censor_target_test = 0.35,
                                           censor_sigma = 0.6,
                                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  src_seed <- if (is.null(seed)) NULL else seed + 1
  trn_seed <- if (is.null(seed)) NULL else seed + 2
  tst_seed <- if (is.null(seed)) NULL else seed + 3

  src_obj <- simulate_surv_data7(
    n = n_source,
    dataset = "source",
    truth_par = truth_source,
    time_par = time_source,
    cov_shift = cov_shift,
    target_censor_rate = censor_source,
    censor_sigma = censor_sigma,
    seed = src_seed
  )

  trn_obj <- simulate_surv_data7(
    n = n_target_train,
    dataset = "target",
    truth_par = truth_target,
    time_par = time_target,
    cov_shift = cov_shift,
    target_censor_rate = censor_target_train,
    censor_sigma = censor_sigma,
    seed = trn_seed
  )

  tst_obj <- simulate_surv_data7(
    n = n_target_test,
    dataset = "target",
    truth_par = truth_target,
    time_par = time_target,
    cov_shift = cov_shift,
    target_censor_rate = censor_target_test,
    censor_sigma = censor_sigma,
    seed = tst_seed
  )

  list(
    source_data = src_obj$data,
    target_train = trn_obj$data,
    target_test = tst_obj$data,

    source_info = src_obj,
    target_train_info = trn_obj,
    target_test_info = tst_obj,

    covariates = paste0("X", 1:7),

    truth_source = truth_source,
    truth_target = truth_target,

    time_source = time_source,
    time_target = time_target,

    cov_shift = cov_shift,
    censor_source = censor_source,
    censor_target_train = censor_target_train,
    censor_target_test = censor_target_test,
    censor_sigma = censor_sigma
  )
}


################################################################################
# 6. Transfer truth constructor
################################################################################

make_transfer_truth7 <- function(source_truth = default_truth7(TRUE),
                                 delta_linear = c(0.45, -0.40, 0.25, 0, 0, 0, 0),
                                 delta_nonlinear = c(I_X4_2 = 0,
                                                     I_X1X3 = 0,
                                                     I_sinX3 = 0),
                                 shift_scale = 1) {
  src <- source_truth

  beta_linear_target <- src$beta_linear + shift_scale * delta_linear

  if (src$use_nonlinear_in_truth) {
    src_nl <- src$beta_nonlinear

    nm <- union(names(src_nl), names(delta_nonlinear))

    src_pad <- rep(0, length(nm))
    names(src_pad) <- nm

    del_pad <- rep(0, length(nm))
    names(del_pad) <- nm

    src_pad[names(src_nl)] <- src_nl
    del_pad[names(delta_nonlinear)] <- delta_nonlinear

    beta_nonlinear_target <- src_pad + shift_scale * del_pad
  } else {
    beta_nonlinear_target <- numeric(0)
  }

  default_truth7(
    use_nonlinear_in_truth = src$use_nonlinear_in_truth,
    intercept = unname(src$beta_true["Intercept"]),
    beta_linear = beta_linear_target,
    beta_nonlinear = beta_nonlinear_target
  )
}


################################################################################
# 7. Baseline estimators and prediction helpers
################################################################################

estimate_beta_target_only <- function(target_train,
                                      t,
                                      covariates = paste0("X", 1:7),
                                      eps_G = 1e-3,
                                      w_cap = 20) {
  estimate_betaS(
    data = target_train,
    t = t,
    Y_col = "Y",
    delta_col = "delta",
    covariates = covariates,
    eps_G = eps_G,
    w_cap = w_cap
  )
}


estimate_beta_pooled <- function(source_data,
                                 target_train,
                                 t,
                                 covariates = paste0("X", 1:7),
                                 eps_G = 1e-3,
                                 w_cap = 20) {
  pooled <- rbind(source_data, target_train)

  estimate_betaS(
    data = pooled,
    t = t,
    Y_col = "Y",
    delta_col = "delta",
    covariates = covariates,
    eps_G = eps_G,
    w_cap = w_cap
  )
}


estimate_beta_target_only_cox <- function(target_train,
                                          covariates = paste0("X", 1:7)) {
  fml <- as.formula(
    paste("Surv(Y, delta) ~", paste(covariates, collapse = " + "))
  )

  fit <- survival::coxph(fml, data = target_train, ties = "breslow")
  beta <- stats::coef(fit)
  beta_out <- c(Intercept = 0, beta)

  list(
    beta = beta_out,
    fit = fit
  )
}


predict_logistic_risk_t <- function(beta,
                                    newdata,
                                    covariates = paste0("X", 1:7)) {
  lp <- linear_lp_from_beta(
    beta = beta,
    newdata = newdata,
    covariates = covariates
  )

  list(
    lp = lp,
    p_event_t = stats::plogis(lp)
  )
}


predict_cox_risk_t <- function(cox_fit, newdata, t) {
  lp <- as.vector(stats::predict(cox_fit, newdata = newdata, type = "lp"))

  bh <- survival::basehaz(cox_fit, centered = FALSE)
  idx <- findInterval(t, bh$time)
  H0_t <- if (idx == 0) 0 else bh$hazard[idx]

  p_event_t <- 1 - exp(-H0_t * exp(lp))

  list(
    lp = lp,
    p_event_t = p_event_t
  )
}


predict_oracle_risk_t <- function(sim_data, t) {
  eta_true <- sim_data$target_test_info$eta_true
  time_par <- sim_data$time_target

  if (time_par$time_model == "weibull") {
    p_event_t <- 1 - exp(
      - exp(eta_true) * (t / time_par$weib_scale)^time_par$weib_shape
    )
    risk_score <- eta_true
  } else if (time_par$time_model == "aft_logistic") {
    p_event_t <- stats::plogis(
      (log(t / time_par$t0) + eta_true) / time_par$aft_logistic_scale
    )
    risk_score <- eta_true
  } else {
    stop("Unsupported time model in predict_oracle_risk_t().")
  }

  list(
    lp = risk_score,
    p_event_t = p_event_t
  )
}


################################################################################
# 8. Evaluation
################################################################################

evaluate_surv_metrics_at_t <- function(target_test,
                                       t,
                                       risk_score,
                                       p_event_t,
                                       eps_G = 1e-3) {
  D_t <- as.integer(target_test$Y <= t & target_test$delta == 1)
  risk_score <- as.numeric(risk_score)
  p_event_t <- as.numeric(p_event_t)

  n_case <- sum(D_t == 1, na.rm = TRUE)
  n_ctrl <- sum(D_t == 0, na.rm = TRUE)

  w_ipcw <- estimate_weights(
    data = target_test,
    t = t,
    Y_col = "Y",
    delta_col = "delta",
    eps_G = eps_G,
    w_cap = Inf
  )

  brier_t <- mean(w_ipcw * (D_t - p_event_t)^2, na.rm = TRUE)
  if (!is.finite(brier_t)) brier_t <- NA_real_

  p_clip <- pmin(pmax(p_event_t, 1e-6), 1 - 1e-6)
  log_loss_t <- mean(
    w_ipcw * (-(D_t * log(p_clip) + (1 - D_t) * log(1 - p_clip))),
    na.rm = TRUE
  )
  if (!is.finite(log_loss_t)) log_loss_t <- NA_real_

  risk_sd <- if (sum(is.finite(risk_score)) >= 2) {
    stats::sd(risk_score[is.finite(risk_score)])
  } else {
    NA_real_
  }

  if (n_case < 5 || n_ctrl < 5 || !is.finite(risk_sd) || risk_sd < 1e-10 ||
      any(!is.finite(risk_score))) {
    auc_t <- NA_real_
  } else {
    auc_t <- tryCatch({
      auc_obj <- timeROC::timeROC(
        T = target_test$Y,
        delta = target_test$delta,
        marker = risk_score,
        cause = 1,
        weighting = "marginal",
        times = t,
        iid = FALSE
      )
      as.numeric(tail(auc_obj$AUC, 1))
    }, error = function(e) NA_real_)
  }

  harrell_c <- tryCatch({
    harrell_fit <- survival::concordance(
      survival::Surv(Y, delta) ~ risk_score,
      data = target_test,
      reverse = TRUE,
      timewt = "n"
    )
    as.numeric(harrell_fit$concordance)
  }, error = function(e) NA_real_)

  uno_c <- tryCatch({
    uno_fit <- survival::concordance(
      survival::Surv(Y, delta) ~ risk_score,
      data = target_test,
      reverse = TRUE,
      timewt = "n/G2"
    )
    as.numeric(uno_fit$concordance)
  }, error = function(e) NA_real_)

  list(
    auc_t = auc_t,
    brier_t = brier_t,
    log_loss_t = log_loss_t,
    uno_c = uno_c,
    harrell_c = harrell_c,
    n_case = n_case,
    n_ctrl = n_ctrl
  )
}


beta_l2_error <- function(beta_hat, beta_true) {
  nm <- union(names(beta_hat), names(beta_true))

  b1 <- rep(0, length(nm))
  names(b1) <- nm

  b2 <- rep(0, length(nm))
  names(b2) <- nm

  b1[names(beta_hat)] <- beta_hat
  b2[names(beta_true)] <- beta_true

  sqrt(sum((b1 - b2)^2))
}


################################################################################
# 9. Method runner
################################################################################

run_one_method7 <- function(method,
                            sim_data,
                            t,
                            source_covariates = paste0("X", 1:7),
                            target_covariates = paste0("X", 1:7),
                            common_cov = paste0("X", 1:7),
                            target_only_cov = character(0),
                            lambda_w = 1,
                            eps_G = 1e-3,
                            w_cap = 20,
                            max_folds = 5,
                            beta_target_ref = NULL,
                            gamma_source = 1.4,
                            gamma_target = 1.4,
                            select_source = FALSE,
                            select_target = FALSE,
                            use_interaction = FALSE,
                            include_binary = TRUE) {
  source_data  <- sim_data$source_data
  target_train <- sim_data$target_train
  target_test  <- sim_data$target_test

  # Backward-compatible alias.
  if (method == "trans_logistic") {
    method <- "trans_linear_linear"
  }

  if (is.null(beta_target_ref)) {
    beta_target_ref <- c(
      Intercept = 0,
      setNames(rep(0, length(target_covariates)), target_covariates)
    )
  }

  if (method %in% c("source_only", "source_only_logistic")) {
    beta_hat <- estimate_betaS(
      data = source_data,
      t = t,
      covariates = source_covariates,
      eps_G = eps_G,
      w_cap = w_cap
    )

    pred <- predict_logistic_risk_t(beta_hat, target_test, target_covariates)
    met <- evaluate_surv_metrics_at_t(target_test, t, pred$lp, pred$p_event_t, eps_G)

    return(list(
      method = method,
      beta = beta_hat,
      eta = NULL,
      source_beta = beta_hat,
      beta_l2 = beta_l2_error(beta_hat, beta_target_ref),
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  if (method %in% c("target_only", "target_only_logistic")) {
    beta_hat <- estimate_beta_target_only(
      target_train = target_train,
      t = t,
      covariates = target_covariates,
      eps_G = eps_G,
      w_cap = w_cap
    )

    pred <- predict_logistic_risk_t(beta_hat, target_test, target_covariates)
    met <- evaluate_surv_metrics_at_t(target_test, t, pred$lp, pred$p_event_t, eps_G)

    return(list(
      method = method,
      beta = beta_hat,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = beta_l2_error(beta_hat, beta_target_ref),
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  if (method %in% c("source_only_gam", "target_only_gam", "pooled_gam")) {
    gam_train <- if (method == "source_only_gam") {
      source_data
    } else if (method == "target_only_gam") {
      target_train
    } else {
      rbind(source_data, target_train)
    }
    gam_gamma <- if (method == "source_only_gam") gamma_source else gamma_target
    gam_select <- if (method == "source_only_gam") select_source else select_target
    fit_obj <- fit_target_only_gam_model(
      target_train = gam_train,
      t = t,
      cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
      use_interaction = use_interaction,
      include_binary = include_binary,
      eps_G = eps_G,
      w_cap = w_cap,
      gamma = gam_gamma,
      select = gam_select
    )

    pred <- predict_target_only_gam_risk_t(
      object = fit_obj,
      newdata = target_test
    )

    met <- evaluate_surv_metrics_at_t(
      target_test = target_test,
      t = t,
      risk_score = pred$lp,
      p_event_t = pred$p_event_t,
      eps_G = eps_G
    )

    return(list(
      method = method,
      fit_obj = fit_obj,
      beta = NULL,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = NA_real_,
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  if (method %in% c("pooled", "pooled_logistic")) {
    beta_hat <- estimate_beta_pooled(
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = target_covariates,
      eps_G = eps_G,
      w_cap = w_cap
    )

    pred <- predict_logistic_risk_t(beta_hat, target_test, target_covariates)
    met <- evaluate_surv_metrics_at_t(target_test, t, pred$lp, pred$p_event_t, eps_G)

    return(list(
      method = method,
      beta = beta_hat,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = beta_l2_error(beta_hat, beta_target_ref),
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  if (method %in% c("source_only_rsf",
                    "target_only_rsf",
                    "pooled_rsf",
                    "source_only_ipcw_xgb",
                    "target_only_ipcw_xgb",
                    "pooled_ipcw_xgb")) {
    train_data <- if (startsWith(method, "source_only")) {
      source_data
    } else if (startsWith(method, "target_only")) {
      target_train
    } else {
      rbind(source_data, target_train)
    }
    covariates_use <- if (startsWith(method, "source_only")) {
      source_covariates
    } else {
      target_covariates
    }

    if (grepl("_rsf$", method)) {
      fit_obj <- fit_rsf_survival_model(
        train_data = train_data,
        covariates = covariates_use,
        seed = 202605 + nrow(train_data)
      )
      pred <- predict_rsf_risk_t(fit_obj, target_test, t)
    } else {
      fit_obj <- fit_ipcw_xgb_model(
        train_data = train_data,
        t = t,
        covariates = covariates_use,
        eps_G = eps_G,
        w_cap = w_cap,
        seed = 202605 + nrow(train_data)
      )
      pred <- predict_ipcw_xgb_risk_t(fit_obj, target_test)
    }

    met <- evaluate_surv_metrics_at_t(
      target_test = target_test,
      t = t,
      risk_score = pred$lp,
      p_event_t = pred$p_event_t,
      eps_G = eps_G
    )

    return(list(
      method = method,
      fit_obj = fit_obj,
      beta = NULL,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = NA_real_,
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  if (method %in% c("transcox_bic", "transcox_fixed")) {
    fit_obj <- fit_transcox_model(
      source_data = source_data,
      target_train = target_train,
      source_covariates = source_covariates,
      target_covariates = target_covariates,
      select_bic = method == "transcox_bic"
    )
    pred <- predict_transcox_risk_t(fit_obj, target_test, t)
    met <- evaluate_surv_metrics_at_t(target_test, t, pred$lp, pred$p_event_t, eps_G)

    return(list(
      method = method,
      fit_obj = fit_obj,
      beta = NULL,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = NA_real_,
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  if (method %in% c("trans_linear_linear",
                    "trans_linear_gam",
                    "trans_linear_gam_calibrated",
                    "trans_gam_linear",
                    "trans_gam_linear_calibrated",
                    "trans_gam_auto_calibrated",
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
                    "trans_gam_gam_auto_cal",
                    "trans_gam_gam",
                    "trans_gam_gam_tp",
                    "trans_gam_gam_source_tp",
                    "trans_gam_gam_all_tp",
                    "trans_gam_gam_select",
                    "trans_gam_gam_calibrated",
                    "trans_gam_gam_no_tilt",
                    "trans_gam_gam_tp_no_tilt",
                    "trans_gam_gam_source_tp_no_tilt",
                    "trans_gam_gam_all_tp_no_tilt",
                    "trans_gam_gam_select_no_tilt",
                    "trans_gam_gam_calibrated_no_tilt",
                    "trans_gam_gam_partial_cal_no_tilt",
                    "trans_gam_gam_penalized_cal_no_tilt",
                    "trans_gam_gam_partial_cal",
                    "trans_gam_gam_penalized_cal",
                    "source_gam_recal_intercept",
                    "source_gam_recal_slope")) {
    fit_obj <- fit_transfer_method(
      method = method,
      source_data = source_data,
      target_train = target_train,
      t = t,
      covariates = target_covariates,
      cont_covariates = c("X3", "X4", "X5", "X6", "X7"),
      use_interaction = use_interaction,
      include_binary = include_binary,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      alpha_eta = 1,
      lambda_choice = "lambda.1se",
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target
    )

    pred <- predict_transfer_risk_t(
      object = fit_obj,
      newdata = target_test
    )

    met <- evaluate_surv_metrics_at_t(
      target_test = target_test,
      t = t,
      risk_score = pred$lp,
      p_event_t = pred$p_event_t,
      eps_G = eps_G
    )

    beta_l2_value <- NA_real_
    if (method == "trans_linear_linear") {
      beta_l2_value <- beta_l2_error(fit_obj$beta, beta_target_ref)
    }

    source_alpha_value <- NA_real_
    if (exists("extract_source_alpha", mode = "function")) {
      source_alpha_value <- extract_source_alpha(fit_obj)
    }

    return(list(
      method = method,
      fit_obj = fit_obj,
      beta = if (method == "trans_linear_linear") fit_obj$beta else NULL,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = beta_l2_value,
      source_alpha = source_alpha_value,
      partial_weight = if (is.null(fit_obj$partial_weight)) NA_real_ else fit_obj$partial_weight,
      alpha_penalty_estimated = if (is.null(fit_obj$alpha_penalty_estimated)) NA_real_ else fit_obj$alpha_penalty_estimated,
      selected_shift = if (is.null(fit_obj$selected_shift)) NA_character_ else fit_obj$selected_shift,
      selected_target = if (is.null(fit_obj$selected_target)) NA_character_ else fit_obj$selected_target,
      selected_calibration = if (is.null(fit_obj$selected_calibration)) NA_character_ else fit_obj$selected_calibration,
      auto_score_source = if (is.null(fit_obj$auto_scores)) NA_real_ else fit_obj$auto_scores$score_mean[fit_obj$auto_scores$target_model == "source"],
      auto_score_linear = if (is.null(fit_obj$auto_scores)) NA_real_ else fit_obj$auto_scores$score_mean[fit_obj$auto_scores$target_model == "linear"],
      auto_score_gam = if (is.null(fit_obj$auto_scores)) NA_real_ else fit_obj$auto_scores$score_mean[fit_obj$auto_scores$target_model == "gam"],
      auto_correction_gain = if (is.null(fit_obj$auto_correction_gain)) NA_real_ else fit_obj$auto_correction_gain,
      auto_gam_gain = if (is.null(fit_obj$auto_gam_gain)) NA_real_ else fit_obj$auto_gam_gain,
      auto_gam_gain_se = if (is.null(fit_obj$auto_gam_gain_se)) NA_real_ else fit_obj$auto_gam_gain_se,
      auto_score_no_cal = if (is.null(fit_obj$auto_score_no_cal)) NA_real_ else fit_obj$auto_score_no_cal,
      auto_score_cal = if (is.null(fit_obj$auto_score_cal)) NA_real_ else fit_obj$auto_score_cal,
      auto_cal_loss_diff = if (is.null(fit_obj$auto_cal_loss_diff)) NA_real_ else fit_obj$auto_cal_loss_diff,
      auto_cal_loss_diff_se = if (is.null(fit_obj$auto_cal_loss_diff_se)) NA_real_ else fit_obj$auto_cal_loss_diff_se,
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  if (method %in% c("source_only_cox", "target_only_cox", "pooled_cox")) {
    cox_train <- if (method == "source_only_cox") {
      source_data
    } else if (method == "target_only_cox") {
      target_train
    } else {
      rbind(source_data, target_train)
    }
    cox_covariates <- if (method == "source_only_cox") {
      source_covariates
    } else {
      target_covariates
    }
    cox_obj <- estimate_beta_target_only_cox(
      target_train = cox_train,
      covariates = cox_covariates
    )

    pred <- predict_cox_risk_t(cox_obj$fit, target_test, t)
    met <- evaluate_surv_metrics_at_t(target_test, t, pred$lp, pred$p_event_t, eps_G)

    return(list(
      method = method,
      beta = cox_obj$beta,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = beta_l2_error(cox_obj$beta, beta_target_ref),
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c,
      extra = cox_obj
    ))
  }

  if (method == "oracle") {
    pred <- predict_oracle_risk_t(sim_data = sim_data, t = t)
    met <- evaluate_surv_metrics_at_t(target_test, t, pred$lp, pred$p_event_t, eps_G)

    return(list(
      method = method,
      beta = beta_target_ref,
      eta = NULL,
      source_beta = NULL,
      beta_l2 = 0,
      auc_t = met$auc_t,
      brier_t = met$brier_t,
      log_loss_t = met$log_loss_t,
      uno_c = met$uno_c,
      harrell_c = met$harrell_c
    ))
  }

  stop("Unknown method: ", method)
}


run_all_methods_one_setting7 <- function(sim_data,
                                         t,
                                         methods = c("source_only_logistic",
                                                     "target_only_logistic",
                                                     "pooled_logistic",
                                                     "source_only_gam",
                                                     "target_only_gam",
                                                     "pooled_gam",
                                                     "source_only_rsf",
                                                     "target_only_rsf",
                                                     "pooled_rsf",
                                                     "source_only_ipcw_xgb",
                                                     "target_only_ipcw_xgb",
                                                     "pooled_ipcw_xgb",
                                                     "source_only_cox",
                                                     "target_only_cox",
                                                     "pooled_cox",
                                                     "source_gam_recal_intercept",
                                                     "source_gam_recal_slope",
                                                     "trans_linear_linear",
                                                     "trans_linear_gam",
                                                     "trans_linear_gam_calibrated",
                                                     "trans_gam_linear",
                                                     "trans_gam_linear_calibrated",
                                                     "trans_gam_auto_calibrated",
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
                                                     "trans_gam_gam_auto_cal",
                                                     "trans_gam_gam_no_tilt",
                                                     "trans_gam_gam_tp_no_tilt",
                                                     "trans_gam_gam_source_tp_no_tilt",
                                                     "trans_gam_gam_all_tp_no_tilt",
                                                     "trans_gam_gam_select_no_tilt",
                                                     "trans_gam_gam_calibrated_no_tilt",
                                                     "trans_gam_gam_partial_cal_no_tilt",
                                                     "trans_gam_gam_penalized_cal_no_tilt",
                                                     "trans_gam_gam",
                                                     "trans_gam_gam_tp",
                                                     "trans_gam_gam_source_tp",
                                                     "trans_gam_gam_all_tp",
                                                     "trans_gam_gam_select",
                                                     "trans_gam_gam_calibrated",
                                                     "trans_gam_gam_partial_cal",
                                                     "trans_gam_gam_penalized_cal",
                                                     "oracle"),
                                         source_covariates = paste0("X", 1:7),
                                         target_covariates = paste0("X", 1:7),
                                         common_cov = paste0("X", 1:7),
                                         target_only_cov = character(0),
                                         lambda_w = 1,
                                         eps_G = 1e-3,
                                         w_cap = 20,
                                         max_folds = 5,
                                         beta_target_ref = NULL,
                                         gamma_source = 1.4,
                                         gamma_target = 1.4,
                                         select_source = FALSE,
                                         select_target = FALSE,
                                         use_interaction = FALSE,
                                         include_binary = TRUE) {
  out <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    out[[i]] <- run_one_method7(
      method = methods[i],
      sim_data = sim_data,
      t = t,
      source_covariates = source_covariates,
      target_covariates = target_covariates,
      common_cov = common_cov,
      target_only_cov = target_only_cov,
      lambda_w = lambda_w,
      eps_G = eps_G,
      w_cap = w_cap,
      max_folds = max_folds,
      beta_target_ref = beta_target_ref,
      gamma_source = gamma_source,
      gamma_target = gamma_target,
      select_source = select_source,
      select_target = select_target,
      use_interaction = use_interaction,
      include_binary = include_binary
    )
  }

  out
}


################################################################################
# 10. Setting generator
################################################################################

make_recommended_setting7 <- function(n_source = 1000,
                                      n_target_train = 200,
                                      n_target_test = 1000,
                                      shift_level = c("moderate", "strong"),
                                      use_nonlinear_in_truth = TRUE,
                                      time_model = c("weibull", "aft_logistic"),
                                      shift_scale = 1,
                                      censor_source = 0.35,
                                      censor_target_train = 0.35,
                                      censor_target_test = 0.35,
                                      censor_sigma = 0.6) {
  shift_level <- match.arg(shift_level)
  time_model <- match.arg(time_model)

  cov_shift <- default_cov_shift7(shift_level)

  truth_source <- default_truth7(
    use_nonlinear_in_truth = use_nonlinear_in_truth
  )

  truth_target <- make_transfer_truth7(
    source_truth = truth_source,
    delta_linear = c(0.45, -0.40, 0.25, 0, 0, 0, 0),
    delta_nonlinear = c(I_X4_2 = 0, I_X1X3 = 0, I_sinX3 = 0),
    shift_scale = shift_scale
  )

  time_par <- default_time7(time_model)

  list(
    n_source = n_source,
    n_target_train = n_target_train,
    n_target_test = n_target_test,

    cov_shift = cov_shift,

    truth_source = truth_source,
    truth_target = truth_target,

    time_source = time_par,
    time_target = time_par,

    censor_source = censor_source,
    censor_target_train = censor_target_train,
    censor_target_test = censor_target_test,
    censor_sigma = censor_sigma
  )
}


generate_from_setting7 <- function(setting, seed = NULL) {
  generate_train_test_surv_data7(
    n_source = setting$n_source,
    n_target_train = setting$n_target_train,
    n_target_test = setting$n_target_test,

    truth_source = setting$truth_source,
    truth_target = setting$truth_target,

    time_source = setting$time_source,
    time_target = setting$time_target,

    cov_shift = setting$cov_shift,

    censor_source = setting$censor_source,
    censor_target_train = setting$censor_target_train,
    censor_target_test = setting$censor_target_test,
    censor_sigma = ifelse(is.null(setting$censor_sigma), 0.6, setting$censor_sigma),

    seed = seed
  )
}


################################################################################
# 11. AFT-logistic nonlinear 2 x 2 setting constructor
#
# Four settings:
#   - AFT_NL_few_terms_small_shift
#   - AFT_NL_few_terms_large_shift
#   - AFT_NL_many_terms_small_shift
#   - AFT_NL_many_terms_large_shift
#
# Factor 1: nonlinear complexity
#   few  = only X4^2 active in the source nonlinear truth
#   many = X4^2 and sin(X3) active, with no interaction term
#
# Factor 2: source-target discrepancy
#   small = source score plus a small linear target residual
#   large = source score plus a larger linear target residual
################################################################################

make_aft_nonlinear_setting7 <- function(setting_name,
                                        nonlinear_type = c("few", "many"),
                                        shift_size = c("small", "large"),
                                        n_source = 1000,
                                        n_target_train = 150,
                                        n_target_test = 1000,
                                        shift_level = "moderate",
                                        censor_source = 0.35,
                                        censor_target_train = 0.30,
                                        censor_target_test = 0.30,
                                        censor_sigma = 0.6,
                                        t0 = 10,
                                        aft_logistic_scale = 0.7) {
  nonlinear_type <- match.arg(nonlinear_type)
  shift_size <- match.arg(shift_size)

  setting <- make_recommended_setting7(
    n_source = n_source,
    n_target_train = n_target_train,
    n_target_test = n_target_test,
    shift_level = shift_level,
    use_nonlinear_in_truth = TRUE,
    time_model = "aft_logistic",
    shift_scale = 1.0,
    censor_source = censor_source,
    censor_target_train = censor_target_train,
    censor_target_test = censor_target_test,
    censor_sigma = censor_sigma
  )

  # Explicitly set AFT-logistic time parameters.
  setting$time_source <- default_time7(
    time_model = "aft_logistic",
    t0 = t0,
    aft_logistic_scale = aft_logistic_scale
  )
  setting$time_target <- default_time7(
    time_model = "aft_logistic",
    t0 = t0,
    aft_logistic_scale = aft_logistic_scale
  )

  # Source truth:
  # Match the Weibull nonlinear extended design. The true DGP intentionally
  # excludes interactions, matching the main-effect GAM/linear comparison.
  beta_linear_source <- c(
    0.95,   # X1
    -0.80,  # X2
    0.55,   # X3
    -0.45,  # X4
    0.35,   # X5
    0.00,   # X6
    0.00    # X7
  )

  beta_nl_source <- if (nonlinear_type == "few") {
    c(I_X4_2 = 0.70, I_X1X3 = 0.00, I_sinX3 = 0.00)
  } else {
    c(I_X4_2 = 0.70, I_X1X3 = 0.00, I_sinX3 = 0.60)
  }

  setting$truth_source <- default_truth7(
    use_nonlinear_in_truth = TRUE,
    intercept = 0,
    beta_linear = beta_linear_source,
    beta_nonlinear = beta_nl_source
  )

  # Target truth:
  #   eta_target = alpha_source * eta_source + linear_residual(X).
  # This creates the same qualitative regimes as the Weibull extended design:
  # small shifts are nearly source-score sufficient, while large shifts require
  # a target correction and often benefit from the flexible GAM correction.
  if (shift_size == "small") {
    alpha_source <- 0.80
    linear_residual <- c(0.10, -0.08, 0.08, -0.06, 0.06, 0.20, -0.20)
    shift_design <- "source_score_plus_small_linear_residual"
  } else {
    alpha_source <- 0.55
    linear_residual <- c(0.18, -0.14, 0.12, -0.10, 0.10, 0.65, -0.60)
    shift_design <- "source_score_plus_large_linear_residual"
  }

  beta_linear_target <- alpha_source * beta_linear_source + linear_residual
  beta_nl_target <- alpha_source * beta_nl_source

  delta_linear <- beta_linear_target - beta_linear_source

  nl_names <- union(names(beta_nl_source), names(beta_nl_target))
  src_nl_pad <- rep(0, length(nl_names))
  names(src_nl_pad) <- nl_names
  src_nl_pad[names(beta_nl_source)] <- beta_nl_source

  tgt_nl_pad <- rep(0, length(nl_names))
  names(tgt_nl_pad) <- nl_names
  tgt_nl_pad[names(beta_nl_target)] <- beta_nl_target

  delta_nonlinear <- tgt_nl_pad - src_nl_pad

  setting$truth_target <- make_transfer_truth7(
    source_truth = setting$truth_source,
    delta_linear = delta_linear,
    delta_nonlinear = delta_nonlinear,
    shift_scale = 1.0
  )

  setting$setting_name <- setting_name
  setting$nonlinear_type <- nonlinear_type
  setting$shift_size <- shift_size
  setting$delta_linear <- delta_linear
  setting$delta_nonlinear <- delta_nonlinear
  setting$beta_linear_source <- beta_linear_source
  setting$beta_linear_target <- beta_linear_target
  setting$beta_nonlinear_source <- beta_nl_source
  setting$beta_nonlinear_target <- beta_nl_target
  setting$alpha_source <- alpha_source
  setting$linear_residual <- linear_residual
  setting$shift_design <- shift_design
  setting$t0 <- t0
  setting$aft_logistic_scale <- aft_logistic_scale

  setting
}


make_setting_AFT_NL_FewSmall <- function(...) {
  make_aft_nonlinear_setting7(
    setting_name = "AFT_NL_few_terms_small_shift",
    nonlinear_type = "few",
    shift_size = "small",
    ...
  )
}


make_setting_AFT_NL_FewLarge <- function(...) {
  make_aft_nonlinear_setting7(
    setting_name = "AFT_NL_few_terms_large_shift",
    nonlinear_type = "few",
    shift_size = "large",
    ...
  )
}


make_setting_AFT_NL_ManySmall <- function(...) {
  make_aft_nonlinear_setting7(
    setting_name = "AFT_NL_many_terms_small_shift",
    nonlinear_type = "many",
    shift_size = "small",
    ...
  )
}


make_setting_AFT_NL_ManyLarge <- function(...) {
  make_aft_nonlinear_setting7(
    setting_name = "AFT_NL_many_terms_large_shift",
    nonlinear_type = "many",
    shift_size = "large",
    ...
  )
}
