cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1])) else normalizePath("scripts/main_simulation_server.R", mustWork = FALSE)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
r_dir <- file.path(repo_root, "R")
options(transgam.r_dir = r_dir)

source(file.path(r_dir, "server_common.R"))
source(file.path(r_dir, "transfer_gam_functions.R"))

suppressPackageStartupMessages({
  library(MASS)
  library(survival)
  library(glmnet)
  library(parallel)
})

parse_cli <- function(defaults) {
  raw <- commandArgs(trailingOnly = TRUE)
  if (any(raw %in% c("-h", "--help"))) {
    cat("Usage: Rscript scripts/main_simulation_server.R [options]\n\n")
    cat("Common options:\n")
    cat("  --model weibull|aftlogistic|both\n")
    cat("  --methods stat,transfer,gam_gam_basis,oracle\n")
    cat("  --settings all|aligned_linear|aligned_nonlinear|misaligned_linear|misaligned_nonlinear\n")
    cat("  --n-rep 100 --rep-start 1 --rep-end 10 --n-s 2000 --n-t 300 --n-test 1000\n")
    cat("  --tq \"{0.25,0.50,0.75}\" --out-root results/run --cores 8\n")
    quit(save = "no", status = 0)
  }

  out <- defaults
  bool_flags <- c("no_parallel", "save_full_objects", "select_source", "select_target", "use_interaction")
  numeric_fields <- c(
    "n_rep", "rep_start", "rep_end", "n_s", "n_t", "n_test", "cores",
    "max_folds", "rsf_trees", "rsf_min_node_size", "rsf_threads",
    "xgb_threads", "gamma_source", "gamma_target", "alpha_penalty",
    "transcox_lambda1", "transcox_lambda2", "transcox_learning_rate",
    "transcox_nsteps", "base_seed"
  )

  i <- 1L
  while (i <= length(raw)) {
    token <- raw[i]
    if (!grepl("^--", token)) stop("Unexpected command-line token: ", token)
    key <- sub("^--", "", token)
    value <- NULL
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      value <- paste(parts[-1], collapse = "=")
    }
    key <- gsub("-", "_", key, fixed = TRUE)

    if (key %in% bool_flags) {
      out[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (is.null(value)) {
      if (i == length(raw)) stop("Missing value for --", gsub("_", "-", key))
      value <- raw[i + 1L]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
    if (!key %in% names(out)) stop("Unknown option --", gsub("_", "-", key))
    out[[key]] <- if (key %in% numeric_fields) as.numeric(value) else value
  }

  for (nm in c("n_rep", "rep_start", "rep_end", "n_s", "n_t", "n_test", "cores",
               "max_folds", "rsf_trees", "rsf_min_node_size", "rsf_threads",
               "xgb_threads", "transcox_nsteps", "base_seed")) {
    out[[nm]] <- as.integer(out[[nm]])
  }
  out
}

args <- parse_cli(list(
  model = "weibull",
  methods = "all",
  settings = "all",
  n_rep = 50L,
  rep_start = 1L,
  rep_end = -1L,
  n_s = 1000L,
  n_t = 150L,
  n_test = 1000L,
  tq = "{0.50}",
  out_root = file.path(repo_root, "results"),
  out_dir = "",
  cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")),
  no_parallel = FALSE,
  save_full_objects = FALSE,
  max_folds = 5L,
  rsf_trees = 300L,
  rsf_min_node_size = 20L,
  rsf_threads = 1L,
  xgb_threads = 1L,
  gamma_source = 1.4,
  gamma_target = 1.4,
  alpha_penalty = 10,
  transcox_lambda1 = 0.1,
  transcox_lambda2 = 0.1,
  transcox_lambda_grid = "{0.1,0.5,1,2}",
  transcox_learning_rate = 0.004,
  transcox_nsteps = 100L,
  transcox_backend = "r",
  select_source = FALSE,
  select_target = FALSE,
  use_interaction = FALSE,
  base_seed = 20260509L
))

if (args$rep_end < 0L) {
  args$rep_end <- args$n_rep
}
if (args$rep_start < 1L) stop("--rep-start must be >= 1.")
if (args$rep_end < args$rep_start) stop("--rep-end must be >= --rep-start.")

options(
  translogistic.rsf_trees = args$rsf_trees,
  translogistic.rsf_min_node_size = args$rsf_min_node_size,
  translogistic.rsf_threads = max(1L, args$rsf_threads),
  translogistic.xgb_threads = max(1L, args$xgb_threads),
  translogistic.alpha_penalty = args$alpha_penalty
)

methods_to_run <- parse_methods(args$methods, include_oracle = TRUE)
tq_list <- parse_tq_list(args$tq)
transcox_lambda_grid <- as.numeric(parse_token_list(args$transcox_lambda_grid))
if (length(transcox_lambda_grid) == 0 || any(!is.finite(transcox_lambda_grid)) || any(transcox_lambda_grid < 0)) {
  stop("--transcox-lambda-grid must contain nonnegative numeric values.")
}

options(
  translogistic.transcox_lambda1 = args$transcox_lambda1,
  translogistic.transcox_lambda2 = args$transcox_lambda2,
  translogistic.transcox_lambda1_vec = transcox_lambda_grid,
  translogistic.transcox_lambda2_vec = transcox_lambda_grid,
  translogistic.transcox_learning_rate = args$transcox_learning_rate,
  translogistic.transcox_nsteps = args$transcox_nsteps,
  translogistic.transcox_backend = args$transcox_backend
)

select_setting_keys <- function(x) {
  tokens <- parse_token_list(x)
  if (length(tokens) == 0 || any(tokens == "all")) {
    return(c("aligned_linear", "aligned_nonlinear", "misaligned_linear", "misaligned_nonlinear"))
  }
  aliases <- c(
    fewsmall = "few_small",
    fewlarge = "few_large",
    manysmall = "many_small",
    manylarge = "many_large",
    alignedlinear = "aligned_linear",
    alignednonlinear = "aligned_nonlinear",
    misalignedlinear = "misaligned_linear",
    misalignednonlinear = "misaligned_nonlinear",
    slope_linear = "misaligned_linear",
    slope_nonlinear = "misaligned_nonlinear",
    legacy_all = "legacy_all"
  )
  tokens <- ifelse(tokens %in% names(aliases), aliases[tokens], tokens)
  if (any(tokens == "legacy_all")) {
    tokens <- c(setdiff(tokens, "legacy_all"), "few_small", "few_large", "many_small", "many_large")
  }
  valid <- c(
    "aligned_linear", "aligned_nonlinear",
    "misaligned_linear", "misaligned_nonlinear",
    "few_small", "few_large", "many_small", "many_large"
  )
  bad <- setdiff(tokens, valid)
  if (length(bad) > 0) stop("Unknown setting key(s): ", paste(bad, collapse = ", "))
  unique(tokens)
}

setting_keys <- select_setting_keys(args$settings)

extract_linear_beta_ref <- function(truth_obj) {
  c(
    Intercept = unname(truth_obj$beta_true["Intercept"]),
    setNames(truth_obj$beta_linear, paste0("X", 1:7))
  )
}

make_adaptive_dgp_setting7 <- function(setting_name,
                                       residual_type = c("linear", "nonlinear"),
                                       source_alignment = c("aligned", "misaligned"),
                                       n_source,
                                       n_target_train,
                                       n_target_test,
                                       time_model = c("weibull", "aft_logistic"),
                                       censor_source = 0.35,
                                       censor_target_train = 0.30,
                                       censor_target_test = 0.30) {
  residual_type <- match.arg(residual_type)
  source_alignment <- match.arg(source_alignment)
  time_model <- match.arg(time_model)

  setting <- make_recommended_setting7(
    n_source = n_source,
    n_target_train = n_target_train,
    n_target_test = n_target_test,
    shift_level = "moderate",
    use_nonlinear_in_truth = TRUE,
    time_model = time_model,
    shift_scale = 1.0,
    censor_source = censor_source,
    censor_target_train = censor_target_train,
    censor_target_test = censor_target_test
  )

  beta_linear_source <- c(0.95, -0.80, 0.55, -0.45, 0.35, 0.00, 0.00)
  beta_nl_source <- c(I_X4_2 = 0.70, I_X1X3 = 0.00, I_sinX3 = 0.60)

  if (time_model == "aft_logistic") {
    # AFT-logistic needs a slightly stronger smooth target residual than the
    # Weibull DGP; otherwise linear Cox scores can explain too much of the
    # target ranking. Keep the same 2 x 2 setting names and CLI interface.
    beta_linear_source <- c(0.85, -0.70, 0.45, -0.35, 0.25, 0.00, 0.00)
    beta_nl_source <- c(I_X4_2 = 0.65, I_X1X3 = 0.00, I_sinX3 = 0.55)

    if (source_alignment == "aligned" && residual_type == "linear") {
      alpha_source <- 0.88
      linear_residual <- c(-0.18, 0.14, 0.10, -0.08, 0.06, 0.24, -0.20)
      nonlinear_residual <- c(I_X4_2 = 0.00, I_X1X3 = 0.00, I_sinX3 = 0.00)
    } else if (source_alignment == "aligned" && residual_type == "nonlinear") {
      alpha_source <- 0.95
      linear_residual <- c(-0.06, 0.04, 0.04, -0.03, 0.03, 0.10, -0.08)
      nonlinear_residual <- c(I_X4_2 = -0.45, I_X1X3 = 0.00, I_sinX3 = 0.78)
    } else if (source_alignment == "misaligned" && residual_type == "linear") {
      alpha_source <- 0.42
      linear_residual <- c(-0.34, 0.30, 0.12, -0.10, 0.08, 0.46, -0.40)
      nonlinear_residual <- c(I_X4_2 = 0.00, I_X1X3 = 0.00, I_sinX3 = 0.00)
    } else {
      alpha_source <- 0.50
      linear_residual <- c(-0.24, 0.20, -0.10, 0.08, -0.06, 0.18, -0.16)
      nonlinear_residual <- c(I_X4_2 = -0.55, I_X1X3 = 0.00, I_sinX3 = 1.05)
    }
  } else if (source_alignment == "aligned" && residual_type == "linear") {
    alpha_source <- 0.75
    linear_residual <- c(-0.25, 0.20, 0.12, -0.08, 0.06, 0.32, -0.28)
    nonlinear_residual <- c(I_X4_2 = 0.00, I_X1X3 = 0.00, I_sinX3 = 0.00)
  } else if (source_alignment == "aligned" && residual_type == "nonlinear") {
    alpha_source <- 0.80
    linear_residual <- c(-0.25, 0.20, 0.12, -0.08, 0.06, 0.32, -0.28)
    # Main-effect nonlinear residuals match the GAM working model; no interaction term.
    nonlinear_residual <- c(I_X4_2 = -0.22, I_X1X3 = 0.00, I_sinX3 = 0.35)
  } else if (source_alignment == "misaligned" && residual_type == "linear") {
    alpha_source <- 0.40
    linear_residual <- c(-0.38, 0.32, 0.14, -0.09, 0.07, 0.50, -0.45)
    nonlinear_residual <- c(I_X4_2 = 0.00, I_X1X3 = 0.00, I_sinX3 = 0.00)
  } else {
    alpha_source <- 0.55
    linear_residual <- c(-0.50, 0.42, -0.25, 0.20, -0.15, 0.25, -0.22)
    # Weak linear target signal plus strong smooth residual keeps Cox from dominating.
    nonlinear_residual <- c(I_X4_2 = -0.10, I_X1X3 = 0.00, I_sinX3 = 1.10)
  }

  setting$truth_source <- default_truth7(
    use_nonlinear_in_truth = TRUE,
    intercept = 0,
    beta_linear = beta_linear_source,
    beta_nonlinear = beta_nl_source
  )

  beta_linear_target <- alpha_source * beta_linear_source + linear_residual
  beta_nl_target <- alpha_source * beta_nl_source + nonlinear_residual

  delta_linear <- beta_linear_target - beta_linear_source
  nl_names <- union(names(beta_nl_source), names(beta_nl_target))
  src_nl_pad <- setNames(rep(0, length(nl_names)), nl_names)
  tgt_nl_pad <- setNames(rep(0, length(nl_names)), nl_names)
  src_nl_pad[names(beta_nl_source)] <- beta_nl_source
  tgt_nl_pad[names(beta_nl_target)] <- beta_nl_target

  setting$truth_target <- make_transfer_truth7(
    source_truth = setting$truth_source,
    delta_linear = delta_linear,
    delta_nonlinear = tgt_nl_pad - src_nl_pad,
    shift_scale = 1.0
  )

  setting$setting_name <- setting_name
  setting$nonlinear_type <- residual_type
  setting$shift_size <- source_alignment
  setting$dgp_type <- paste(source_alignment, residual_type, sep = "_")
  setting$alpha_source <- alpha_source
  setting$linear_residual <- linear_residual
  setting$nonlinear_residual <- nonlinear_residual
  setting$shift_design <- paste0(
    "target_eta=", alpha_source, "*source_eta+",
    residual_type, "_residual_with_linear_component"
  )
  setting
}

make_weibull_nonlinear_setting7 <- function(setting_name,
                                            nonlinear_type = c("few", "many"),
                                            shift_size = c("small", "large"),
                                            n_source,
                                            n_target_train,
                                            n_target_test,
                                            censor_source = 0.35,
                                            censor_target_train = 0.30,
                                            censor_target_test = 0.30) {
  nonlinear_type <- match.arg(nonlinear_type)
  shift_size <- match.arg(shift_size)

  setting <- make_recommended_setting7(
    n_source = n_source,
    n_target_train = n_target_train,
    n_target_test = n_target_test,
    shift_level = "moderate",
    use_nonlinear_in_truth = TRUE,
    time_model = "weibull",
    shift_scale = 1.0,
    censor_source = censor_source,
    censor_target_train = censor_target_train,
    censor_target_test = censor_target_test
  )

  beta_linear_source <- c(0.95, -0.80, 0.55, -0.45, 0.35, 0.00, 0.00)
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
  src_nl_pad <- setNames(rep(0, length(nl_names)), nl_names)
  tgt_nl_pad <- setNames(rep(0, length(nl_names)), nl_names)
  src_nl_pad[names(beta_nl_source)] <- beta_nl_source
  tgt_nl_pad[names(beta_nl_target)] <- beta_nl_target

  setting$truth_target <- make_transfer_truth7(
    source_truth = setting$truth_source,
    delta_linear = delta_linear,
    delta_nonlinear = tgt_nl_pad - src_nl_pad,
    shift_scale = 1.0
  )

  setting$setting_name <- setting_name
  setting$nonlinear_type <- nonlinear_type
  setting$shift_size <- shift_size
  setting$alpha_source <- alpha_source
  setting$linear_residual <- linear_residual
  setting$shift_design <- shift_design
  setting
}

make_setting_list <- function(model) {
  if (model == "weibull") {
    source(file.path(r_dir, "sim_setup_weibull.R"))
    constructors <- list(
      aligned_linear = function() make_adaptive_dgp_setting7(
        "W_ADAPT_aligned_linear_residual", "linear", "aligned",
        args$n_s, args$n_t, args$n_test, time_model = "weibull"),
      aligned_nonlinear = function() make_adaptive_dgp_setting7(
        "W_ADAPT_aligned_nonlinear_residual", "nonlinear", "aligned",
        args$n_s, args$n_t, args$n_test, time_model = "weibull"),
      misaligned_linear = function() make_adaptive_dgp_setting7(
        "W_ADAPT_misaligned_linear_residual", "linear", "misaligned",
        args$n_s, args$n_t, args$n_test, time_model = "weibull"),
      misaligned_nonlinear = function() make_adaptive_dgp_setting7(
        "W_ADAPT_misaligned_nonlinear_residual", "nonlinear", "misaligned",
        args$n_s, args$n_t, args$n_test, time_model = "weibull"),
      few_small = function() make_weibull_nonlinear_setting7(
        "W_NL_few_terms_small_shift", "few", "small", args$n_s, args$n_t, args$n_test),
      few_large = function() make_weibull_nonlinear_setting7(
        "W_NL_few_terms_large_shift", "few", "large", args$n_s, args$n_t, args$n_test),
      many_small = function() make_weibull_nonlinear_setting7(
        "W_NL_many_terms_small_shift", "many", "small", args$n_s, args$n_t, args$n_test),
      many_large = function() make_weibull_nonlinear_setting7(
        "W_NL_many_terms_large_shift", "many", "large", args$n_s, args$n_t, args$n_test)
    )
    return(lapply(setting_keys, function(k) constructors[[k]]()))
  }

  if (model == "aftlogistic") {
    source(file.path(r_dir, "sim_setup_aftlogistic.R"))
    constructors <- list(
      aligned_linear = function() make_adaptive_dgp_setting7(
        "AFT_ADAPT_aligned_linear_residual", "linear", "aligned",
        args$n_s, args$n_t, args$n_test, time_model = "aft_logistic"),
      aligned_nonlinear = function() make_adaptive_dgp_setting7(
        "AFT_ADAPT_aligned_nonlinear_residual", "nonlinear", "aligned",
        args$n_s, args$n_t, args$n_test, time_model = "aft_logistic"),
      misaligned_linear = function() make_adaptive_dgp_setting7(
        "AFT_ADAPT_misaligned_linear_residual", "linear", "misaligned",
        args$n_s, args$n_t, args$n_test, time_model = "aft_logistic"),
      misaligned_nonlinear = function() make_adaptive_dgp_setting7(
        "AFT_ADAPT_misaligned_nonlinear_residual", "nonlinear", "misaligned",
        args$n_s, args$n_t, args$n_test, time_model = "aft_logistic"),
      few_small = function() make_setting_AFT_NL_FewSmall(
        n_source = args$n_s, n_target_train = args$n_t, n_target_test = args$n_test),
      few_large = function() make_setting_AFT_NL_FewLarge(
        n_source = args$n_s, n_target_train = args$n_t, n_target_test = args$n_test),
      many_small = function() make_setting_AFT_NL_ManySmall(
        n_source = args$n_s, n_target_train = args$n_t, n_target_test = args$n_test),
      many_large = function() make_setting_AFT_NL_ManyLarge(
        n_source = args$n_s, n_target_train = args$n_t, n_target_test = args$n_test)
    )
    return(lapply(setting_keys, function(k) constructors[[k]]()))
  }

  stop("Unknown simulation model: ", model)
}

truth_eta_on_covariates <- function(dfX, truth_obj) {
  X <- make_design_features7(dfX, use_nonlinear = truth_obj$use_nonlinear_in_truth)
  X <- cbind(Intercept = 1, X)
  beta <- truth_obj$beta_true
  as.vector(X[, names(beta), drop = FALSE] %*% beta)
}

make_eta_alignment_one_setting <- function(setting, model, n_check = 5000, seed = 202605) {
  set.seed(seed)
  x_target <- gen_covariates7(
    n = n_check,
    dataset = "target",
    cov_shift = setting$cov_shift
  )
  eta_source_on_target <- truth_eta_on_covariates(x_target, setting$truth_source)
  eta_target_on_target <- truth_eta_on_covariates(x_target, setting$truth_target)

  data.frame(
    sim_model = model,
    setting_name = setting$setting_name,
    dgp_type = if (is.null(setting$dgp_type)) NA_character_ else setting$dgp_type,
    nonlinear_type = setting$nonlinear_type,
    shift_size = setting$shift_size,
    true_alpha_source = if (is.null(setting$alpha_source)) NA_real_ else setting$alpha_source,
    shift_design = if (is.null(setting$shift_design)) NA_character_ else setting$shift_design,
    eta_cor_source_target_on_target_x = stats::cor(eta_source_on_target, eta_target_on_target),
    eta_source_mean = mean(eta_source_on_target),
    eta_target_mean = mean(eta_target_on_target),
    eta_source_sd = stats::sd(eta_source_on_target),
    eta_target_sd = stats::sd(eta_target_on_target),
    eta_source_q05 = unname(stats::quantile(eta_source_on_target, 0.05)),
    eta_source_q50 = unname(stats::quantile(eta_source_on_target, 0.50)),
    eta_source_q95 = unname(stats::quantile(eta_source_on_target, 0.95)),
    eta_target_q05 = unname(stats::quantile(eta_target_on_target, 0.05)),
    eta_target_q50 = unname(stats::quantile(eta_target_on_target, 0.50)),
    eta_target_q95 = unname(stats::quantile(eta_target_on_target, 0.95)),
    stringsAsFactors = FALSE
  )
}

make_eta_alignment_table <- function(setting_list, model) {
  rows <- lapply(seq_along(setting_list), function(i) {
    make_eta_alignment_one_setting(
      setting = setting_list[[i]],
      model = model,
      seed = args$base_seed + 700000 + i
    )
  })
  do.call(rbind, rows)
}

make_error_row <- function(model, setting, rep_id, tq, t_eval, method, error_message,
                           sim_data) {
  D_train <- as.integer(sim_data$target_train$Y <= t_eval & sim_data$target_train$delta == 1)
  D_test <- as.integer(sim_data$target_test$Y <= t_eval & sim_data$target_test$delta == 1)
  data.frame(
    sim_model = model,
    setting_name = setting$setting_name,
    dgp_type = if (is.null(setting$dgp_type)) NA_character_ else setting$dgp_type,
    nonlinear_type = setting$nonlinear_type,
    shift_size = setting$shift_size,
    true_alpha_source = if (is.null(setting$alpha_source)) NA_real_ else setting$alpha_source,
    rep = rep_id,
    tq = tq,
    t_eval_rule = format_t_eval_rule(tq),
    method = method,
    source_alpha = NA_real_,
    selected_shift = NA_character_,
    selected_target = NA_character_,
    selected_calibration = NA_character_,
    auto_score_source = NA_real_,
    auto_score_linear = NA_real_,
    auto_score_gam = NA_real_,
    auto_correction_gain = NA_real_,
    auto_gam_gain = NA_real_,
    auto_gam_gain_se = NA_real_,
    auto_score_no_cal = NA_real_,
    auto_score_cal = NA_real_,
    auto_cal_loss_diff = NA_real_,
    auto_cal_loss_diff_se = NA_real_,
    auc_t = NA_real_,
    brier_t = NA_real_,
    log_loss_t = NA_real_,
    uno_c = NA_real_,
    harrell_c = NA_real_,
    t_eval = t_eval,
    n_source = setting$n_source,
    n_target_train = setting$n_target_train,
    n_target_test = setting$n_target_test,
    n_event_train = sum(D_train == 1),
    n_nonevent_train = sum(D_train == 0),
    n_event_test = sum(D_test == 1),
    n_nonevent_test = sum(D_test == 0),
    gamma_source = args$gamma_source,
    gamma_target = args$gamma_target,
    select_source = args$select_source,
    select_target = args$select_target,
    use_interaction = args$use_interaction,
    include_binary = TRUE,
    error = error_message,
    stringsAsFactors = FALSE
  )
}

make_result_row <- function(model, setting, rep_id, tq, t_eval, method, fit_obj,
                            sim_data) {
  row <- make_error_row(model, setting, rep_id, tq, t_eval, method, NA_character_, sim_data)
  row$source_alpha <- if (is.null(fit_obj$source_alpha)) NA_real_ else fit_obj$source_alpha
  row$selected_shift <- if (is.null(fit_obj$selected_shift)) NA_character_ else fit_obj$selected_shift
  row$selected_target <- if (is.null(fit_obj$selected_target)) NA_character_ else fit_obj$selected_target
  row$selected_calibration <- if (is.null(fit_obj$selected_calibration)) NA_character_ else fit_obj$selected_calibration
  row$auto_score_source <- if (is.null(fit_obj$auto_score_source)) NA_real_ else fit_obj$auto_score_source
  row$auto_score_linear <- if (is.null(fit_obj$auto_score_linear)) NA_real_ else fit_obj$auto_score_linear
  row$auto_score_gam <- if (is.null(fit_obj$auto_score_gam)) NA_real_ else fit_obj$auto_score_gam
  row$auto_correction_gain <- if (is.null(fit_obj$auto_correction_gain)) NA_real_ else fit_obj$auto_correction_gain
  row$auto_gam_gain <- if (is.null(fit_obj$auto_gam_gain)) NA_real_ else fit_obj$auto_gam_gain
  row$auto_gam_gain_se <- if (is.null(fit_obj$auto_gam_gain_se)) NA_real_ else fit_obj$auto_gam_gain_se
  row$auto_score_no_cal <- if (is.null(fit_obj$auto_score_no_cal)) NA_real_ else fit_obj$auto_score_no_cal
  row$auto_score_cal <- if (is.null(fit_obj$auto_score_cal)) NA_real_ else fit_obj$auto_score_cal
  row$auto_cal_loss_diff <- if (is.null(fit_obj$auto_cal_loss_diff)) NA_real_ else fit_obj$auto_cal_loss_diff
  row$auto_cal_loss_diff_se <- if (is.null(fit_obj$auto_cal_loss_diff_se)) NA_real_ else fit_obj$auto_cal_loss_diff_se
  row$auc_t <- if (is.null(fit_obj$auc_t)) NA_real_ else fit_obj$auc_t
  row$brier_t <- if (is.null(fit_obj$brier_t)) NA_real_ else fit_obj$brier_t
  row$log_loss_t <- if (is.null(fit_obj$log_loss_t)) NA_real_ else fit_obj$log_loss_t
  row$uno_c <- if (is.null(fit_obj$uno_c)) NA_real_ else fit_obj$uno_c
  row$harrell_c <- if (is.null(fit_obj$harrell_c)) NA_real_ else fit_obj$harrell_c
  row$error <- NA_character_
  row
}

run_one_sim_task <- function(task, model) {
  setting <- task$setting
  rep_id <- task$rep_id
  sim_data <- generate_from_setting7(setting = setting, seed = args$base_seed + 100000 * task$setting_id + rep_id)
  beta_target_ref <- extract_linear_beta_ref(setting$truth_target)

  rows <- list()
  row_id <- 1L
  for (tq in tq_list) {
    t_eval <- as.numeric(stats::quantile(sim_data$target_train$Y, probs = tq, na.rm = TRUE))
    for (method in methods_to_run) {
      fit_obj <- tryCatch({
        run_one_method7(
          method = method,
          sim_data = sim_data,
          t = t_eval,
          source_covariates = paste0("X", 1:7),
          target_covariates = paste0("X", 1:7),
          common_cov = paste0("X", 1:7),
          target_only_cov = character(0),
          lambda_w = 1,
          eps_G = 1e-3,
          w_cap = 20,
          max_folds = args$max_folds,
          beta_target_ref = beta_target_ref,
          gamma_source = args$gamma_source,
          gamma_target = args$gamma_target,
          select_source = args$select_source,
          select_target = args$select_target,
          use_interaction = args$use_interaction,
          include_binary = TRUE
        )
      }, error = function(e) {
        structure(list(error_message = conditionMessage(e)), class = "simulation_error")
      })

      rows[[row_id]] <- if (inherits(fit_obj, "simulation_error")) {
        make_error_row(model, setting, rep_id, tq, t_eval, method, fit_obj$error_message, sim_data)
      } else {
        make_result_row(model, setting, rep_id, tq, t_eval, method, fit_obj, sim_data)
      }
      row_id <- row_id + 1L
    }
  }

  object <- NULL
  if (isTRUE(args$save_full_objects)) {
    object <- list(setting = setting, rep = rep_id, sim_data = sim_data)
  }

  list(task_id = task$task_id, summary_df = do.call(rbind, rows), object = object)
}

make_auto_choice_summary <- function(df) {
  auto <- df[df$method %in% AUTO_METHODS, , drop = FALSE]
  if (nrow(auto) == 0) return(data.frame())
  aggregate(
    rep ~ sim_model + setting_name + tq + t_eval_rule + method + selected_shift + selected_target + selected_calibration,
    data = auto,
    FUN = length
  )
}

run_model <- function(model) {
  setting_list <- make_setting_list(model)
  eta_alignment <- make_eta_alignment_table(setting_list, model)
  rep_ids <- seq.int(args$rep_start, args$rep_end)
  task_grid <- expand.grid(setting_id = seq_along(setting_list), rep_id = rep_ids)
  task_list <- lapply(seq_len(nrow(task_grid)), function(i) {
    setting_id <- task_grid$setting_id[i]
    list(
      task_id = i,
      setting_id = setting_id,
      rep_id = task_grid$rep_id[i],
      setting = setting_list[[setting_id]]
    )
  })

  out_dir <- if (nzchar(args$out_dir)) {
    if (exists("models", inherits = TRUE) && length(models) > 1) file.path(args$out_dir, model) else args$out_dir
  } else {
    file.path(args$out_root, paste0("simulation_", model))
  }
  out_dir <- prepare_out_dir(out_dir, overwrite = TRUE)
  write_run_metadata(out_dir, args, methods_to_run, tq_list)

  n_cores <- max(1L, min(args$cores, length(task_list)))
  use_parallel <- !isTRUE(args$no_parallel) && n_cores > 1L && .Platform$OS.type != "windows"

  cat("====================================================\n")
  cat("Simulation model:", model, "\n")
  cat("Output:", out_dir, "\n")
  cat("Settings:", paste(vapply(setting_list, `[[`, character(1), "setting_name"), collapse = ", "), "\n")
  cat("Eta alignment check:\n")
  print(eta_alignment)
  cat("Methods:", paste(methods_to_run, collapse = ", "), "\n")
  cat("tq:", paste(tq_list, collapse = ", "), "\n")
  cat("Rep ids:", args$rep_start, "to", args$rep_end, "\n")
  cat("RSF trees:", args$rsf_trees, "| RSF min.node.size:", args$rsf_min_node_size,
      "| RSF threads:", args$rsf_threads, "| XGBoost threads:", args$xgb_threads, "\n")
  cat("Tasks:", length(task_list), "| Parallel:", use_parallel, "| Cores:", n_cores, "\n")
  cat("====================================================\n")

  task_results <- if (use_parallel) {
    parallel::mclapply(task_list, run_one_sim_task, model = model, mc.cores = n_cores)
  } else {
    lapply(task_list, function(task) run_one_sim_task(task, model))
  }
  task_results <- task_results[order(vapply(task_results, `[[`, integer(1), "task_id"))]

  all_results_df <- do.call(rbind, lapply(task_results, `[[`, "summary_df"))
  summary_table <- make_metric_summary(
    all_results_df,
    group_cols = c(
      "sim_model", "setting_name", "dgp_type", "nonlinear_type",
      "shift_size", "true_alpha_source", "tq", "t_eval_rule", "method"
    )
  )
  summary_table <- summary_table[order(summary_table$setting_name, summary_table$tq, -summary_table$auc_t_mean), ]

  auto_choice_summary <- make_auto_choice_summary(all_results_df)

  saveRDS(all_results_df, file.path(out_dir, paste0("all_results_", model, "_nonlinear.rds")))
  write.csv(all_results_df, file.path(out_dir, paste0("all_results_", model, "_nonlinear.csv")), row.names = FALSE)
  saveRDS(summary_table, file.path(out_dir, paste0("summary_table_", model, "_nonlinear.rds")))
  write.csv(summary_table, file.path(out_dir, paste0("summary_table_", model, "_nonlinear.csv")), row.names = FALSE)
  write.csv(eta_alignment, file.path(out_dir, paste0("eta_alignment_", model, ".csv")), row.names = FALSE)
  write.csv(auto_choice_summary, file.path(out_dir, paste0("auto_choice_summary_", model, "_nonlinear.csv")), row.names = FALSE)

  if (isTRUE(args$save_full_objects)) {
    all_objects <- lapply(task_results, `[[`, "object")
    saveRDS(all_objects, file.path(out_dir, paste0("all_objects_", model, "_nonlinear.rds")))
  }

  cat("\n==================== FINAL SUMMARY:", model, "====================\n")
  print(summary_table)
  cat("Saved files to:", out_dir, "\n")
}

model_arg <- tolower(args$model)
models <- switch(
  model_arg,
  weibull = "weibull",
  aft = "aftlogistic",
  aftlogistic = "aftlogistic",
  both = c("weibull", "aftlogistic"),
  stop("Unknown --model: ", args$model)
)

for (model in models) {
  run_model(model)
}
