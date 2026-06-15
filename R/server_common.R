local_r_lib <- file.path(getwd(), "_r_lib")
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

STAT_BASELINE_METHODS <- c(
  "source_only_logistic",
  "target_only_logistic",
  "pooled_logistic",
  "source_only_gam",
  "target_only_gam",
  "pooled_gam",
  "source_only_cox",
  "target_only_cox",
  "pooled_cox"
)

ML_BASELINE_METHODS <- c(
  "source_only_rsf",
  "target_only_rsf",
  "pooled_rsf",
  "source_only_ipcw_xgb",
  "target_only_ipcw_xgb",
  "pooled_ipcw_xgb"
)

RSF_BASELINE_METHODS <- c("source_only_rsf", "target_only_rsf", "pooled_rsf")
XGB_BASELINE_METHODS <- c("source_only_ipcw_xgb", "target_only_ipcw_xgb", "pooled_ipcw_xgb")

NONTRANSFER_BASELINE_METHODS <- c(STAT_BASELINE_METHODS, ML_BASELINE_METHODS)

TRANSFER_FIXED_METHODS <- c(
  "source_gam_recal_intercept",
  "source_gam_recal_slope",
  "trans_linear_linear",
  "trans_linear_gam",
  "trans_linear_gam_calibrated",
  "trans_gam_linear",
  "trans_gam_linear_calibrated",
  "trans_gam_gam_no_tilt",
  "trans_gam_gam_calibrated_no_tilt",
  "trans_gam_gam_partial_cal_no_tilt",
  "trans_gam_gam_penalized_cal_no_tilt",
  "trans_gam_gam",
  "trans_gam_gam_calibrated",
  "trans_gam_gam_partial_cal",
  "trans_gam_gam_penalized_cal"
)

SCHEME1_METHODS <- c(
  "trans_gam_gam_select_no_tilt",
  "trans_gam_gam_select"
)

GAM_GAM_TP_METHODS <- c(
  "trans_gam_gam_tp_no_tilt",
  "trans_gam_gam_tp"
)

GAM_GAM_BASIS_METHODS <- c(
  "trans_gam_gam_tp_no_tilt",
  "trans_gam_gam_tp",
  "trans_gam_gam_source_tp_no_tilt",
  "trans_gam_gam_source_tp",
  "trans_gam_gam_all_tp_no_tilt",
  "trans_gam_gam_all_tp"
)

AUTO_METHODS <- c(
  "trans_gam_auto_calibrated",
  "trans_gam_auto_calibrated_tilt",
  "trans_gam_auto_no_cal",
  "trans_gam_auto_no_cal_tilt",
  "trans_gam_auto_1se_no_cal",
  "trans_gam_auto_1se_calibrated",
  "trans_gam_auto_two_layer",
  "trans_gam_linear_auto_cal",
  "trans_gam_gam_auto_cal"
)

GAM_GAM_CALIBRATION_METHODS <- c(
  "trans_gam_gam_no_tilt",
  "trans_gam_gam_select_no_tilt",
  "trans_gam_gam_calibrated_no_tilt",
  "trans_gam_gam_partial_cal_no_tilt",
  "trans_gam_gam_penalized_cal_no_tilt",
  "trans_gam_gam_auto_cal",
  "trans_gam_gam",
  "trans_gam_gam_select",
  "trans_gam_gam_calibrated",
  "trans_gam_gam_partial_cal",
  "trans_gam_gam_penalized_cal"
)

GAM_GAM_PENALIZED_CAL_METHODS <- c(
  "trans_gam_gam_penalized_cal_no_tilt",
  "trans_gam_gam_penalized_cal"
)

GAM_GAM_CONVEX_CAL_METHODS <- c(
  "trans_gam_gam_partial_cal_no_tilt",
  "trans_gam_gam_partial_cal"
)

TRANSCox_METHODS <- c("transcox_bic", "transcox_fixed")

SIMULATION_METHODS_ALL <- c(
  NONTRANSFER_BASELINE_METHODS,
  TRANSFER_FIXED_METHODS,
  SCHEME1_METHODS,
  GAM_GAM_BASIS_METHODS,
  AUTO_METHODS,
  "oracle"
)

REALDATA_METHODS_ALL <- c(
  NONTRANSFER_BASELINE_METHODS,
  TRANSFER_FIXED_METHODS,
  SCHEME1_METHODS,
  GAM_GAM_BASIS_METHODS,
  AUTO_METHODS
)

SIMULATION_METHODS_VALID <- c(SIMULATION_METHODS_ALL, TRANSCox_METHODS)
REALDATA_METHODS_VALID <- c(REALDATA_METHODS_ALL, TRANSCox_METHODS)

parse_token_list <- function(x) {
  x <- trimws(paste(x, collapse = ","))
  x <- chartr("{}[]()", "      ", x)
  x <- gsub("\uFF0C", ",", x, fixed = TRUE)
  x <- gsub("\uFF1B", ";", x, fixed = TRUE)
  pieces <- trimws(unlist(strsplit(x, "[,;[:space:]]+")))
  pieces[nzchar(pieces)]
}

parse_tq_list <- function(x) {
  pieces <- parse_token_list(x)
  tq <- suppressWarnings(as.numeric(pieces))

  if (length(tq) == 0 || any(!is.finite(tq))) {
    stop("tq must contain numeric quantiles, e.g. {0.25}, {0.5,0.75}.")
  }
  if (any(tq <= 0 | tq >= 1)) {
    stop("tq values must be strictly between 0 and 1.")
  }

  unique(tq)
}

format_t_eval_rule <- function(tq) {
  paste0("target_q", formatC(round(100 * tq), width = 2, flag = "0", format = "d"))
}

parse_methods <- function(x, include_oracle = FALSE) {
  tokens <- parse_token_list(x)
  if (length(tokens) == 0 || any(tokens == "all")) {
    return(if (include_oracle) SIMULATION_METHODS_ALL else REALDATA_METHODS_ALL)
  }

  out <- character(0)
  for (token in tokens) {
    expanded <- switch(
      token,
      baseline = NONTRANSFER_BASELINE_METHODS,
      nontransfer = NONTRANSFER_BASELINE_METHODS,
      stat = STAT_BASELINE_METHODS,
      standard = STAT_BASELINE_METHODS,
      ml = ML_BASELINE_METHODS,
      machine_learning = ML_BASELINE_METHODS,
      rsf = RSF_BASELINE_METHODS,
      xgboost = XGB_BASELINE_METHODS,
      xgb = XGB_BASELINE_METHODS,
      transfer = TRANSFER_FIXED_METHODS,
      fixed_transfer = TRANSFER_FIXED_METHODS,
      scheme1 = SCHEME1_METHODS,
      selective_gam = SCHEME1_METHODS,
      gam_gam_select = SCHEME1_METHODS,
      gam_gam_tp = GAM_GAM_TP_METHODS,
      gam_gam_basis = GAM_GAM_BASIS_METHODS,
      basis_ablation = GAM_GAM_BASIS_METHODS,
      tp = GAM_GAM_TP_METHODS,
      calibration = GAM_GAM_CALIBRATION_METHODS,
      cal = GAM_GAM_CALIBRATION_METHODS,
      gam_gam_calibration = GAM_GAM_CALIBRATION_METHODS,
      penalized_calibration = GAM_GAM_PENALIZED_CAL_METHODS,
      penalized_auto = GAM_GAM_PENALIZED_CAL_METHODS,
      convex_calibration = GAM_GAM_CONVEX_CAL_METHODS,
      convex_cal = GAM_GAM_CONVEX_CAL_METHODS,
      transcox = "transcox_bic",
      transcox_bic = "transcox_bic",
      transcox_fixed = "transcox_fixed",
      auto = AUTO_METHODS,
      ablation = c(TRANSFER_FIXED_METHODS, AUTO_METHODS),
      oracle = if (include_oracle) "oracle" else stop("oracle is not available for real data."),
      token
    )
    out <- c(out, expanded)
  }

  valid <- if (include_oracle) SIMULATION_METHODS_VALID else REALDATA_METHODS_VALID
  bad <- setdiff(out, valid)
  if (length(bad) > 0) {
    stop("Unknown method(s): ", paste(bad, collapse = ", "))
  }

  unique(out)
}

prepare_out_dir <- function(out_dir, overwrite = TRUE) {
  if (!nzchar(out_dir)) stop("out_dir cannot be empty.")

  normalized <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (normalized %in% c("/", cwd, paste0(cwd, "/"))) {
    stop("Refusing to overwrite root/current working directory: ", normalized,
         ". Use a dedicated output subdirectory.")
  }

  if (dir.exists(out_dir) && isTRUE(overwrite)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(out_dir, winslash = "/", mustWork = TRUE)
}

summarise_metric <- function(df, group_cols, metric) {
  groups <- unique(df[, group_cols, drop = FALSE])
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, , drop = FALSE]
    keep <- rep(TRUE, nrow(df))
    for (col in group_cols) keep <- keep & df[[col]] == g[[col]]

    x <- df[[metric]][keep]
    x <- x[is.finite(x)]

    out <- g
    out[[paste0(metric, "_mean")]] <- if (length(x) > 0) mean(x) else NA_real_
    out[[paste0(metric, "_sd")]] <- if (length(x) > 1) stats::sd(x) else NA_real_
    out
  })

  do.call(rbind, rows)
}

make_metric_summary <- function(df, group_cols,
                                metrics = c("auc_t", "brier_t", "log_loss_t", "uno_c", "harrell_c")) {
  metrics <- metrics[metrics %in% names(df)]
  metric_tables <- lapply(metrics, function(m) summarise_metric(df, group_cols, m))
  Reduce(function(x, y) merge(x, y, by = group_cols, all = TRUE), metric_tables)
}

write_run_metadata <- function(out_dir, args, methods, tq_list) {
  writeLines(methods, con = file.path(out_dir, "methods_used.txt"))
  writeLines(as.character(tq_list), con = file.path(out_dir, "tq_used.txt"))
  capture.output(str(args), file = file.path(out_dir, "args_used.txt"))
}
