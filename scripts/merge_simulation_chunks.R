cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_file <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1])) else normalizePath("scripts/merge_simulation_chunks.R", mustWork = FALSE)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "R", "server_common.R"))

parse_cli <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  if (any(raw %in% c("-h", "--help")) || length(raw) == 0) {
    cat("Usage: Rscript scripts/merge_simulation_chunks.R --chunk-root results/weibull_300/chunks --out-dir results/weibull_300/combined\n")
    quit(save = "no", status = if (length(raw) == 0) 1 else 0)
  }
  out <- list(chunk_root = "", out_dir = "")
  i <- 1L
  while (i <= length(raw)) {
    key <- raw[i]
    if (!grepl("^--", key)) stop("Unexpected token: ", key)
    key <- gsub("-", "_", sub("^--", "", key), fixed = TRUE)
    if (!key %in% names(out)) stop("Unknown option --", gsub("_", "-", key))
    if (i == length(raw)) stop("Missing value for --", gsub("_", "-", key))
    out[[key]] <- raw[i + 1L]
    i <- i + 2L
  }
  if (!nzchar(out$chunk_root)) stop("--chunk-root is required.")
  if (!nzchar(out$out_dir)) stop("--out-dir is required.")
  out
}

args <- parse_cli()
if (!dir.exists(args$chunk_root)) stop("Cannot find chunk root: ", args$chunk_root)
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(args$chunk_root, pattern = "^all_results_.*_nonlinear\\.csv$",
                    recursive = TRUE, full.names = TRUE)
if (length(files) == 0) {
  stop("No all_results_*_nonlinear.csv files found under ", args$chunk_root)
}

pieces <- lapply(files, read.csv, stringsAsFactors = FALSE)
all_results <- do.call(rbind, pieces)
all_results <- all_results[order(all_results$sim_model, all_results$setting_name,
                                 all_results$rep, all_results$tq, all_results$method), ]

model_tag <- if (length(unique(all_results$sim_model)) == 1) unique(all_results$sim_model) else "combined"
out_csv <- file.path(args$out_dir, paste0("all_results_", model_tag, "_nonlinear.csv"))
out_rds <- file.path(args$out_dir, paste0("all_results_", model_tag, "_nonlinear.rds"))
write.csv(all_results, out_csv, row.names = FALSE)
saveRDS(all_results, out_rds)

summary_table <- make_metric_summary(
  all_results,
  group_cols = c(
    "sim_model", "setting_name", "dgp_type", "nonlinear_type",
    "shift_size", "true_alpha_source", "tq", "t_eval_rule", "method"
  )
)
summary_csv <- file.path(args$out_dir, paste0("summary_table_", model_tag, "_nonlinear.csv"))
summary_rds <- file.path(args$out_dir, paste0("summary_table_", model_tag, "_nonlinear.rds"))
write.csv(summary_table, summary_csv, row.names = FALSE)
saveRDS(summary_table, summary_rds)

cat("Merged", length(files), "chunk file(s).\n")
cat("Wrote:", out_csv, "\n")
cat("Wrote:", summary_csv, "\n")
