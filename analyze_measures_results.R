## Analyze measures_results.csv from run_simulation_parallel.R.
##
## Output:
##   measures_analysis_summary.csv
##
## The summary is grouped by scale, method, contrast, and measure.
## - scale == "log": log(RR/OR/SR) scale; this is the primary scale because
##   se, ci_low, and ci_up in measures_results.csv are stored on this scale.
## - scale == "exp": RR/OR/SR scale; SE is computed by the delta method:
##   se_exp = estimate * se_log.

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, cmd_args, fixed = TRUE)

  if (length(match) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", cmd_args[match[1]]),
                                 winslash = "/", mustWork = TRUE)))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

first_existing_file <- function(paths) {
  paths <- paths[nzchar(paths)]
  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

ensure_output_file <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  out_dir <- dirname(path)

  ok <- dir.exists(out_dir) || suppressWarnings(dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  ))

  if (ok) {
    test_file <- tempfile("write_test_", tmpdir = out_dir)
    ok <- tryCatch({
      writeLines("ok", test_file)
      unlink(test_file)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
  }

  if (isTRUE(ok)) return(path)

  fallback_dir <- file.path(Sys.getenv("TEMP", unset = tempdir()),
                            "gop_two_modules_parallel")
  fallback_dir <- normalizePath(fallback_dir, winslash = "/", mustWork = FALSE)
  ok <- dir.exists(fallback_dir) || dir.create(fallback_dir, recursive = TRUE)
  if (!ok) stop("Cannot create output directory: ", fallback_dir, call. = FALSE)

  fallback <- file.path(fallback_dir, basename(path))
  warning(
    "Cannot write to requested output file: ", path,
    ". Using temporary output file instead: ", fallback,
    call. = FALSE
  )
  fallback
}

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x)
}

get_true_value <- function(x) {
  x <- x[is.finite(x)]
  ux <- unique(round(x, 12))
  if (length(ux) == 0) return(NA_real_)
  if (length(ux) > 1) {
    warning("Multiple true values found in one group; using their mean.",
            call. = FALSE)
    return(mean(x))
  }
  x[1]
}

est_result <- function(df, estimate_col, se_col, low_col, up_col, true_col,
                       p_col = "p_value") {
  est <- df[[estimate_col]]
  se <- df[[se_col]]
  low <- df[[low_col]]
  up <- df[[up_col]]
  para_true <- get_true_value(df[[true_col]])

  keep <- is.finite(est) & is.finite(se) & is.finite(low) & is.finite(up) &
    is.finite(para_true)
  est <- est[keep]
  se <- se[keep]
  low <- low[keep]
  up <- up[keep]

  n <- length(est)
  mean_est <- safe_mean(est)
  mean_se <- safe_mean(se)
  mcsd <- safe_sd(est)

  p_value <- if (p_col %in% names(df)) df[[p_col]][keep] else rep(NA_real_, n)

  data.frame(
    n_sim = n,
    true_value = para_true,
    mean_estimate = mean_est,
    bias = mean_est - para_true,
    mean_se = mean_se,
    mean_se_over_sqrt_n = mean_se / sqrt(n),
    mcsd = mcsd,
    se_mcsd_ratio = mean_se / mcsd,
    coverage = safe_mean((low < para_true) & (up > para_true)),
    p_le_0.05 = safe_mean(p_value <= 0.05),
    row.names = NULL
  )
}

make_exp_scale_columns <- function(df) {
  df$se_exp <- df$estimate * df$se
  df
}

summarize_measures <- function(df) {
  required <- c(
    "method", "contrast", "measure",
    "log_estimate", "se", "ci_low", "ci_up", "true_log_effect",
    "estimate", "ci_low_exp", "ci_up_exp", "true_effect", "p_value"
  )
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("measures_results.csv is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  df <- make_exp_scale_columns(df)
  groups <- split(df, interaction(df$method, df$contrast, df$measure,
                                  drop = TRUE, lex.order = TRUE))

  out <- lapply(groups, function(g) {
    group_info <- data.frame(
      method = g$method[1],
      contrast = g$contrast[1],
      measure = g$measure[1],
      convergence_rate = if ("converged" %in% names(g)) {
        safe_mean(to_logical(g$converged))
      } else {
        NA_real_
      },
      row.names = NULL
    )

    log_result <- est_result(
      g,
      estimate_col = "log_estimate",
      se_col = "se",
      low_col = "ci_low",
      up_col = "ci_up",
      true_col = "true_log_effect"
    )
    log_result <- cbind(scale = "log", group_info, log_result)

    exp_result <- est_result(
      g,
      estimate_col = "estimate",
      se_col = "se_exp",
      low_col = "ci_low_exp",
      up_col = "ci_up_exp",
      true_col = "true_effect"
    )
    exp_result <- cbind(scale = "exp", group_info, exp_result)

    rbind(log_result, exp_result)
  })

  summary <- do.call(rbind, out)
  rownames(summary) <- NULL
  summary[order(summary$scale, summary$method, summary$contrast), ]
}

run_measures_analysis <- function(input_file = NULL, output_file = NULL) {
  project_dir <- get_script_dir()

  input_file_env <- Sys.getenv("GOP_MEASURES_FILE", unset = "")
  output_file_env <- Sys.getenv("GOP_MEASURES_ANALYSIS_OUTPUT", unset = "")

  candidate_input_files <- c(
    input_file,
    input_file_env,
    file.path(project_dir, "simulation_output_parallel", "measures_results.csv"),
    file.path(Sys.getenv("GOP_OUTPUT_DIR", unset = ""),
              "measures_results.csv"),
    file.path(Sys.getenv("TEMP", unset = tempdir()),
              "gop_two_modules_parallel", "measures_results.csv")
  )

  input_file <- first_existing_file(candidate_input_files)
  if (is.na(input_file)) {
    stop(
      "Cannot find measures_results.csv. Set GOP_MEASURES_FILE to its full path.",
      call. = FALSE
    )
  }

  if (is.null(output_file) || !nzchar(output_file)) {
    output_file <- if (nzchar(output_file_env)) {
      output_file_env
    } else {
      file.path(dirname(input_file), "measures_analysis_summary.csv")
    }
  }
  output_file <- ensure_output_file(output_file)

  cat("Reading: ", input_file, "\n", sep = "")
  measures_results <- read.csv(input_file, stringsAsFactors = FALSE)
  analysis_summary <- summarize_measures(measures_results)

  write.csv(analysis_summary, output_file, row.names = FALSE)

  cat("Wrote: ", output_file, "\n", sep = "")
  print(analysis_summary)
  invisible(analysis_summary)
}

if (!identical(Sys.getenv("GOP_RUN_ANALYSIS_MAIN", unset = "TRUE"), "FALSE")) {
  run_measures_analysis()
}
