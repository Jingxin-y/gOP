## Entry point for the two-time-point gOP RR/OR/OR functions.
## Keep using source("gop_two.R"); implementation lives in gop_two_modules/.

.gop_two_this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = FALSE),
  error = function(e) NA_character_
)
if (length(.gop_two_this_file) != 1 || is.na(.gop_two_this_file)) {
  .gop_two_this_file <- normalizePath("gop_two.R", winslash = "/", mustWork = FALSE)
}

.gop_two_module_dir <- file.path(dirname(.gop_two_this_file), "gop_two_modules")
if (!dir.exists(.gop_two_module_dir)) {
  .gop_two_module_dir <- file.path(getwd(), "gop_two_modules")
}

for (.gop_two_module in c(
  "01_utils.R",
  "02_probabilities.R",
  "03_likelihood.R",
  "04_variance.R",
  "05_estimator.R",
  "op_gop_measure_mle.R"
)) {
  source(file.path(.gop_two_module_dir, .gop_two_module))
}

rm(.gop_two_this_file, .gop_two_module_dir, .gop_two_module)
