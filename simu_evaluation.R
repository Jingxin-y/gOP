`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}


standardize_sim_data <- function(dat) {
  dat <- as.data.frame(dat)
  
  if (!"y" %in% names(dat) && "Y" %in% names(dat)) dat$y <- dat$Y
  if (!"a0" %in% names(dat) && "A" %in% names(dat)) dat$a0 <- dat$A
  if (!"a1" %in% names(dat) && "B" %in% names(dat)) dat$a1 <- dat$B
  
  if (!"z1" %in% names(dat) && "Xc1" %in% names(dat)) dat$z1 <- dat$Xc1
  if (!"z2" %in% names(dat) && "Xb" %in% names(dat)) dat$z2 <- dat$Xb
  if (!"z3" %in% names(dat) && "Xc2" %in% names(dat)) dat$z3 <- dat$Xc2
  
  needed <- c("y", "a0", "a1", "z1", "z2", "z3", "p00", "p10", "p01", "p11")
  missing <- setdiff(needed, names(dat))
  if (length(missing) > 0) {
    stop("Simulation data are missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  
  if (!"id" %in% names(dat)) dat$id <- seq_len(nrow(dat))
  dat
}

effect_from_pair <- function(p_num, p_den, measure, eps = 1e-8) {
  if (exists("measure_eta", mode = "function", inherits = TRUE)) {
    return(measure_eta(p_num, p_den, measure, eps = eps))
  }
  
  measure <- toupper(measure)
  
  if (measure == "OR") return(logit_local(p_num, eps = eps) - logit_local(p_den, eps = eps))
  if (measure == "RR") return(log(p_num / p_den))
  if (measure == "SR") return(log((1 - p_num) / (1 - p_den)))
  
  stop("Unsupported measure: ", measure, call. = FALSE)
}


selected_effects_from_prob <- function(prob, measures = c("RR", "OR", "OR"),
                                       eps = 1e-8) {
  prob <- as.matrix(prob)
  colnames(prob) <- c("p00", "p10", "p01", "p11")
  measures <- normalize_measures(measures)
  
  out <- cbind(
    eta_10_00 = effect_from_pair(prob[, "p10"], prob[, "p00"],
                                 measures[["p10_p00"]], eps = eps),
    eta_01_00 = effect_from_pair(prob[, "p01"], prob[, "p00"],
                                 measures[["p01_p00"]], eps = eps),
    eta_11_10 = effect_from_pair(prob[, "p11"], prob[, "p10"],
                                 measures[["p11_p10"]], eps = eps)
  )
  attr(out, "measures") <- measures
  out
}

make_prediction_object <- function(cell_prob, measures = c("RR", "OR", "OR"),
                                   eps = 1e-8) {
  cell_prob <- as.matrix(cell_prob)
  colnames(cell_prob) <- c("p00", "p10", "p01", "p11")
  measures <- normalize_measures(measures)
  
  list(
    cell_prob = cell_prob,
    selected_effects = selected_effects_from_prob(cell_prob, measures = measures,
                                                  eps = eps),
    measures = measures
  )
}

normalize_measures <- function(measures) {
  normalize.measures.Two(measures)
}


fit_rbc_method <- function(dat, eps = 1e-8) {
  if (!requireNamespace("rbc", quietly = TRUE)) {
    stop("The rbc package is required. Install it with install.packages('rbc').",
         call. = FALSE)
  }
  
  dat <- standardize_sim_data(dat)
  fit <- rbc::rbc(
    formula = y ~ 1 + z1 + z2 + z3 | 0 + a0 | 0 + a1 + a0:a1,
    init = rbc::Bernoulli(
      prob = mean(dat$y, na.rm = TRUE)
    ),
    flows = list(
      rbc::BinomialGLM(link = "logit"),
      rbc::ScaleRisk1,
      rbc::ScaleOdds
    ),
    data = dat,
    hessian = TRUE
  )
  
  fit_obj <- fit$fits$fits[[fit$fits$n_fits]]
  list(
    method = "RBC",
    object = fit,
    convergence = is.null(fit_obj$convergence) || fit_obj$convergence == 0,
    value = as.numeric(stats::logLik(fit))
  )
}

predict_rbc_method <- function(fit, newdata, measures = c("RR", "OR", "OR"),
                               eps = 1e-8) {
  newdata <- standardize_sim_data(newdata)
  nd00 <- nd10 <- nd01 <- nd11 <- newdata
  nd00$a0 <- 0; nd00$a1 <- 0
  nd10$a0 <- 1; nd10$a1 <- 0
  nd01$a0 <- 0; nd01$a1 <- 1
  nd11$a0 <- 1; nd11$a1 <- 1
  
  prob <- cbind(
    p00 = as.numeric(stats::predict(fit$object, newdata = nd00)),
    p10 = as.numeric(stats::predict(fit$object, newdata = nd10)),
    p01 = as.numeric(stats::predict(fit$object, newdata = nd01)),
    p11 = as.numeric(stats::predict(fit$object, newdata = nd11))
  )
  make_prediction_object(prob, measures = measures, eps = eps)
}

gop_formulas <- function() {
  list(
    f0 = ~ 1,
    f1 = ~ 1,
    f2 = ~ 1,
    f3 = ~ 1 + z1 + z2 + z3
  )
}

fit_gop_method <- function(dat, measures = c("RR", "OR", "OR"),
                           max.step = 1000, thres = 1e-8) {
  dat <- standardize_sim_data(dat)
  measures <- normalize_measures(measures)
  f <- gop_formulas()
  
  fit <- MLEst.TwoByMeasures(
    y = "y",
    a0 = "a0",
    a1 = "a1",
    x0 = f$f0,
    x1 = f$f1,
    x2 = f$f2,
    x3 = f$f3,
    data = dat,
    measures = measures,
    max.step = max.step,
    thres = thres
  )
  
  list(
    method = "gOP",
    object = fit,
    measures = measures,
    formulas = f,
    convergence = isTRUE(fit$convergence),
    value = fit$fit$value %||% NA_real_
  )
}
predict_gop_method <- function(fit, newdata, eps = 1e-8) {
  newdata <- standardize_sim_data(newdata)
  f <- fit$formulas
  x0 <- stats::model.matrix(f$f0, newdata)
  x1 <- stats::model.matrix(f$f1, newdata)
  x2 <- stats::model.matrix(f$f2, newdata)
  x3 <- stats::model.matrix(f$f3, newdata)
  
  point <- fit$object$point.est
  q0 <- ncol(x0)
  q1 <- ncol(x1)
  q2 <- ncol(x2)
  q3 <- ncol(x3)
  
  b0 <- point[seq_len(q0)]
  b1 <- point[q0 + seq_len(q1)]
  b2 <- point[q0 + q1 + seq_len(q2)]
  b3 <- point[q0 + q1 + q2 + seq_len(q3)]
  
  theta0 <- as.vector(x0 %*% b0)
  theta1 <- as.vector(x1 %*% b1)
  theta2 <- as.vector(x2 %*% b2)
  phi <- as.vector(x3 %*% b3)
  prob <- getTwo.by.measures(theta0, theta1, theta2, phi,
                             measures = fit$measures)
  colnames(prob) <- c("p00", "p10", "p01", "p11")
}


evalu_gOP_CM <- function(fit){
  point <- fit$object$point.est
  se <- fit$object$se.est
  CI.low <- fit$object$conf.lower
  CI.up <- fit$object$conf.upper
  p <- fit$object$p.value
  measure <- list(
    mea.point = point[1:3],
    mea.se = se[1:3],
    mea.CI.low = CI.low[1:3],
    mea.CI.up = CI.up[1:3],
    mea.p = p[1:3]
  )
  cellprob <- fit$object$prob.est
  list(measures = measure,
       cellprob = cellprob)
}

make_exp_ci <- function(log_est, var_log_est, conf = 0.95) {
  se_log <- sqrt(var_log_est)
  z <- qnorm(1 - (1 - conf) / 2)
  p.temp <- stats::pnorm(log_est / se_log, 0, 1)
  p.value <- 2 * pmin(p.temp, 1 - p.temp)
  data.frame(
    log_est = log_est,
    se_log = se_log,
    estimate = exp(log_est),
    lower = log_est - z * se_log,
    upper = log_est + z * se_log,
    p = p.value,
    var_effect = exp(log_est)^2 * var_log_est
  )
}

evalu_rbc_CM <- function(fit,pred){
  
  b <- coef(fit$object)
  V <- vcov(fit$object)
  
  idx_a0  <- 5
  idx_a1  <- 6
  idx_int <- 7
  
  log_RR_10_00 <- b[idx_a0]
  var_log_RR_10_00 <- V[idx_a0, idx_a0]
  
  log_OR_01_00 <- b[idx_a1]
  var_log_OR_01_00 <- V[idx_a1, idx_a1]
  
  log_OR_11_10 <- b[idx_a1] + b[idx_int]
  var_log_OR_11_10 <- 
    V[idx_a1, idx_a1] +
    V[idx_int, idx_int] +
    2 * V[idx_a1, idx_int]
  
  RR_p10_p00 = make_exp_ci(log_RR_10_00, var_log_RR_10_00)
  OR_p01_p00 = make_exp_ci(log_OR_01_00, var_log_OR_01_00)
  OR_p11_p10 = make_exp_ci(log_OR_11_10, var_log_OR_11_10)
  
  measure <- list (
      mea.point =c(RR_p10_p00$log_est,OR_p01_p00$log_est,OR_p11_p10$log_est),
      mea.se = c(RR_p10_p00$se_log,OR_p01_p00$se_log,OR_p11_p10$se_log),
      mea.CI.low = c(RR_p10_p00$lower,OR_p01_p00$lower,OR_p11_p10$lower),
      mea.CI.up = c(RR_p10_p00$upper,OR_p01_p00$upper,OR_p11_p10$upper),
      mea.p = c(RR_p10_p00$p,OR_p01_p00$p,OR_p11_p10$p)
  )
  cellprob <- pred$cell_prob
  list(measures = measure,
       cellprob = cellprob)
}