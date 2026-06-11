## ============================================================
## OP / gOP + three-measure MLE for binary outcome probabilities
##
## p_ab(x) = P(Y = 1 | A = a, B = b, X = x), a,b in {0,1}
##
## Model parameters:
##   nuisance: log OP or log gOP
##   measure10: p10 vs p00
##   measure01: p01 vs p00
##   measure11: p11 vs p10
##
## Supported measures:
##   OR  : logit(p1) - logit(p0)
##   RR  : log(p1 / p0)
##   SR  : log((1 - p1) / (1 - p0))
##   RD  : atanh(p1 - p0)
##   CMH : log{-log(1 - p1)} - log{-log(1 - p0)}
##         i.e. cumulative-hazard / cloglog contrast
##
## Supported nuisance:
##   nuisance = "gOP": logit(p00)+logit(p10)+logit(p01)+logit(p11)
##   nuisance = "OP" or "OP10": logit(p00)+logit(p10)
##
## Data can be individual-level or grouped binomial data.
## Required columns:
##   y : binary outcome or number of events
##   n : optional denominator; omit for individual-level binary data
##   A : binary 0/1
##   B : binary 0/1
## ============================================================

opm_clip01 <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}

opm_clip_interval <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

opm_normalize_measure <- function(x) {
  x <- toupper(x)
  ok <- c("OR", "RR", "SR", "RD", "CMH")
  if (!all(x %in% ok)) {
    stop("measure must be one of: OR, RR, SR, RD, CMH")
  }
  x
}

opm_normalize_measures <- function(measures) {
  target_names <- c("p10_p00", "p01_p00", "p11_p10")

  if (length(measures) == 1) {
    measures <- rep(measures, 3)
  }

  if (length(measures) != 3) {
    stop("measures must have length 3, ordered as c('p10_p00', 'p01_p00', 'p11_p10').")
  }

  if (is.null(names(measures)) || any(names(measures) == "")) {
    names(measures) <- target_names
  } else {
    missing_names <- setdiff(target_names, names(measures))
    if (length(missing_names) > 0) {
      stop("Named measures must include: p10_p00, p01_p00, p11_p10")
    }
    measures <- measures[target_names]
  }

  nm <- names(measures)
  measures <- opm_normalize_measure(measures)
  names(measures) <- nm
  measures
}

opm_normalize_nuisance <- function(nuisance) {
  nuisance <- toupper(nuisance[1])
  if (nuisance %in% c("GOP", "GENERALIZED_OP", "GENERALIZED ODDS PRODUCT")) {
    return("GOP")
  }
  if (nuisance %in% c("OP", "OP10", "PAIR_OP", "PAIRWISE_OP", "ODDS PRODUCT")) {
    return("OP10")
  }
  stop("nuisance must be 'gOP' or 'OP'/'OP10'.")
}

opm_ch <- function(p, eps = 1e-10) {
  -log1p(-opm_clip01(p, eps))
}

opm_pair_prob <- function(pref, eta, measure, eps = 1e-10) {
  measure <- toupper(measure)

  if (!is.finite(pref) || !is.finite(eta)) {
    return(NA_real_)
  }

  if (measure == "OR") {
    return(plogis(qlogis(opm_clip01(pref, eps)) + eta))
  }

  if (measure == "RR") {
    return(pref * exp(eta))
  }

  if (measure == "SR") {
    return(1 - (1 - pref) * exp(eta))
  }

  if (measure == "RD") {
    return(pref + tanh(eta))
  }

  if (measure == "CMH") {
    return(1 - exp(-opm_ch(pref, eps) * exp(eta)))
  }

  stop("Unknown measure.")
}

## For a fixed eta and measure, return the interval of reference risks
## pref such that pair_prob(pref, eta, measure) is in (eps, 1 - eps).
opm_ref_interval <- function(eta, measure, eps = 1e-10) {
  measure <- toupper(measure)
  lower <- eps
  upper <- 1 - eps

  if (!is.finite(eta)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  if (measure == "OR") {
    lower <- max(lower, stats::plogis(stats::qlogis(eps) - eta))
    upper <- min(upper, stats::plogis(stats::qlogis(1 - eps) - eta))
  } else if (measure == "RR") {
    r <- exp(eta)
    lower <- max(lower, eps / r)
    upper <- min(upper, (1 - eps) / r)
  } else if (measure == "SR") {
    s <- exp(eta)
    lower <- max(lower, 1 - (1 - eps) / s)
    upper <- min(upper, 1 - eps / s)
  } else if (measure == "RD") {
    d <- tanh(eta)
    lower <- max(lower, eps - d)
    upper <- min(upper, 1 - eps - d)
  } else if (measure == "CMH") {
    h <- exp(eta)
    lower <- max(lower, 1 - exp(-opm_ch(eps, eps) / h))
    upper <- min(upper, 1 - exp(-opm_ch(1 - eps, eps) / h))
  } else {
    stop("Unknown measure.")
  }

  c(lower = lower, upper = upper)
}

## Inverse of pair_prob in the first argument.
opm_inv_pair_prob <- function(p1, eta, measure, eps = 1e-10) {
  measure <- toupper(measure)
  p1 <- opm_clip01(p1, eps)

  if (!is.finite(eta)) {
    return(NA_real_)
  }

  if (measure == "OR") {
    return(stats::plogis(stats::qlogis(p1) - eta))
  }

  if (measure == "RR") {
    return(p1 / exp(eta))
  }

  if (measure == "SR") {
    return(1 - (1 - p1) / exp(eta))
  }

  if (measure == "RD") {
    return(p1 - tanh(eta))
  }

  if (measure == "CMH") {
    return(1 - exp(-opm_ch(p1, eps) / exp(eta)))
  }

  stop("Unknown measure.")
}

## Analytic valid interval for p00, including constraints from p10, p01, and p11.
opm_p00_domain <- function(eta10, eta01, eta11, measures, eps = 1e-10) {
  lower <- eps
  upper <- 1 - eps

  int10 <- opm_ref_interval(eta10, measures["p10_p00"], eps = eps)
  int01 <- opm_ref_interval(eta01, measures["p01_p00"], eps = eps)
  int11_ref_p10 <- opm_ref_interval(eta11, measures["p11_p10"], eps = eps)

  lower <- max(lower, int10["lower"], int01["lower"], na.rm = TRUE)
  upper <- min(upper, int10["upper"], int01["upper"], na.rm = TRUE)

  ## p11 validity imposes that p10 must be inside int11_ref_p10.
  ## Since p10 = pair_prob(p00, eta10, measure10) is monotone increasing,
  ## we map this p10 interval back to p00.
  back_lower <- opm_inv_pair_prob(int11_ref_p10["lower"], eta10, measures["p10_p00"], eps = eps)
  back_upper <- opm_inv_pair_prob(int11_ref_p10["upper"], eta10, measures["p10_p00"], eps = eps)
  back <- sort(c(back_lower, back_upper))

  lower <- max(lower, back[1], na.rm = TRUE)
  upper <- min(upper, back[2], na.rm = TRUE)

  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  c(lower = lower, upper = upper)
}

opm_measure_eta <- function(p1, p0, measure, eps = 1e-10) {
  measure <- toupper(measure)
  p1 <- opm_clip01(p1, eps)
  p0 <- opm_clip01(p0, eps)

  if (measure == "OR") {
    return(qlogis(p1) - qlogis(p0))
  }

  if (measure == "RR") {
    return(log(p1 / p0))
  }

  if (measure == "SR") {
    return(log((1 - p1) / (1 - p0)))
  }

  if (measure == "RD") {
    rd <- opm_clip_interval(p1 - p0, -1 + eps, 1 - eps)
    return(atanh(rd))
  }

  if (measure == "CMH") {
    return(log(opm_ch(p1, eps) / opm_ch(p0, eps)))
  }

  stop("Unknown measure.")
}

opm_nuisance_eta <- function(p, nuisance = "GOP", eps = 1e-10) {
  nuisance <- opm_normalize_nuisance(nuisance)
  p <- opm_clip01(p, eps)

  if (nuisance == "GOP") {
    return(sum(qlogis(p[c("p00", "p10", "p01", "p11")])) )
  }

  if (nuisance == "OP10") {
    return(qlogis(p["p00"]) + qlogis(p["p10"]))
  }

  stop("Unknown nuisance.")
}

opm_eval_p_from_p00 <- function(p00, eta10, eta01, eta11, measures, eps = 1e-10) {
  p10 <- opm_pair_prob(p00, eta10, measures["p10_p00"], eps = eps)
  p01 <- opm_pair_prob(p00, eta01, measures["p01_p00"], eps = eps)
  p11 <- opm_pair_prob(p10, eta11, measures["p11_p10"], eps = eps)

  c(p00 = p00, p10 = p10, p01 = p01, p11 = p11)
}

opm_valid_p <- function(p, eps = 1e-10) {
  all(is.finite(p)) && all(p > eps) && all(p < 1 - eps)
}

opm_prob_one <- function(eta_nuisance,
                         eta10,
                         eta01,
                         eta11,
                         measures,
                         nuisance = "GOP",
                         eps = 1e-10,
                         root_grid = 81) {

  nuisance <- opm_normalize_nuisance(nuisance)

  hfun <- function(p00) {
    p <- opm_eval_p_from_p00(
      p00 = p00,
      eta10 = eta10,
      eta01 = eta01,
      eta11 = eta11,
      measures = measures,
      eps = eps
    )

    if (!opm_valid_p(p, eps = eps)) {
      return(NA_real_)
    }

    opm_nuisance_eta(p, nuisance = nuisance, eps = eps) - eta_nuisance
  }

  valid_fun <- function(p00) {
    p <- opm_eval_p_from_p00(
      p00 = p00,
      eta10 = eta10,
      eta01 = eta01,
      eta11 = eta11,
      measures = measures,
      eps = eps
    )
    opm_valid_p(p, eps = eps)
  }

  dom <- opm_p00_domain(eta10, eta01, eta11, measures, eps = eps)
  lower <- dom["lower"]
  upper <- dom["upper"]

  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    return(c(p00 = NA_real_, p10 = NA_real_, p01 = NA_real_, p11 = NA_real_))
  }

  width <- upper - lower
  lower <- lower + 1e-8 * width
  upper <- upper - 1e-8 * width

  if (!valid_fun(lower) || !valid_fun(upper)) {
    lower <- lower + 1e-6 * width
    upper <- upper - 1e-6 * width
  }

  h_lower <- hfun(lower)
  h_upper <- hfun(upper)

  ## In the gOP case, the analytic domain usually brackets every real nuisance eta.
  ## OP10 can fail when p01/p11 constraints make the implied probabilities invalid;
  ## then we use a grid inside the analytic domain to check for a crossing.
  if (!is.finite(h_lower) || !is.finite(h_upper) || h_lower * h_upper > 0) {
    grid <- seq(lower, upper, length.out = root_grid)
    h_grid <- vapply(grid, hfun, numeric(1))
    ok <- is.finite(h_grid)

    if (!any(ok)) {
      return(c(p00 = NA_real_, p10 = NA_real_, p01 = NA_real_, p11 = NA_real_))
    }

    if (any(abs(h_grid[ok]) < 1e-10)) {
      root <- grid[which(ok)[which.min(abs(h_grid[ok]))]]
      return(opm_eval_p_from_p00(root, eta10, eta01, eta11, measures, eps = eps))
    }

    sg <- sign(h_grid)
    cross <- which(ok[-length(ok)] & ok[-1] & sg[-length(sg)] * sg[-1] <= 0)
    if (length(cross) == 0) {
      return(c(p00 = NA_real_, p10 = NA_real_, p01 = NA_real_, p11 = NA_real_))
    }
    lower <- grid[cross[1]]
    upper <- grid[cross[1] + 1]
  }

  root <- tryCatch(
    stats::uniroot(hfun, interval = c(lower, upper), tol = eps)$root,
    error = function(e) NA_real_
  )

  if (!is.finite(root)) {
    return(c(p00 = NA_real_, p10 = NA_real_, p01 = NA_real_, p11 = NA_real_))
  }

  opm_eval_p_from_p00(root, eta10, eta01, eta11, measures, eps = eps)
}

opm_build_X <- function(formulas, data) {
  list(
    Xn = stats::model.matrix(formulas$f_nuisance, data),
    X10 = stats::model.matrix(formulas$f10, data),
    X01 = stats::model.matrix(formulas$f01, data),
    X11 = stats::model.matrix(formulas$f11, data)
  )
}

opm_unpack <- function(par, lens) {
  i1 <- lens[1]
  i2 <- i1 + lens[2]
  i3 <- i2 + lens[3]
  i4 <- i3 + lens[4]

  list(
    alpha_nuisance = par[seq_len(i1)],
    beta10 = par[(i1 + 1):i2],
    beta01 = par[(i2 + 1):i3],
    beta11 = par[(i3 + 1):i4]
  )
}

opm_prob_matrix <- function(par,
                            Xs,
                            lens,
                            measures,
                            nuisance = "GOP",
                            eps = 1e-10,
                            root_grid = 81) {
  b <- opm_unpack(par, lens)

  eta_n <- drop(Xs$Xn %*% b$alpha_nuisance)
  eta10 <- drop(Xs$X10 %*% b$beta10)
  eta01 <- drop(Xs$X01 %*% b$beta01)
  eta11 <- drop(Xs$X11 %*% b$beta11)

  n_obs <- length(eta_n)
  P <- matrix(NA_real_, nrow = n_obs, ncol = 4)
  colnames(P) <- c("p00", "p10", "p01", "p11")

  for (i in seq_len(n_obs)) {
    P[i, ] <- opm_prob_one(
      eta_nuisance = eta_n[i],
      eta10 = eta10[i],
      eta01 = eta01[i],
      eta11 = eta11[i],
      measures = measures,
      nuisance = nuisance,
      eps = eps,
      root_grid = root_grid
    )
  }

  P
}

opm_safe_vcov <- function(H) {
  H <- (H + t(H)) / 2

  V <- tryCatch(
    solve(H),
    error = function(e) NULL
  )

  if (is.null(V) || any(!is.finite(V))) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Hessian is singular. Install MASS or provide a better starting value.")
    }
    V <- MASS::ginv(H)
  }

  (V + t(V)) / 2
}

opm_wald_table <- function(est, se, conf.level = 0.95) {
  zcrit <- stats::qnorm(1 - (1 - conf.level) / 2)
  zval <- est / se
  pval <- 2 * stats::pnorm(abs(zval), lower.tail = FALSE)

  data.frame(
    term = names(est),
    estimate = as.numeric(est),
    se = as.numeric(se),
    z = as.numeric(zval),
    p.value = as.numeric(pval),
    conf.low = as.numeric(est - zcrit * se),
    conf.high = as.numeric(est + zcrit * se),
    row.names = NULL
  )
}

opm_transform_table <- function(est, se, type, conf.level = 0.95) {
  zcrit <- stats::qnorm(1 - (1 - conf.level) / 2)

  if (type == "RD") {
    trans_est <- tanh(est)
    trans_se <- (1 - trans_est^2) * se
    low <- tanh(est - zcrit * se)
    high <- tanh(est + zcrit * se)
  } else {
    trans_est <- exp(est)
    trans_se <- exp(est) * se
    low <- exp(est - zcrit * se)
    high <- exp(est + zcrit * se)
  }

  data.frame(
    term = names(est),
    scale = type,
    estimate = as.numeric(trans_est),
    se.delta = as.numeric(trans_se),
    conf.low = as.numeric(low),
    conf.high = as.numeric(high),
    row.names = NULL
  )
}

opm_fd_jacobian <- function(fun, x, rel_step = 1e-5) {
  f0 <- fun(x)
  J <- matrix(NA_real_, nrow = length(f0), ncol = length(x))

  for (j in seq_along(x)) {
    h <- rel_step * (abs(x[j]) + 1)
    xp <- x
    xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h

    fp <- fun(xp)
    fm <- fun(xm)
    J[, j] <- (fp - fm) / (2 * h)
  }

  rownames(J) <- names(f0)
  colnames(J) <- names(x)
  J
}

opm_initial_start <- function(yv, nv, Av, Bv, measures, nuisance, Xs, lens, par_names, eps = 1e-10) {
  cell_p <- function(a, b) {
    idx <- Av == a & Bv == b
    if (!any(idx)) {
      return((sum(yv) + 0.5) / (sum(nv) + 1))
    }
    (sum(yv[idx]) + 0.5) / (sum(nv[idx]) + 1)
  }

  p_hat <- c(
    p00 = cell_p(0, 0),
    p10 = cell_p(1, 0),
    p01 = cell_p(0, 1),
    p11 = cell_p(1, 1)
  )

  eta_n_hat <- opm_nuisance_eta(p_hat, nuisance = nuisance, eps = eps)
  eta10_hat <- opm_measure_eta(p_hat["p10"], p_hat["p00"], measures["p10_p00"], eps = eps)
  eta01_hat <- opm_measure_eta(p_hat["p01"], p_hat["p00"], measures["p01_p00"], eps = eps)
  eta11_hat <- opm_measure_eta(p_hat["p11"], p_hat["p10"], measures["p11_p10"], eps = eps)

  start <- rep(0, sum(lens))
  names(start) <- par_names

  offset_n <- 0
  offset10 <- lens[1]
  offset01 <- lens[1] + lens[2]
  offset11 <- lens[1] + lens[2] + lens[3]

  if ("(Intercept)" %in% colnames(Xs$Xn)) {
    start[offset_n + match("(Intercept)", colnames(Xs$Xn))] <- eta_n_hat
  }
  if ("(Intercept)" %in% colnames(Xs$X10)) {
    start[offset10 + match("(Intercept)", colnames(Xs$X10))] <- eta10_hat
  }
  if ("(Intercept)" %in% colnames(Xs$X01)) {
    start[offset01 + match("(Intercept)", colnames(Xs$X01))] <- eta01_hat
  }
  if ("(Intercept)" %in% colnames(Xs$X11)) {
    start[offset11 + match("(Intercept)", colnames(Xs$X11))] <- eta11_hat
  }

  start
}

fit_op_measure_mle <- function(data,
                               y,
                               A,
                               B,
                               n = NULL,
                               f_nuisance = ~ 1,
                               f10 = ~ 1,
                               f01 = ~ 1,
                               f11 = ~ 1,
                               measures = c("OR", "RR", "RR"),
                               nuisance = "gOP",
                               start = NULL,
                               eps = 1e-10,
                               root_grid = 81,
                               conf.level = 0.95,
                               control_nm = list(maxit = 3000, reltol = 1e-8),
                               control_bfgs = list(maxit = 3000, reltol = 1e-9)) {

  measures <- opm_normalize_measures(measures)
  nuisance <- opm_normalize_nuisance(nuisance)

  yv <- data[[y]]
  nv <- if (is.null(n)) rep(1, length(yv)) else data[[n]]
  Av <- as.integer(data[[A]])
  Bv <- as.integer(data[[B]])

  if (!all(Av %in% c(0, 1))) stop("A must be coded as 0/1.")
  if (!all(Bv %in% c(0, 1))) stop("B must be coded as 0/1.")
  if (any(yv < 0 | yv > nv)) stop("y must satisfy 0 <= y <= n.")

  formulas <- list(
    f_nuisance = f_nuisance,
    f10 = f10,
    f01 = f01,
    f11 = f11
  )

  Xs <- opm_build_X(formulas, data)

  lens <- c(
    nuisance = ncol(Xs$Xn),
    p10_p00 = ncol(Xs$X10),
    p01_p00 = ncol(Xs$X01),
    p11_p10 = ncol(Xs$X11)
  )

  nuisance_label <- if (nuisance == "GOP") "log_gOP" else "log_OP10"

  par_names <- c(
    paste0(nuisance_label, ":", colnames(Xs$Xn)),
    paste0("eta_", measures["p10_p00"], "_p10_vs_p00:", colnames(Xs$X10)),
    paste0("eta_", measures["p01_p00"], "_p01_vs_p00:", colnames(Xs$X01)),
    paste0("eta_", measures["p11_p10"], "_p11_vs_p10:", colnames(Xs$X11))
  )

  if (is.null(start)) {
    start <- opm_initial_start(
      yv = yv,
      nv = nv,
      Av = Av,
      Bv = Bv,
      measures = measures,
      nuisance = nuisance,
      Xs = Xs,
      lens = lens,
      par_names = par_names,
      eps = eps
    )
  } else {
    if (length(start) != sum(lens)) {
      stop("Length of start is not equal to the number of model parameters.")
    }
    names(start) <- par_names
  }

  observed_prob_from_par <- function(par) {
    P <- opm_prob_matrix(
      par = par,
      Xs = Xs,
      lens = lens,
      measures = measures,
      nuisance = nuisance,
      eps = eps,
      root_grid = root_grid
    )

    observed.prob(P, Av, Bv)
  }

  negloglik <- function(par) {
    if (any(!is.finite(par))) {
      return(1e80)
    }

    pobs <- observed_prob_from_par(par)

    if (any(!is.finite(pobs)) || any(pobs <= eps) || any(pobs >= 1 - eps)) {
      return(1e80 + sum(par^2, na.rm = TRUE))
    }

    -sum(stats::dbinom(yv, size = nv, prob = pobs, log = TRUE))
  }

  fit0 <- tryCatch(
    stats::optim(
      par = start,
      fn = negloglik,
      method = "Nelder-Mead",
      control = control_nm
    ),
    error = function(e) NULL
  )

  start_bfgs <- if (is.null(fit0)) start else fit0$par

  fit1 <- tryCatch(
    stats::optim(
      par = start_bfgs,
      fn = negloglik,
      method = "BFGS",
      hessian = TRUE,
      control = control_bfgs
    ),
    error = function(e) NULL
  )

  if (is.null(fit1) && is.null(fit0)) {
    stop("optim failed from the supplied start values.")
  }

  if (is.null(fit1)) {
    fit <- fit0
    fit$hessian <- stats::optimHess(fit$par, negloglik)
  } else if (!is.null(fit0) && is.finite(fit0$value) && fit0$value < fit1$value) {
    fit <- fit0
    fit$hessian <- stats::optimHess(fit$par, negloglik)
  } else {
    fit <- fit1
  }

  names(fit$par) <- par_names

  V <- opm_safe_vcov(fit$hessian)
  rownames(V) <- colnames(V) <- par_names

  se <- sqrt(pmax(diag(V), 0))
  names(se) <- par_names

  coef_table <- opm_wald_table(fit$par, se, conf.level = conf.level)

  idx_n <- seq_len(lens[1])
  idx10 <- lens[1] + seq_len(lens[2])
  idx01 <- lens[1] + lens[2] + seq_len(lens[3])
  idx11 <- lens[1] + lens[2] + lens[3] + seq_len(lens[4])

  nuisance_table <- opm_transform_table(
    est = fit$par[idx_n],
    se = se[idx_n],
    type = "OP/GOP",
    conf.level = conf.level
  )

  measure_table <- rbind(
    opm_transform_table(fit$par[idx10], se[idx10], measures["p10_p00"], conf.level),
    opm_transform_table(fit$par[idx01], se[idx01], measures["p01_p00"], conf.level),
    opm_transform_table(fit$par[idx11], se[idx11], measures["p11_p10"], conf.level)
  )

  fit_info <- data.frame(
    logLik = -fit$value,
    AIC = 2 * length(fit$par) + 2 * fit$value,
    convergence = fit$convergence,
    row.names = NULL
  )

  out <- list(
    call = match.call(),
    nuisance = nuisance,
    measures = measures,
    coefficients = fit$par,
    vcov = V,
    coef_table = coef_table,
    nuisance_table = nuisance_table,
    measure_table = measure_table,
    fit_info = fit_info,
    convergence = fit$convergence,
    message = fit$message,
    formulas = formulas,
    lens = lens,
    x_colnames = list(
      Xn = colnames(Xs$Xn),
      X10 = colnames(Xs$X10),
      X01 = colnames(Xs$X01),
      X11 = colnames(Xs$X11)
    ),
    data = data,
    y = y,
    n = n,
    A = A,
    B = B,
    eps = eps,
    root_grid = root_grid,
    negloglik = negloglik
  )

  class(out) <- "op_measure_mle"
  out
}

predict.op_measure_mle <- function(object,
                                   newdata = NULL,
                                   se.fit = FALSE,
                                   conf.level = 0.95,
                                   rel_step = 1e-5,
                                   ...) {
  if (is.null(newdata)) {
    newdata <- object$data
  }

  Xs <- opm_build_X(object$formulas, newdata)

  Xs$Xn <- Xs$Xn[, object$x_colnames$Xn, drop = FALSE]
  Xs$X10 <- Xs$X10[, object$x_colnames$X10, drop = FALSE]
  Xs$X01 <- Xs$X01[, object$x_colnames$X01, drop = FALSE]
  Xs$X11 <- Xs$X11[, object$x_colnames$X11, drop = FALSE]

  prob_fun <- function(par) {
    P <- opm_prob_matrix(
      par = par,
      Xs = Xs,
      lens = object$lens,
      measures = object$measures,
      nuisance = object$nuisance,
      eps = object$eps,
      root_grid = object$root_grid
    )
    v <- as.vector(t(P))
    names(v) <- as.vector(t(outer(seq_len(nrow(P)), colnames(P), paste, sep = ":")))
    v
  }

  P <- opm_prob_matrix(
    par = object$coefficients,
    Xs = Xs,
    lens = object$lens,
    measures = object$measures,
    nuisance = object$nuisance,
    eps = object$eps,
    root_grid = object$root_grid
  )

  ans <- data.frame(
    row = seq_len(nrow(P)),
    p00 = P[, "p00"],
    p10 = P[, "p10"],
    p01 = P[, "p01"],
    p11 = P[, "p11"],
    row.names = NULL
  )

  if (se.fit) {
    J <- opm_fd_jacobian(prob_fun, object$coefficients, rel_step = rel_step)
    VV <- J %*% object$vcov %*% t(J)
    se_vec <- sqrt(pmax(diag(VV), 0))
    se_mat <- matrix(se_vec, nrow = nrow(P), ncol = 4, byrow = TRUE)
    colnames(se_mat) <- paste0(colnames(P), ".se")

    zcrit <- stats::qnorm(1 - (1 - conf.level) / 2)
    low <- pmax(0, P - zcrit * se_mat)
    high <- pmin(1, P + zcrit * se_mat)
    colnames(low) <- paste0(colnames(P), ".low")
    colnames(high) <- paste0(colnames(P), ".high")

    ans <- cbind(ans, as.data.frame(se_mat), as.data.frame(low), as.data.frame(high))
  }

  ans
}

print.op_measure_mle <- function(x, ...) {
  cat("\nNuisance parameterization:\n")
  cat("  ", x$nuisance, "\n", sep = "")
  cat("\nMeasures, in order p10/p00, p01/p00, p11/p10:\n")
  print(x$measures)
  cat("\nFit information:\n")
  print(x$fit_info)
  cat("\nCoefficient table, model scale:\n")
  print(x$coef_table)
  invisible(x)
}

## Convenience wrapper: p10,p00 OR; p01,p00 RR; p11,p10 RR.
fit_OP_OR_RR_RR <- function(data,
                            y,
                            A,
                            B,
                            n = NULL,
                            nuisance = "gOP",
                            f_nuisance = ~ 1,
                            f10 = ~ 1,
                            f01 = ~ 1,
                            f11 = ~ 1,
                            ...) {
  fit_op_measure_mle(
    data = data,
    y = y,
    A = A,
    B = B,
    n = n,
    f_nuisance = f_nuisance,
    f10 = f10,
    f01 = f01,
    f11 = f11,
    measures = c(p10_p00 = "OR", p01_p00 = "RR", p11_p10 = "RR"),
    nuisance = nuisance,
    ...
  )
}

## Convenience wrapper: p10,p00 SR; p01,p00 OR; p11,p10 OR.
fit_OP_SR_OR_OR <- function(data,
                            y,
                            A,
                            B,
                            n = NULL,
                            nuisance = "gOP",
                            f_nuisance = ~ 1,
                            f10 = ~ 1,
                            f01 = ~ 1,
                            f11 = ~ 1,
                            ...) {
  fit_op_measure_mle(
    data = data,
    y = y,
    A = A,
    B = B,
    n = n,
    f_nuisance = f_nuisance,
    f10 = f10,
    f01 = f01,
    f11 = f11,
    measures = c(p10_p00 = "SR", p01_p00 = "OR", p11_p10 = "OR"),
    nuisance = nuisance,
    ...
  )
}

## Convenience wrapper: p10,p00 OR; p01,p00 SR; p11,p10 SR.
fit_OP_OR_SR_SR <- function(data,
                            y,
                            A,
                            B,
                            n = NULL,
                            nuisance = "gOP",
                            f_nuisance = ~ 1,
                            f10 = ~ 1,
                            f01 = ~ 1,
                            f11 = ~ 1,
                            ...) {
  fit_op_measure_mle(
    data = data,
    y = y,
    A = A,
    B = B,
    n = n,
    f_nuisance = f_nuisance,
    f10 = f10,
    f01 = f01,
    f11 = f11,
    measures = c(p10_p00 = "OR", p01_p00 = "SR", p11_p10 = "SR"),
    nuisance = nuisance,
    ...
  )
}
