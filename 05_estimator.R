### MlEst-style wrapper and result object

WrapResults.TwoByMeasures <- function(point.est, cov, x0, x1, x2, x3, a0, a1,
                                      converged, fit = NULL, ci.est = NULL,
                                      measures = c("RR", "OR", "OR")) {
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  x3 <- as.matrix(x3)
  measures <- normalize.measures.Two(measures)
  measure.key <- measure.key.Two(measures)
  
  name <- coef.names(x0, x1, x2, x3)
  se.est <- sqrt(diag(cov))
  
  if (is.null(ci.est)) {
    conf.lower <- point.est + stats::qnorm(0.025) * se.est
    conf.upper <- point.est + stats::qnorm(0.975) * se.est
    p.temp <- stats::pnorm(point.est / se.est, 0, 1)
    p.value <- 2 * pmin(p.temp, 1 - p.temp)
  } else {
    conf.lower <- ci.est$low
    conf.upper <- ci.est$up
    p.value <- ci.est$p
  }
  
  names(point.est) <- names(se.est) <- names(conf.lower) <-
    names(conf.upper) <- names(p.value) <- name
  rownames(cov) <- colnames(cov) <- name
  
  q0 <- ncol(x0)
  q1 <- ncol(x1)
  q2 <- ncol(x2)
  q3 <- ncol(x3)
  b0 <- point.est[seq_len(q0)]
  b1 <- point.est[q0 + seq_len(q1)]
  b2 <- point.est[q0 + q1 + seq_len(q2)]
  b3 <- point.est[q0 + q1 + q2 + seq_len(q3)]
  
  theta0 <- as.vector(x0 %*% b0)
  theta1 <- as.vector(x1 %*% b1)
  theta2 <- as.vector(x2 %*% b2)
  phi <- as.vector(x3 %*% b3)
  prob.est <- getTwo.by.measures(
    theta0, theta1, theta2, phi,
    measures = measures
  )
  colnames(prob.est) <- c("p00", "p10", "p01", "p11")
  
  param.est <- cbind(
    exp(theta0),
    exp(theta1),
    exp(theta2),
    gOP = exp(phi)
  )
  colnames(param.est) <- c(
    paste0(measures[["p10_p00"]], "_p10_p00"),
    paste0(measures[["p01_p00"]], "_p01_p00"),
    paste0(measures[["p11_p10"]], "_p11_p10"),
    "gOP"
  )
  
  coefficients <- cbind(point.est, se.est, conf.lower, conf.upper, p.value)
  
  sol <- list(
    param = paste0("Two", measure.key),
    measures = measures,
    point.est = point.est,
    se.est = se.est,
    cov = cov,
    conf.lower = conf.lower,
    conf.upper = conf.upper,
    p.value = p.value,
    coefficients = coefficients,
    param.est = param.est,
    prob.est = prob.est,
    pA = observed.prob(prob.est, a0, a1),
    x0 = x0,
    x1 = x1,
    x2 = x2,
    x3 = x3,
    a0 = a0,
    a1 = a1,
    convergence = converged,
    fit = fit
  )
  class(sol) <- c("gop_two_mle", "list")
  attr(sol, "hidden") <- c(
    "param", "measures", "se.est", "cov", "conf.lower", "conf.upper",
    "p.value", "param.est", "prob.est", "pA", "x0", "x1", "x2",
    "x3", "a0", "a1", "converged", "convergence", "fit"
  )
  sol
}

WrapResults.TwoRROROR <- function(...) {
  WrapResults.TwoByMeasures(..., measures = c("RR", "OR", "OR"))
}

WrapResults.TwoORRRRR <- function(...) {
  WrapResults.TwoByMeasures(..., measures = c("OR", "RR", "RR"))
}

WrapResults.TwoSROROR <- function(...) {
  WrapResults.TwoByMeasures(..., measures = c("SR", "OR", "OR"))
}

WrapResults.TwoORSRSR <- function(...) {
  WrapResults.TwoByMeasures(..., measures = c("OR", "SR", "SR"))
}

MLEst.TwoByMeasures <- function(y = NULL, a0 = NULL, a1 = NULL,
                                x0 = NULL, x1 = NULL, x2 = NULL, x3 = NULL,
                                data = NULL,
                                measures = c("RR", "OR", "OR"),
                                weight = NULL,
                                max.step = 1000, thres = 1e-8,
                                b0.start = NULL, b1.start = NULL,
                                b2.start = NULL, b3.start = NULL,
                                CI = "wald", conf.level = 0.95,
                                na.omit = TRUE) {
  CI <- tolower(CI)
  measures <- normalize.measures.Two(measures)
  if (!identical(CI, "wald")) {
    stop("Only Wald confidence intervals are currently implemented.", call. = FALSE)
  }
  
  get.vector <- function(z, data) {
    if (!is.null(data) && is.character(z) && length(z) == 1) {
      return(data[[z]])
    }
    if (!is.null(data) && inherits(z, "formula")) {
      mf <- stats::model.frame(z, data = data, na.action = stats::na.pass)
      return(stats::model.response(mf))
    }
    z
  }
  
  get.matrix <- function(z, data) {
    if (is.null(z)) {
      stop("x0, x1, x2, and x3 must all be supplied.", call. = FALSE)
    }
    if (!is.null(data) && inherits(z, "formula")) {
      mf <- stats::model.frame(z, data = data, na.action = stats::na.pass)
      return(stats::model.matrix(z, data = mf))
    }
    if (!is.null(data) && is.character(z)) {
      return(as.matrix(data[, z, drop = FALSE]))
    }
    as.matrix(z)
  }
  
  y <- as.vector(get.vector(y, data))
  a0 <- as.vector(get.vector(a0, data))
  a1 <- as.vector(get.vector(a1, data))
  if (is.null(y) || is.null(a0) || is.null(a1)) {
    stop("y, a0, and a1 must be supplied as vectors or data column names.", call. = FALSE)
  }
  x0 <- get.matrix(x0, data)
  x1 <- get.matrix(x1, data)
  x2 <- get.matrix(x2, data)
  x3 <- get.matrix(x3, data)
  if (is.null(weight)) {
    weight <- rep(1, length(y))
  } else {
    weight <- as.vector(get.vector(weight, data))
  }
  
  n <- length(y)
  if (length(a0) != n || length(a1) != n || length(weight) != n ||
      nrow(x0) != n || nrow(x1) != n || nrow(x2) != n || nrow(x3) != n) {
    stop("y, a0, a1, weight, x0, x1, x2, and x3 must have compatible lengths.", call. = FALSE)
  }
  
  keep <- stats::complete.cases(y, a0, a1, weight) &
    stats::complete.cases(x0) &
    stats::complete.cases(x1) &
    stats::complete.cases(x2) &
    stats::complete.cases(x3)
  if (any(!keep)) {
    if (!na.omit) {
      stop("Missing values detected. Set na.omit = TRUE to drop incomplete rows.", call. = FALSE)
    }
    y <- y[keep]
    a0 <- a0[keep]
    a1 <- a1[keep]
    weight <- weight[keep]
    x0 <- x0[keep, , drop = FALSE]
    x1 <- x1[keep, , drop = FALSE]
    x2 <- x2[keep, , drop = FALSE]
    x3 <- x3[keep, , drop = FALSE]
  }
  
  if (is.null(b0.start)) {
    b0.start <- rep(0, ncol(x0))
  }
  if (is.null(b1.start)) {
    b1.start <- rep(0, ncol(x1))
  }
  if (is.null(b2.start)) {
    b2.start <- rep(0, ncol(x2))
  }
  if (is.null(b3.start)) {
    b3.start <- rep(0, ncol(x3))
    intercept <- which(apply(x3, 2, function(z) all(z == 1)))
    if (length(intercept) > 0) {
      pbar <- clip_prob(weighted.mean(y, weight))
      b3.start[intercept[1]] <- 4 * logit(pbar)
    }
  }
  
  mle <- maxlike.TwoByMeasures(
    y = y, a0 = a0, a1 = a1,
    x0 = x0, x1 = x1, x2 = x2, x3 = x3,
    b0.start = b0.start, b1.start = b1.start,
    b2.start = b2.start, b3.start = b3.start,
    weight = weight, measures = measures,
    max.step = max.step, thres = thres
  )
  
  point.est <- mle$par
  q0 <- ncol(x0)
  q1 <- ncol(x1)
  q2 <- ncol(x2)
  q3 <- ncol(x3)
  b0.ml <- point.est[seq_len(q0)]
  b1.ml <- point.est[q0 + seq_len(q1)]
  b2.ml <- point.est[q0 + q1 + seq_len(q2)]
  b3.ml <- point.est[q0 + q1 + q2 + seq_len(q3)]
  
  cov <- var.mle.by.measures(
    a0 = a0, a1 = a1,
    b0.ml = b0.ml, b1.ml = b1.ml,
    b2.ml = b2.ml, b3.ml = b3.ml,
    x0 = x0, x1 = x1, x2 = x2, x3 = x3,
    weight = weight,
    measures = measures
  )
  se.est <- sqrt(diag(cov))
  alpha <- 1 - conf.level
  ci.est <- list(
    low = point.est + stats::qnorm(alpha / 2) * se.est,
    up = point.est + stats::qnorm(1 - alpha / 2) * se.est,
    p = {
      p.temp <- stats::pnorm(point.est / se.est, 0, 1)
      2 * pmin(p.temp, 1 - p.temp)
    }
  )
  
  sol <- WrapResults.TwoByMeasures(
    point.est = point.est, cov = cov,
    x0 = x0, x1 = x1, x2 = x2, x3 = x3,
    a0 = a0, a1 = a1,
    converged = mle$convergence,
    fit = mle,
    ci.est = ci.est,
    measures = measures
  )
  sol$weight <- weight
  sol$n <- length(y)
  sol$na.action <- which(!keep)
  sol$conf.level <- conf.level
  sol
}

MLEst.TwoRROROR <- function(...) {
  MLEst.TwoByMeasures(..., measures = c("RR", "OR", "OR"))
}

MLEst.TwoORRRRR <- function(...) {
  MLEst.TwoByMeasures(..., measures = c("OR", "RR", "RR"))
}

MLEst.TwoSROROR <- function(...) {
  MLEst.TwoByMeasures(..., measures = c("SR", "OR", "OR"))
}

MLEst.TwoORSRSR <- function(...) {
  MLEst.TwoByMeasures(..., measures = c("OR", "SR", "SR"))
}

MlEst.TwoRROROR <- MLEst.TwoRROROR
MlEst.TwoORRRRR <- MLEst.TwoORRRRR
MlEst.TwoSROROR <- MLEst.TwoSROROR
MlEst.TwoORSRSR <- MLEst.TwoORSRSR
MlEst.TwoByMeasures <- MLEst.TwoByMeasures

print.gop_two_mle <- function(x, ...) {
  measures <- if (!is.null(x$measures)) {
    paste(unname(x$measures), collapse = "/")
  } else {
    "RR/OR/OR"
  }
  cat("Two-time-point gOP MLE (", measures, ")\n", sep = "")
  cat("  n:", x$n, "\n")
  cat("  converged:", x$converged, "\n")
  if (!is.null(x$fit$value)) {
    cat("  negative log-likelihood:", x$fit$value, "\n")
  }
  print(x$coefficients)
  invisible(x)
}
