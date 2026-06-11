### Fisher-information variance estimate

coef.names <- function(x0, x1, x2, x3) {
  make.names.block <- function(prefix, x) {
    x <- as.matrix(x)
    cn <- colnames(x)
    if (is.null(cn)) {
      cn <- seq_len(ncol(x))
    }
    paste0(prefix, ".", cn)
  }
  
  c(
    make.names.block("alpha1", x0),
    make.names.block("alpha2", x1),
    make.names.block("alpha3", x2),
    make.names.block("beta", x3)
  )
}

coef.names.TwoRROROR <- coef.names
coef.names.TwoORRRRR <- coef.names
coef.names.TwoSROROR <- coef.names
coef.names.TwoORSRSR <- coef.names

var.mle.TwoRROROR <- function(a0, a1, b0.ml, b1.ml, b2.ml, b3.ml,
                              x0, x1, x2, x3,
                              weight = rep(1, length(a0))) {
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  x3 <- as.matrix(x3)
  a0 <- as.vector(a0)
  a1 <- as.vector(a1)
  weight <- as.vector(weight)
  
  theta0 <- x0 %*% b0.ml
  theta1 <- x1 %*% b1.ml
  theta2 <- x2 %*% b2.ml
  phi <- x3 %*% b3.ml
  
  pa0a1 <- getTwoRROROR(theta0, theta1, theta2, phi)
  p00 <- pa0a1[, 1]
  p10 <- pa0a1[, 2]
  pA <- observed.prob(pa0a1, a0, a1)
  n <- nrow(x0)
  
  aa <- (1 - p00) / (1 - p10)
  bb <- 1 / (1 - p10)
  h <- 1 / (2 * (1 + aa))
  k <- bb / (1 + aa)
  a.tilde <- aa / (2 * (1 + aa))
  
  D00 <- cbind(-k * x0, -h * x1, -h * x2, h * x3)
  D10 <- cbind(k * x0, -a.tilde * x1, -a.tilde * x2, a.tilde * x3)
  D01 <- cbind(-k * x0, (1 - h) * x1, -h * x2, h * x3)
  D11 <- cbind(k * x0, -a.tilde * x1, (1 - a.tilde) * x2, a.tilde * x3)
  
  cell00 <- a0 == 0 & a1 == 0
  cell10 <- a0 == 1 & a1 == 0
  cell01 <- a0 == 0 & a1 == 1
  cell11 <- a0 == 1 & a1 == 1
  
  tmp <- matrix(NA_real_, n, ncol(D00))
  tmp[cell00, ] <- D00[cell00, , drop = FALSE]
  tmp[cell10, ] <- D10[cell10, , drop = FALSE]
  tmp[cell01, ] <- D01[cell01, , drop = FALSE]
  tmp[cell11, ] <- D11[cell11, , drop = FALSE]
  
  if (any(!is.finite(pA)) || any(!is.finite(tmp))) {
    stop("Non-finite fitted probabilities or derivatives in Fisher information.", call. = FALSE)
  }
  
  info.weight <- pA * (1 - pA)
  fisher.info <- t(info.weight * weight * tmp) %*% tmp
  var.mat <- if (requireNamespace("MASS", quietly = TRUE)) {
    MASS::ginv(fisher.info)
  } else {
    solve(fisher.info)
  }
  
  par.names <- coef.names(x0, x1, x2, x3)
  if (length(par.names) == ncol(var.mat) && all(!is.na(par.names))) {
    rownames(var.mat) <- par.names
    colnames(var.mat) <- par.names
  }
  var.mat
}

var.mle.TwoORRRRR <- function(a0, a1,
                              b0.ml, b1.ml, b2.ml, b3.ml,
                              x0, x1, x2, x3,
                              weight = rep(1, length(a0)),
                              eps = 1e-10) {
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  x3 <- as.matrix(x3)
  
  a0 <- as.vector(a0)
  a1 <- as.vector(a1)
  weight <- as.vector(weight)
  
  n <- length(a0)
  
  if (length(a1) != n || length(weight) != n) {
    stop("a0, a1, and weight must have the same length.", call. = FALSE)
  }
  
  if (nrow(x0) != n || nrow(x1) != n || nrow(x2) != n || nrow(x3) != n) {
    stop("x0, x1, x2, and x3 must have the same number of rows as length(a0).",
         call. = FALSE)
  }
  
  theta0 <- as.vector(x0 %*% b0.ml)
  theta1 <- as.vector(x1 %*% b1.ml)
  theta2 <- as.vector(x2 %*% b2.ml)
  phi    <- as.vector(x3 %*% b3.ml)
  
  pa0a1 <- getTwoORRRRR(theta0, theta1, theta2, phi, eps = eps)
  
  if (is.null(colnames(pa0a1))) {
    colnames(pa0a1) <- c("p00", "p10", "p01", "p11")
  }
  
  p00 <- pa0a1[, "p00"]
  p10 <- pa0a1[, "p10"]
  p01 <- pa0a1[, "p01"]
  p11 <- pa0a1[, "p11"]
  
  pA <- observed.prob(pa0a1, a0, a1)
  pA <- clip_prob(pA, eps = eps)
  
  if (any(!is.finite(pa0a1)) || any(!is.finite(pA))) {
    stop("Non-finite fitted probabilities in Fisher information.",
         call. = FALSE)
  }
  
  if (any(pa0a1 <= eps | pa0a1 >= 1 - eps, na.rm = TRUE)) {
    warning("Some fitted cell probabilities are close to 0 or 1; variance may be unstable.")
  }
  
  q00 <- 1 - p00
  q10 <- 1 - p10
  q01 <- 1 - p01
  q11 <- 1 - p11
  
  if (any(q01 <= eps | q11 <= eps, na.rm = TRUE)) {
    stop("q01 or q11 is too close to 0; derivatives are unstable.",
         call. = FALSE)
  }
  
  A <- q00 / q01
  B <- q10 / q11
  S <- 2 + A + B
  
  if (any(!is.finite(A)) || any(!is.finite(B)) || any(!is.finite(S)) ||
      any(S <= 0)) {
    stop("Non-finite or invalid derivative coefficients.",
         call. = FALSE)
  }
  
  D00 <- cbind(
    -(1 + B) / S * x0,
    -1 / (S * q01) * x1,
    -1 / (S * q11) * x2,
    1 / S * x3
  )
  
  D10 <- cbind(
    (1 + A) / S * x0,
    -1 / (S * q01) * x1,
    -1 / (S * q11) * x2,
    1 / S * x3
  )
  
  D01 <- cbind(
    -A * (1 + B) / S * x0,
    (S - A) / (S * q01) * x1,
    -A / (S * q11) * x2,
    A / S * x3
  )
  
  D11 <- cbind(
    B * (1 + A) / S * x0,
    -B / (S * q01) * x1,
    (S - B) / (S * q11) * x2,
    B / S * x3
  )
  
  cell00 <- a0 == 0 & a1 == 0
  cell10 <- a0 == 1 & a1 == 0
  cell01 <- a0 == 0 & a1 == 1
  cell11 <- a0 == 1 & a1 == 1
  
  tmp <- matrix(NA_real_, n, ncol(D00))
  
  tmp[cell00, ] <- D00[cell00, , drop = FALSE]
  tmp[cell10, ] <- D10[cell10, , drop = FALSE]
  tmp[cell01, ] <- D01[cell01, , drop = FALSE]
  tmp[cell11, ] <- D11[cell11, , drop = FALSE]
  
  if (any(!is.finite(tmp))) {
    stop("Non-finite derivatives in Fisher information.",
         call. = FALSE)
  }
  
  info.weight <- pA * (1 - pA)
  
  fisher.info <- t(info.weight * weight * tmp) %*% tmp
  
  var.mat <- if (requireNamespace("MASS", quietly = TRUE)) {
    MASS::ginv(fisher.info)
  } else {
    solve(fisher.info)
  }
  
  par.names <- coef.names(x0, x1, x2, x3)
  
  if (length(par.names) == ncol(var.mat) && all(!is.na(par.names))) {
    rownames(var.mat) <- par.names
    colnames(var.mat) <- par.names
  }
  
  var.mat
}

var.mle.TwoSROROR <- function(a0, a1, b0.ml, b1.ml, b2.ml, b3.ml,
                              x0, x1, x2, x3,
                              weight = rep(1, length(a0))) {
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  x3 <- as.matrix(x3)
  
  a0 <- as.vector(a0)
  a1 <- as.vector(a1)
  weight <- as.vector(weight)
  
  theta0 <- as.vector(x0 %*% b0.ml)
  theta1 <- as.vector(x1 %*% b1.ml)
  theta2 <- as.vector(x2 %*% b2.ml)
  phi    <- as.vector(x3 %*% b3.ml)
  
  pa0a1 <- getTwoSROROR(theta0, theta1, theta2, phi)
  
  p00 <- pa0a1[, 1]
  p10 <- pa0a1[, 2]
  
  pA <- observed.prob(pa0a1, a0, a1)
  
  n <- nrow(x0)
  
  ## SROROR 的核心系数
  den <- p00 + p10
  
  k <- 1 / den
  h <- p10 / (2 * den)
  a.tilde <- p00 / (2 * den)
  
  ## D00 = d logit(p00) / d beta
  D00 <- cbind(
    k * x0,
    -h * x1,
    -h * x2,
    h * x3
  )
  
  ## D10 = d logit(p10) / d beta
  D10 <- cbind(
    -k * x0,
    -a.tilde * x1,
    -a.tilde * x2,
    a.tilde * x3
  )
  
  ## D01 = d logit(p01) / d beta
  ## because logit(p01) = logit(p00) + theta1
  D01 <- cbind(
    k * x0,
    (1 - h) * x1,
    -h * x2,
    h * x3
  )
  
  ## D11 = d logit(p11) / d beta
  ## because logit(p11) = logit(p10) + theta2
  D11 <- cbind(
    -k * x0,
    -a.tilde * x1,
    (1 - a.tilde) * x2,
    a.tilde * x3
  )
  
  cell00 <- a0 == 0 & a1 == 0
  cell10 <- a0 == 1 & a1 == 0
  cell01 <- a0 == 0 & a1 == 1
  cell11 <- a0 == 1 & a1 == 1
  
  tmp <- matrix(NA_real_, n, ncol(D00))
  
  tmp[cell00, ] <- D00[cell00, , drop = FALSE]
  tmp[cell10, ] <- D10[cell10, , drop = FALSE]
  tmp[cell01, ] <- D01[cell01, , drop = FALSE]
  tmp[cell11, ] <- D11[cell11, , drop = FALSE]
  
  if (any(!is.finite(pA)) || any(!is.finite(tmp))) {
    stop("Non-finite fitted probabilities or derivatives in Fisher information.",
         call. = FALSE)
  }
  
  info.weight <- pA * (1 - pA)
  
  fisher.info <- t(info.weight * weight * tmp) %*% tmp
  
  var.mat <- if (requireNamespace("MASS", quietly = TRUE)) {
    MASS::ginv(fisher.info)
  } else {
    solve(fisher.info)
  }
  
  par.names <- coef.names(x0, x1, x2, x3)
  
  if (length(par.names) == ncol(var.mat) && all(!is.na(par.names))) {
    rownames(var.mat) <- par.names
    colnames(var.mat) <- par.names
  }
  
  var.mat
}

var.mle.TwoORSRSR <- function(a0, a1,
                              b0.ml, b1.ml, b2.ml, b3.ml,
                              x0, x1, x2, x3,
                              weight = rep(1, length(a0)),
                              eps = 1e-10) {
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  x3 <- as.matrix(x3)
  
  a0 <- as.vector(a0)
  a1 <- as.vector(a1)
  weight <- as.vector(weight)
  
  n <- length(a0)
  
  if (length(a1) != n || length(weight) != n) {
    stop("a0, a1, and weight must have the same length.", call. = FALSE)
  }
  
  if (nrow(x0) != n || nrow(x1) != n || nrow(x2) != n || nrow(x3) != n) {
    stop("x0, x1, x2, and x3 must have the same number of rows as length(a0).",
         call. = FALSE)
  }
  
  theta0 <- as.vector(x0 %*% b0.ml)
  theta1 <- as.vector(x1 %*% b1.ml)
  theta2 <- as.vector(x2 %*% b2.ml)
  phi    <- as.vector(x3 %*% b3.ml)
  
  pa0a1 <- getTwoORSRSR(theta0, theta1, theta2, phi, eps = eps)
  
  if (is.null(colnames(pa0a1))) {
    colnames(pa0a1) <- c("p00", "p10", "p01", "p11")
  }
  
  p00 <- pa0a1[, "p00"]
  p10 <- pa0a1[, "p10"]
  p01 <- pa0a1[, "p01"]
  p11 <- pa0a1[, "p11"]
  
  pA <- observed.prob(pa0a1, a0, a1)
  pA <- clip_prob(pA, eps = eps)
  
  if (any(!is.finite(pa0a1)) || any(!is.finite(pA))) {
    stop("Non-finite fitted probabilities in Fisher information.",
         call. = FALSE)
  }
  
  if (any(pa0a1 <= eps | pa0a1 >= 1 - eps, na.rm = TRUE)) {
    warning("Some fitted cell probabilities are close to 0 or 1; variance may be unstable.")
  }
  
  ## For SR derivatives, denominators are p01 and p11,
  ## not q01 and q11.
  if (any(p01 <= eps | p11 <= eps, na.rm = TRUE)) {
    stop("p01 or p11 is too close to 0; derivatives are unstable.",
         call. = FALSE)
  }
  
  A <- p00 / p01
  B <- p10 / p11
  S <- 2 + A + B
  
  if (any(!is.finite(A)) || any(!is.finite(B)) || any(!is.finite(S)) ||
      any(S <= 0)) {
    stop("Non-finite or invalid derivative coefficients.",
         call. = FALSE)
  }
  
  D00 <- cbind(
    -(1 + B) / S * x0,
    1 / (S * p01) * x1,
    1 / (S * p11) * x2,
    1 / S * x3
  )
  
  D10 <- cbind(
    (1 + A) / S * x0,
    1 / (S * p01) * x1,
    1 / (S * p11) * x2,
    1 / S * x3
  )
  
  D01 <- cbind(
    -A * (1 + B) / S * x0,
    -(S - A) / (S * p01) * x1,
    A / (S * p11) * x2,
    A / S * x3
  )
  
  D11 <- cbind(
    B * (1 + A) / S * x0,
    B / (S * p01) * x1,
    -(S - B) / (S * p11) * x2,
    B / S * x3
  )
  
  cell00 <- a0 == 0 & a1 == 0
  cell10 <- a0 == 1 & a1 == 0
  cell01 <- a0 == 0 & a1 == 1
  cell11 <- a0 == 1 & a1 == 1
  
  tmp <- matrix(NA_real_, n, ncol(D00))
  
  tmp[cell00, ] <- D00[cell00, , drop = FALSE]
  tmp[cell10, ] <- D10[cell10, , drop = FALSE]
  tmp[cell01, ] <- D01[cell01, , drop = FALSE]
  tmp[cell11, ] <- D11[cell11, , drop = FALSE]
  
  if (any(!is.finite(tmp))) {
    stop("Non-finite derivatives in Fisher information.",
         call. = FALSE)
  }
  
  info.weight <- pA * (1 - pA)
  
  fisher.info <- t(info.weight * weight * tmp) %*% tmp
  
  var.mat <- if (requireNamespace("MASS", quietly = TRUE)) {
    MASS::ginv(fisher.info)
  } else {
    solve(fisher.info)
  }
  
  par.names <- coef.names(x0, x1, x2, x3)
  
  if (length(par.names) == ncol(var.mat) && all(!is.na(par.names))) {
    rownames(var.mat) <- par.names
    colnames(var.mat) <- par.names
  }
  
  var.mat
}

var.mle.by.measures <- function(a0, a1,
                                b0.ml, b1.ml, b2.ml, b3.ml,
                                x0, x1, x2, x3,
                                weight = rep(1, length(a0)),
                                measures = c("RR", "OR", "OR"),
                                ...) {
  fun <- measure.function.Two("var.mle.Two", measures)
  fun(
    a0 = a0, a1 = a1,
    b0.ml = b0.ml, b1.ml = b1.ml,
    b2.ml = b2.ml, b3.ml = b3.ml,
    x0 = x0, x1 = x1, x2 = x2, x3 = x3,
    weight = weight,
    ...
  )
}

var.mle.TwoByMeasures <- var.mle.by.measures
