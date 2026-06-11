## ============================================================
## Simulation data generators for comparing RBC-style and gOP-style DGMs
## ============================================================

## ---------- basic utilities ----------

clip_prob <- function(p, eps = 1e-8) {
  p <- as.numeric(p)
  pmin(pmax(p, eps), 1 - eps)
}

safe_logit <- function(p, eps = 1e-8) {
  qlogis(clip_prob(p, eps))
}

## OR flow: logit(p_new) = logit(p_old) + eta
flow_or <- function(p, eta, eps = 1e-8) {
  clip_prob(plogis(safe_logit(p, eps) + eta), eps)
}

## RR flow: p_new = exp(eta) * p_old
flow_rr <- function(p, eta, eps = 1e-8) {
  clip_prob(exp(eta) * p, eps)
}

## SR flow: (1 - p_new) / (1 - p_old) = exp(eta)
## Positive eta means survival is multiplied up, so event risk goes down.
flow_sr <- function(p, eta, eps = 1e-8) {
  p <- clip_prob(p, eps)
  sr <- exp(eta)
  
  ## keep within the probability space
  sr_max <- (1 - eps) / (1 - p)
  sr <- pmin(sr, sr_max)
  
  clip_prob(1 - sr * (1 - p), eps)
}

## RD flow: p_new = p_old + eta
flow_rd <- function(p, eta, eps = 1e-8) {
  clip_prob(p + eta, eps)
}

## Power odds flow: logit(p_new) = gamma * logit(p_old), gamma = exp(eta)
flow_pow_odds <- function(p, eta, eps = 1e-8) {
  gamma <- exp(eta)
  clip_prob(plogis(gamma * safe_logit(p, eps)), eps)
}

## apply a selected effect measure to a reference risk
apply_measure <- function(p_ref, eta, measure, eps = 1e-8) {
  measure <- toupper(measure)
  
  if (measure == "OR") {
    return(flow_or(p_ref, eta, eps))
  }
  
  if (measure == "RR") {
    return(flow_rr(p_ref, eta, eps))
  }
  
  if (measure == "SR") {
    return(flow_sr(p_ref, eta, eps))
  }
  
  if (measure == "RD") {
    return(flow_rd(p_ref, eta, eps))
  }
  
  stop("Unknown measure: ", measure, call. = FALSE)
}

## recover eta from two risks under a specified measure
measure_eta <- function(p_num, p_den, measure, eps = 1e-8) {
  measure <- toupper(measure)
  p_num <- clip_prob(p_num, eps)
  p_den <- clip_prob(p_den, eps)
  
  if (measure == "OR") {
    return(safe_logit(p_num, eps) - safe_logit(p_den, eps))
  }
  
  if (measure == "RR") {
    return(log(p_num / p_den))
  }
  
  if (measure == "SR") {
    return(log((1 - p_num) / (1 - p_den)))
  }
  
  if (measure == "RD") {
    return(p_num - p_den)
  }
  
  stop("Unknown measure: ", measure, call. = FALSE)
}

## treatment assignment
assign_AB <- function(n, Xc1, Xb,
                      design = c("factorial", "A1_only"),
                      randomized = TRUE,
                      piA = 0.5,
                      piB = 0.5) {
  design <- match.arg(design)
  
  if (design == "A1_only") {
    A <- rep(1L, n)
    
    if (randomized) {
      piB_i <- rep(piB, n)
    } else {
      piB_i <- plogis(qlogis(piB) + 0.30 * Xc1 - 0.35 * Xb)
    }
    
    B <- rbinom(n, 1, piB_i)
    
    return(list(
      A = A,
      B = B,
      piA = rep(1, n),
      piB = piB_i
    ))
  }
  
  if (randomized) {
    piA_i <- rep(piA, n)
    piB_i <- rep(piB, n)
  } else {
    piA_i <- plogis(qlogis(piA) + 0.25 * Xc1 - 0.20 * Xb)
    piB_i <- plogis(qlogis(piB) - 0.15 * Xc1 + 0.30 * Xb)
  }
  
  A <- rbinom(n, 1, piA_i)
  B <- rbinom(n, 1, piB_i)
  
  list(
    A = A,
    B = B,
    piA = piA_i,
    piB = piB_i
  )
}

## oracle summary for true risks and true treatment effects
oracle_effect_summary <- function(dat) {
  risk00 <- mean(dat$p00)
  risk10 <- mean(dat$p10)
  risk01 <- mean(dat$p01)
  risk11 <- mean(dat$p11)
  
  data.frame(
    risk00 = risk00,
    risk10 = risk10,
    risk01 = risk01,
    risk11 = risk11,
    
    RD_A_without_B = mean(dat$p10 - dat$p00),
    RD_B_without_A = mean(dat$p01 - dat$p00),
    RD_B_with_A    = mean(dat$p11 - dat$p10),
    
    RR_A_without_B = risk10 / risk00,
    RR_B_without_A = risk01 / risk00,
    RR_B_with_A    = risk11 / risk10,
    
    OR_B_with_A = (risk11 / (1 - risk11)) / (risk10 / (1 - risk10)),
    
    RD_interaction = mean((dat$p11 - dat$p10) - (dat$p01 - dat$p00)),
    RR_interaction = (risk11 / risk10) / (risk01 / risk00)
  )
}

## ---------- RBC-style composition ----------

rbc_compose <- function(p,
                        eta_or  = 0,
                        eta_rr  = 0,
                        eta_sr  = 0,
                        eta_pow = 0,
                        order = c("OR", "RR", "SR", "PO"),
                        eps = 1e-8) {
  p_now <- clip_prob(p, eps)
  order <- toupper(order)
  
  for (op in order) {
    if (op == "OR") {
      p_now <- flow_or(p_now, eta_or, eps)
    } else if (op == "RR") {
      p_now <- flow_rr(p_now, eta_rr, eps)
    } else if (op == "SR") {
      p_now <- flow_sr(p_now, eta_sr, eps)
    } else if (op %in% c("PO", "POWODDS", "POWER_ODDS")) {
      p_now <- flow_pow_odds(p_now, eta_pow, eps)
    } else {
      stop("Unknown RBC flow: ", op, call. = FALSE)
    }
  }
  
  clip_prob(p_now, eps)
}

simulate_rbc_dgm <- function(n,
                             seed = NULL,
                             design = c("factorial", "A1_only"),
                             randomized = TRUE,
                             piA = 0.5,
                             piB = 0.5) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  design <- match.arg(design)
  
  ## Covariates: two continuous and one binary
  Xc1 <- rnorm(n)
  Xc2 <- 0.5 * Xc1 + rnorm(n, sd = 0.8)
  Xb  <- rbinom(n, 1, 0.45)
  
  ## Baseline risk: p00 = risk without AZT and without ddI
  p00 <- plogis(-1.25 + 0.70 * Xc1 - 0.45 * Xb + 0.25 * Xc2)
  
  ## AZT effect without ddI: p10 vs p00
  p10 <- rbc_compose(
    p00,
    eta_or  = -0.18 + 0.10 * Xc1,
    eta_rr  = -0.08 - 0.04 * Xb,
    eta_sr  =  0.04 + 0.02 * Xb,
    eta_pow =  0.04 + 0.02 * Xc2,
    order   = c("OR", "RR", "SR", "PO")
  )
  
  ## ddI effect without AZT: p01 vs p00
  p01 <- rbc_compose(
    p00,
    eta_or  = -0.10 - 0.08 * Xb,
    eta_rr  = -0.04 + 0.03 * Xc1,
    eta_sr  =  0.03,
    eta_pow =  0.02 + 0.02 * Xc2,
    order   = c("OR", "RR", "SR", "PO")
  )
  
  ## ddI effect with AZT: p11 vs p10
  ## This makes ddI's incremental effect different when AZT is present.
  p11 <- rbc_compose(
    p10,
    eta_or  = -0.22 + 0.12 * Xc1 - 0.10 * Xb,
    eta_rr  = -0.06,
    eta_sr  =  0.05 + 0.02 * Xb,
    eta_pow =  0.04 + 0.02 * Xc2,
    order   = c("OR", "RR", "SR", "PO")
  )
  
  trt_assign <- assign_AB(
    n = n,
    Xc1 = Xc1,
    Xb = Xb,
    design = design,
    randomized = randomized,
    piA = piA,
    piB = piB
  )
  
  A <- trt_assign$A
  B <- trt_assign$B
  
  p_obs <- observed.prob(
    cbind(p00 = p00, p10 = p10, p01 = p01, p11 = p11),
    A,
    B
  )
  Y <- rbinom(n, 1, p_obs)
  
  dat <- data.frame(
    id = seq_len(n),
    Y = Y,
    A = A,       # AZT
    B = B,       # ddI
    arm = paste0("A", A, "B", B),
    
    ## for two-arm AZT-only comparison:
    ## trt = 0 means AZT only, trt = 1 means AZT + ddI
    trt = ifelse(A == 1, B, NA_integer_),
    
    Xc1 = Xc1,
    Xc2 = Xc2,
    Xb = Xb,
    
    piA = trt_assign$piA,
    piB = trt_assign$piB,
    
    p00 = p00,
    p10 = p10,
    p01 = p01,
    p11 = p11,
    p_obs = p_obs,
    
    true_RD_B_with_A = p11 - p10,
    true_RR_B_with_A = p11 / p10,
    true_RD_B_without_A = p01 - p00,
    true_RR_B_without_A = p01 / p00
  )
  
  attr(dat, "dgm") <- "RBC_style"
  attr(dat, "oracle") <- oracle_effect_summary(dat)
  
  dat
}

## ---------- gOP-style / measure-composition DGM ----------

make_gop_covariates <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  z1 <- rnorm(n)
  z2 <- rbinom(n, 1, 0.45)
  z3 <- 0.5 * z1 + rnorm(n, sd = 0.8)
  
  data.frame(z1 = z1, z2 = z2, z3 = z3)
}

getTwo_from_measures <- function(theta0, theta1, theta2, phi, measures) {
  if (exists("normalize.measures.Two", mode = "function")) {
    measures <- normalize.measures.Two(measures)
  } else {
    measures <- toupper(unname(measures))
  }
  
  if (exists("getTwo.by.measures", mode = "function")) {
    return(getTwo.by.measures(theta0, theta1, theta2, phi, measures = measures))
  }
  
  key <- paste0(toupper(unname(measures)), collapse = "")
  fname <- paste0("getTwo", key)
  
  if (!exists(fname, mode = "function")) {
    stop("Cannot find ", fname, "(). Please source your getTwo functions first.",
         call. = FALSE)
  }
  
  get(fname, mode = "function")(theta0, theta1, theta2, phi)
}


simulate_gop_dgm_getTwo <- function(n = 1000,
                                    seed = NULL,
                                    covariates = NULL,
                                    measures = c("RR", "OR", "OR"),
                                    
                                    f0 = ~ 1 + z1 + z2,
                                    f1 = ~ 1 + z1 + z2,
                                    f2 = ~ 1 + z1 + z2,
                                    f3 = ~ 1 + z1 + z2 + z3,
                                    
                                    b0.true = c(-0.15,  0.08, -0.05),
                                    b1.true = c(-0.20, -0.05,  0.08),
                                    b2.true = c(-0.25,  0.10, -0.08),
                                    b3.true = c(-6.00,  0.35, -0.20, 0.15),
                                    
                                    design = c("factorial", "A1_only"),
                                    pi0 = 0.5,
                                    pi1 = 0.5,
                                    randomized = TRUE,
                                    eps = 1e-8) {
  if (!is.null(seed)) set.seed(seed)
  
  design <- match.arg(design)
  
  if (is.null(covariates)) {
    dat_x <- make_gop_covariates(n)
  } else {
    dat_x <- as.data.frame(covariates)
    n <- nrow(dat_x)
  }
  
  X0 <- model.matrix(f0, dat_x)
  X1 <- model.matrix(f1, dat_x)
  X2 <- model.matrix(f2, dat_x)
  X3 <- model.matrix(f3, dat_x)
  
  if (length(b0.true) != ncol(X0)) stop("b0.true length does not match X0.")
  if (length(b1.true) != ncol(X1)) stop("b1.true length does not match X1.")
  if (length(b2.true) != ncol(X2)) stop("b2.true length does not match X2.")
  if (length(b3.true) != ncol(X3)) stop("b3.true length does not match X3.")
  
  theta0 <- as.vector(X0 %*% b0.true)
  theta1 <- as.vector(X1 %*% b1.true)
  theta2 <- as.vector(X2 %*% b2.true)
  phi    <- as.vector(X3 %*% b3.true)
  
  pall <- as.matrix(getTwo_from_measures(theta0, theta1, theta2, phi, measures))
  colnames(pall) <- c("p00", "p10", "p01", "p11")
  
  if (any(!is.finite(pall)) || any(pall <= eps) || any(pall >= 1 - eps)) {
    stop("Generated probabilities are invalid. Try smaller coefficients or adjust b3.true.",
         call. = FALSE)
  }
  
  if (design == "factorial") {
    if (randomized) {
      prob0 <- rep(pi0, n)
      prob1 <- rep(pi1, n)
    } else {
      prob0 <- plogis(qlogis(pi0) + 0.25 * dat_x$z1 - 0.20 * dat_x$z2)
      prob1 <- plogis(qlogis(pi1) - 0.15 * dat_x$z1 + 0.25 * dat_x$z2)
    }
    
    a0 <- rbinom(n, 1, prob0)
    a1 <- rbinom(n, 1, prob1)
  } else {
    a0 <- rep(1L, n)
    prob0 <- rep(1, n)
    
    if (randomized) {
      prob1 <- rep(pi1, n)
    } else {
      prob1 <- plogis(qlogis(pi1) + 0.25 * dat_x$z1 - 0.20 * dat_x$z2)
    }
    
    a1 <- rbinom(n, 1, prob1)
  }
  
  p_obs <- observed.prob(pall, a0, a1)
  
  y <- rbinom(n, 1, p_obs)
  
  out <- data.frame(
    id = seq_len(n),
    y = y,
    a0 = a0,
    a1 = a1,
    trt = ifelse(a0 == 1, a1, NA_integer_),
    dat_x,
    theta0 = theta0,
    theta1 = theta1,
    theta2 = theta2,
    phi = phi,
    p00 = pall[, "p00"],
    p10 = pall[, "p10"],
    p01 = pall[, "p01"],
    p11 = pall[, "p11"],
    p_obs = p_obs,
    RD_11_10 = pall[, "p11"] - pall[, "p10"],
    RR_11_10 = pall[, "p11"] / pall[, "p10"]
  )
  
  attr(out, "X0") <- X0
  attr(out, "X1") <- X1
  attr(out, "X2") <- X2
  attr(out, "X3") <- X3
  attr(out, "measures") <- measures
  attr(out, "true_beta") <- list(b0 = b0.true, b1 = b1.true, b2 = b2.true, b3 = b3.true)
  
  out
}
