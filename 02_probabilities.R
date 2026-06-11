### Cell-probability construction for supported two-time-point gOP models

normalize.measures.Two <- function(measures = c("RR", "OR", "OR")) {
  target.names <- c("p10_p00", "p01_p00", "p11_p10")
  if (is.null(measures)) {
    measures <- c("RR", "OR", "OR")
  }
  measures <- as.character(measures)

  if (length(measures) == 1) {
    compact <- toupper(gsub("[^A-Za-z]", "", measures))
    parsed <- regmatches(compact, gregexpr("OR|RR|SR", compact))[[1]]
    if (length(parsed) == 3 && paste0(parsed, collapse = "") == compact) {
      measures <- parsed
    }
  }

  if (length(measures) != 3) {
    stop("measures must have length 3, ordered as c('p10_p00', 'p01_p00', 'p11_p10').",
         call. = FALSE)
  }

  if (is.null(names(measures)) || any(names(measures) == "")) {
    names(measures) <- target.names
  } else {
    missing.names <- setdiff(target.names, names(measures))
    if (length(missing.names) > 0) {
      stop("Named measures must include: p10_p00, p01_p00, p11_p10.",
           call. = FALSE)
    }
    measures <- measures[target.names]
  }

  nm <- names(measures)
  measures <- toupper(unname(measures))
  if (!all(measures %in% c("OR", "RR", "SR"))) {
    stop("measures must use OR, RR, or SR for these closed-form gOP functions.",
         call. = FALSE)
  }
  names(measures) <- nm
  measures
}

measure.key.Two <- function(measures = c("RR", "OR", "OR")) {
  paste0(unname(normalize.measures.Two(measures)), collapse = "")
}

supported.measure.keys.Two <- function() {
  c("RROROR", "ORRRRR", "SROROR", "ORSRSR")
}

supported.measures.Two <- function() {
  list(
    RROROR = c(p10_p00 = "RR", p01_p00 = "OR", p11_p10 = "OR"),
    ORRRRR = c(p10_p00 = "OR", p01_p00 = "RR", p11_p10 = "RR"),
    SROROR = c(p10_p00 = "SR", p01_p00 = "OR", p11_p10 = "OR"),
    ORSRSR = c(p10_p00 = "OR", p01_p00 = "SR", p11_p10 = "SR")
  )
}

measure.suffix.Two <- function(measures = c("RR", "OR", "OR")) {
  key <- measure.key.Two(measures)
  if (!key %in% supported.measure.keys.Two()) {
    supported <- vapply(
      supported.measures.Two(),
      function(x) paste0("c('", paste(unname(x), collapse = "', '"), "')"),
      character(1)
    )
    stop(
      "Unsupported measure combination for the closed-form gOP functions: ",
      paste(unname(normalize.measures.Two(measures)), collapse = "/"),
      ". Supported combinations are: ",
      paste(supported, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  key
}

measure.function.Two <- function(prefix, measures = c("RR", "OR", "OR")) {
  fname <- paste0(prefix, measure.suffix.Two(measures))
  if (!exists(fname, mode = "function")) {
    stop("Function not found: ", fname, call. = FALSE)
  }
  get(fname, mode = "function")
}

getcellRROROR <- function(p00, theta0, theta1, theta2) {
  p10 <- fRR(p00, exp(theta0))
  p01 <- fOR(p00, exp(theta1))
  p11 <- fOR(p10, exp(theta2))
  
  cbind(p00, p10, p01, p11)
}

getcellORRRRR <- function(p00, theta0, theta1, theta2) {
  p10 <- fOR(p00, exp(theta0))  # OR(p10, p00)
  p01 <- fRR(p00, exp(theta1))  # RR(p01, p00)
  p11 <- fRR(p10, exp(theta2))  # RR(p11, p10)
  
  cbind(p00 = p00, p10 = p10, p01 = p01, p11 = p11)
}


getcellSROROR <- function(p00, theta0, theta1, theta2) {
  p10 <- fSR(p00, exp(theta0))
  p01 <- fOR(p00, exp(theta1))
  p11 <- fOR(p10, exp(theta2))
  
  cbind(p00, p10, p01, p11)
}

getcellORSRSR <- function(p00, theta0, theta1, theta2) {
  p00 <- as.vector(p00)
  theta0 <- as.vector(theta0)
  theta1 <- as.vector(theta1)
  theta2 <- as.vector(theta2)
  
  a <- exp(theta0)  # OR(p10, p00)
  b <- exp(theta1)  # SR(p01, p00)
  c <- exp(theta2)  # SR(p11, p10)
  
  p10 <- a * p00 / (1 - p00 + a * p00)
  p01 <- 1 - b * (1 - p00)
  p11 <- 1 - c * (1 - p10)
  
  cbind(
    p00 = p00,
    p10 = p10,
    p01 = p01,
    p11 = p11
  )
}

getTwoRROROR <- function(theta0, theta1, theta2, phi, eps = 1e-10) {
  theta0 <- as.vector(theta0)
  theta1 <- as.vector(theta1)
  theta2 <- as.vector(theta2)
  phi <- as.vector(phi)
  n <- length(theta0)
  if (length(theta1) != n || length(theta2) != n || length(phi) != n) {
    stop("theta0, theta1, theta2, and phi must have the same length.", call. = FALSE)
  }
  
  r <- exp(theta0)
  m <- exp((phi - theta1 - theta2) / 2)
  p00 <- rep(NA_real_, n)
  
  use.limit <- abs(m - 1) < 1e-8
  p00[use.limit] <- 1 / (1 + r[use.limit])
  
  use.quad <- !use.limit
  if (any(use.quad)) {
    disc <- m[use.quad]^2 * (1 + r[use.quad])^2 +
      4 * r[use.quad] * m[use.quad] * (1 - m[use.quad])
    p00[use.quad] <- (-m[use.quad] * (1 + r[use.quad]) +
      sqrt(pmax(disc, 0))) / (2 * r[use.quad] * (1 - m[use.quad]))
  }
  
  bad <- !is.finite(p00) | p00 <= eps | p00 >= 1 - eps |
    r * p00 >= 1 - eps
  if (any(bad)) {
    for (i in which(bad)) {
      lower <- eps
      upper <- 1 - eps
      if (theta0[i] > 0) {
        upper <- min(upper, exp(-theta0[i]) - eps)
      }
      if (!is.finite(upper) || upper <= lower) {
        next
      }
      root_fun <- function(x) {
        sum(logit(getcellRROROR(x, theta0[i], theta1[i], theta2[i]))) - phi[i]
      }
      p00[i] <- tryCatch(
        uniroot(root_fun, interval = c(lower, upper), tol = 1e-10)$root,
        error = function(e) NA_real_
      )
    }
  }
  
  getcellRROROR(p00, theta0, theta1, theta2)
}

getTwoORRRRR <- function(theta0, theta1, theta2, phi,
                         eps = 1e-10, tol = 1e-10) {
  theta0 <- as.vector(theta0)
  theta1 <- as.vector(theta1)
  theta2 <- as.vector(theta2)
  phi    <- as.vector(phi)
  
  n <- length(theta0)
  if (length(theta1) != n || length(theta2) != n || length(phi) != n) {
    stop("theta0, theta1, theta2, and phi must have the same length.",
         call. = FALSE)
  }
  
  p00 <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    a <- exp(theta0[i])  # OR(p10, p00)
    b <- exp(theta1[i])  # RR(p01, p00)
    c <- exp(theta2[i])  # RR(p11, p10)
    
    upper.raw <- 1
    
    # Constraint from p01 = b * p00 < 1
    if (b > 1) {
      upper.raw <- min(upper.raw, 1 / b)
    }
    
    # Constraint from p11 = c * p10 < 1
    if (c > 1) {
      upper.raw <- min(upper.raw, 1 / (1 + a * (c - 1)))
    }
    
    if (!is.finite(upper.raw) || upper.raw <= 0) {
      next
    }
    
    # Avoid evaluating exactly on the boundary
    lower <- min(eps, upper.raw * 1e-8)
    upper <- upper.raw - min(eps, upper.raw * 1e-8)
    
    if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
      next
    }
    
    root_fun <- function(x) {
      pp <- getcellORRRRR(x, theta0[i], theta1[i], theta2[i])
      
      # Use raw logit calculation here; x is chosen inside the valid interval.
      sum(logit(pp)) - phi[i]
    }
    
    p00[i] <- tryCatch(
      uniroot(root_fun, interval = c(lower, upper), tol = tol)$root,
      error = function(e) NA_real_
    )
  }
  
  getcellORRRRR(p00, theta0, theta1, theta2)
}

## 如果你前面已经定义过 fOR / fSR，可以省略这两个
fOR <- function(p0, OR) {
  OR * p0 / (1 - p0 + OR * p0)
}

fSR <- function(p0, SR) {
  1 - SR * (1 - p0)
}


getTwoSROROR <- function(theta0, theta1, theta2, phi, eps = 1e-10) {
  theta0 <- as.vector(theta0)
  theta1 <- as.vector(theta1)
  theta2 <- as.vector(theta2)
  phi    <- as.vector(phi)
  
  n <- length(theta0)
  
  if (length(theta1) != n || length(theta2) != n || length(phi) != n) {
    stop("theta0, theta1, theta2, and phi must have the same length.",
         call. = FALSE)
  }
  
  s <- exp(theta0)
  
  ## m = odds(p00) * odds(p10)
  m <- exp((phi - theta1 - theta2) / 2)
  
  p00 <- rep(NA_real_, n)
  
  ## closed-form stable solution:
  ## p00 = 2*m*s / {1 - s + 2*m*s + sqrt((1-s)^2 + 4*m*s)}
  ok <- is.finite(s) & is.finite(m) & s > 0 & m > 0
  
  if (any(ok)) {
    disc <- (1 - s[ok])^2 + 4 * m[ok] * s[ok]
    denom <- 1 - s[ok] + 2 * m[ok] * s[ok] + sqrt(pmax(disc, 0))
    
    p00[ok] <- 2 * m[ok] * s[ok] / denom
  }
  
  ## 检查 p00 和 p10 是否在合法区间内
  p10 <- fSR(p00, s)
  
  bad <- !is.finite(p00) | p00 <= eps | p00 >= 1 - eps |
    !is.finite(p10) | p10 <= eps | p10 >= 1 - eps
  
  ## 对数值极端的情况，用 uniroot 兜底
  if (any(bad)) {
    for (i in which(bad)) {
      si <- exp(theta0[i])
      
      if (!is.finite(si) || si <= 0) {
        next
      }
      
      ## p00 本身在 (0,1)
      lower <- eps
      upper <- 1 - eps
      
      ## 还要保证 p10 = 1 - s + s*p00 在 (eps, 1-eps)
      lower <- max(lower, (eps + si - 1) / si)
      upper <- min(upper, 1 - eps / si)
      
      if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
        next
      }
      
      root_fun <- function(x) {
        cell <- getcellSROROR(x, theta0[i], theta1[i], theta2[i])
        sum(logit(cell)) - phi[i]
      }
      
      p00[i] <- tryCatch(
        uniroot(root_fun, interval = c(lower, upper), tol = 1e-10)$root,
        error = function(e) NA_real_
      )
    }
  }
  
  getcellSROROR(p00, theta0, theta1, theta2)
}


getTwoORSRSR <- function(theta0, theta1, theta2, phi,
                         eps = 1e-10, tol = 1e-10) {
  theta0 <- as.vector(theta0)
  theta1 <- as.vector(theta1)
  theta2 <- as.vector(theta2)
  phi    <- as.vector(phi)
  
  n <- length(theta0)
  if (length(theta1) != n || length(theta2) != n || length(phi) != n) {
    stop("theta0, theta1, theta2, and phi must have the same length.",
         call. = FALSE)
  }
  
  p00 <- rep(NA_real_, n)
  
  for (i in seq_len(n)) {
    a <- exp(theta0[i])  # OR(p10, p00)
    b <- exp(theta1[i])  # SR(p01, p00)
    c <- exp(theta2[i])  # SR(p11, p10)
    
    if (!is.finite(a) || !is.finite(b) || !is.finite(c) || !is.finite(phi[i])) {
      next
    }
    
    lower.raw <- 0
    upper.raw <- 1
    
    # Constraint from p01 = 1 - b * (1 - p00) > 0
    # Only active when b > 1
    if (b > 1) {
      lower.raw <- max(lower.raw, 1 - 1 / b)
    }
    
    # Constraint from p11 = 1 - c * (1 - p10) > 0
    # Only active when c > 1
    # p10 = a * p00 / (1 - p00 + a * p00)
    # p10 > 1 - 1/c  <=>  p00 > (c - 1) / (a + c - 1)
    if (c > 1) {
      lower.raw <- max(lower.raw, (c - 1) / (a + c - 1))
    }
    
    if (!is.finite(lower.raw) || lower.raw < 0 || lower.raw >= upper.raw) {
      next
    }
    
    # Avoid evaluating exactly on the boundary
    width <- upper.raw - lower.raw
    delta <- min(eps, width * 1e-8)
    
    lower <- lower.raw + delta
    upper <- upper.raw - delta
    
    if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
      next
    }
    
    root_fun <- function(x) {
      pp <- getcellORSRSR(x, theta0[i], theta1[i], theta2[i])
      
      if (any(!is.finite(pp)) || any(pp <= 0 | pp >= 1)) {
        return(NA_real_)
      }
      
      sum(qlogis(pp)) - phi[i]
    }
    
    p00[i] <- tryCatch({
      f.lower <- root_fun(lower)
      f.upper <- root_fun(upper)
      
      if (!is.finite(f.lower) || !is.finite(f.upper)) {
        NA_real_
      } else if (f.lower == 0) {
        lower
      } else if (f.upper == 0) {
        upper
      } else if (f.lower * f.upper > 0) {
        NA_real_
      } else {
        uniroot(root_fun, interval = c(lower, upper), tol = tol)$root
      }
    }, error = function(e) NA_real_)
  }
  
  getcellORSRSR(p00, theta0, theta1, theta2)
}

getcell.by.measures <- function(p00, theta0, theta1, theta2,
                                measures = c("RR", "OR", "OR")) {
  fun <- measure.function.Two("getcell", measures)
  fun(p00, theta0, theta1, theta2)
}

getTwo.by.measures <- function(theta0, theta1, theta2, phi,
                               measures = c("RR", "OR", "OR"),
                               ...) {
  fun <- measure.function.Two("getTwo", measures)
  fun(theta0, theta1, theta2, phi, ...)
}

observed.prob <- function(pa0a1, a0, a1) {
  a0 <- as.vector(a0)
  a1 <- as.vector(a1)
  pA <- rep(NA_real_, length(a0))
  cell00 <- a0 == 0 & a1 == 0
  cell10 <- a0 == 1 & a1 == 0
  cell01 <- a0 == 0 & a1 == 1
  cell11 <- a0 == 1 & a1 == 1
  pA[cell00] <- pa0a1[cell00, 1]
  pA[cell10] <- pa0a1[cell10, 2]
  pA[cell01] <- pa0a1[cell01, 3]
  pA[cell11] <- pa0a1[cell11, 4]
  clip_prob(pA)
}

### Backward-compatible names used by older scripts/check files.
observed.prob.TwoRROROR <- observed.prob


