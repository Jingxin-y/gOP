### Maximum likelihood optimization

maxlike.TwoByMeasures <- function(y, a0, a1, x0, x1, x2, x3,
                                  b0.start, b1.start, b2.start, b3.start,
                                  weight = rep(1, length(y)),
                                  measures = c("RR", "OR", "OR"),
                                  max.step = 1000, thres = 1e-8) {
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  x3 <- as.matrix(x3)
  y <- as.vector(y)
  a0 <- as.vector(a0)
  a1 <- as.vector(a1)
  weight <- as.vector(weight)
  measures <- normalize.measures.Two(measures)
  
  q0 <- length(b0.start)
  q1 <- length(b1.start)
  q2 <- length(b2.start)
  q3 <- length(b3.start)
  
  neg.log.likelihood.current <- function(b0, b1, b2, b3) {
    theta0 <- x0 %*% b0
    theta1 <- x1 %*% b1
    theta2 <- x2 %*% b2
    phi <- x3 %*% b3
    
    pa0a1 <- getTwo.by.measures(
      theta0, theta1, theta2, phi,
      measures = measures
    )
    pA <- observed.prob(pa0a1, a0, a1)
    
    if (any(!is.finite(pA))) {
      return(1e100)
    }
    -sum(weight * (y * log(pA) + (1 - y) * log1p(-pA)))
  }
  
  neg.log.likelihood <- function(pars) {
    b0 <- pars[1:q0]
    b1 <- pars[(q0 + 1):(q0 + q1)]
    b2 <- pars[(q0 + q1 + 1):(q0 + q1 + q2)]
    b3 <- pars[(q0 + q1 + q2 + 1):(q0 + q1 + q2 + q3)]
    neg.log.likelihood.current(b0, b1, b2, b3)
  }
  
  neg.log.likelihood.b0 <- function(b0) neg.log.likelihood.current(b0, b1, b2, b3)
  neg.log.likelihood.b1 <- function(b1) neg.log.likelihood.current(b0, b1, b2, b3)
  neg.log.likelihood.b2 <- function(b2) neg.log.likelihood.current(b0, b1, b2, b3)
  neg.log.likelihood.b3 <- function(b3) neg.log.likelihood.current(b0, b1, b2, b3)
  
  Diff <- function(x, y) sum((x - y)^2) / sum(x^2 + thres)
  b0 <- b0.start
  b1 <- b1.start
  b2 <- b2.start
  b3 <- b3.start
  diff <- thres + 1
  step <- 0
  while (diff > thres & step < max.step) {
    step <- step + 1
    opt0 <- stats::optim(b0, neg.log.likelihood.b0, control = list(maxit = max(100, max.step / 10)))
    diff0 <- Diff(opt0$par, b0)
    b0 <- opt0$par
    opt1 <- stats::optim(b1, neg.log.likelihood.b1, control = list(maxit = max(100, max.step / 10)))
    diff1 <- Diff(opt1$par, b1)
    b1 <- opt1$par
    opt2 <- stats::optim(b2, neg.log.likelihood.b2, control = list(maxit = max(100, max.step / 10)))
    diff2 <- Diff(opt2$par, b2)
    b2 <- opt2$par
    opt3 <- stats::optim(b3, neg.log.likelihood.b3, control = list(maxit = max(100, max.step / 10)))
    diff3 <- Diff(opt3$par, b3)
    b3 <- opt3$par
    diff <- max(diff0, diff1, diff2, diff3)
  }
  
  list(
    par = c(b0, b1, b2, b3),
    measures = measures,
    convergence = (step < max.step),
    value = neg.log.likelihood(c(b0, b1, b2, b3)),
    step = step
  )
}

maxlike.TwoRROROR <- function(...) {
  maxlike.TwoByMeasures(..., measures = c("RR", "OR", "OR"))
}

maxlike.TwoORRRRR <- function(...) {
  maxlike.TwoByMeasures(..., measures = c("OR", "RR", "RR"))
}

maxlike.TwoSROROR <- function(...) {
  maxlike.TwoByMeasures(..., measures = c("SR", "OR", "OR"))
}

maxlike.TwoORSRSR <- function(...) {
  maxlike.TwoByMeasures(..., measures = c("OR", "SR", "SR"))
}

maxlike.TwoMeasure <- maxlike.TwoByMeasures


