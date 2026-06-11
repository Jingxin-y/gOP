### Basic maps and numeric utilities

clip_prob <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}

fRR <- function(p, RR) {
  p * RR
}

fRD <- function(p, RD) {
  p + RD
}

fSR <- function(p, SR) {
  1 - SR + SR * p
}

fOR <- function(p, OR) {
  OR * p / (1 - p + OR * p)
}

fCHR <- function(p, CHR) {
  1 - (1 - p)^CHR
}

logit <- function(p) {
  p <- clip_prob(p)
  log(p) - log1p(-p)
}

expit <- function(x) {
  ifelse(x >= 0, 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))
}

fOP <- function(p) {
  sum(logit(p))
}
