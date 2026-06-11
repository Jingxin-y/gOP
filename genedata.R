## generate RBC data

dat_rbc_sim <- simulate_rbc_dgm(
  n = 1000,
  seed = 1,
  design = "factorial",
  randomized = TRUE
)

head(dat_rbc_sim)

attr(dat_rbc_sim, "oracle")
table(dat_rbc_sim$arm)

## generete gOP data

dat_gop <- simulate_gop_dgm_getTwo(
  n = 1000,
  seed = 2026,
  measures = c("RR", "OR", "OR"),
  design = "factorial"
)

summary(dat_gop[, c("p00", "p10", "p01", "p11")])
mean(dat_gop$p11 - dat_gop$p10)
mean(dat_gop$p11) / mean(dat_gop$p10)