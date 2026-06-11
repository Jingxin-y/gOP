rm(list = ls())
gc

setwd("D:/codexresults/gOP/gop_two_modules")

files <- c(
  "01_utils.R",
  "02_probabilities.R",
  "03_likelihood.R",
  "04_variance.R",
  "05_estimator.R",
  "simuData.R",
  "simu_evaluation.R"
)
for (f in files) source(f)



##simulate gOP data under RROROR model

dataGOP <- simulate_gop_dgm_getTwo(n = 1000,measures = c("RR", "OR", "OR"),f0 = ~ 1,
                        f1 = ~ 1,f2 = ~ 1,f3 = ~ 1 + z1 + z2 + z3,
                        b0.true = c(0.15),b1.true = c(0.20),
                        b2.true = c(0.25),b3.true = c(-6.00,  0.35, -0.20, 0.15),
                        design = c("factorial"),pi0 = 0.5,pi1 = 0.5,
                        randomized = TRUE,eps = 1e-8)



##estimate by RBC model

d1 <- standardize_sim_data(dataGOP)
fit_rbc <- fit_rbc_method(d1)
predict_rbc <- predict_rbc_method(fit_rbc,d1,measures = c("RR", "OR", "OR"))

##estimate by gOP model
fit_gOP <- fit_gop_method(d1)

## compare two methods

for(r in 1000){
  set.seed(r+2026)
  dataGOP <- simulate_gop_dgm_getTwo(n = 1000,measures = c("RR", "OR", "OR"),f0 = ~ 1,
                                     f1 = ~ 1,f2 = ~ 1,f3 = ~ 1 + z1 + z2 + z3,
                                     b0.true = c(0.15),b1.true = c(0.20),
                                     b2.true = c(0.25),b3.true = c(-6.00,  0.35, -0.20, 0.15),
                                     design = c("factorial"),pi0 = 0.5,pi1 = 0.5,
                                     randomized = TRUE,eps = 1e-8)
  d1 <- standardize_sim_data(dataGOP)
  fit_rbc <- fit_rbc_method(d1)
  fit_gOP <- fit_gop_method(d1)
  
  predict_rbc <- predict_rbc_method(fit_rbc,d1,measures = c("RR", "OR", "OR"))
  
  result_rbc <- evalu_rbc_CM(fit_rbc,predict_rbc)
  
  result_gOP <- evalu_gOP_CM(fit_gOP)
  
  result <- list(rbc = result_rbc,
                 gOP = result_gOP)
}
