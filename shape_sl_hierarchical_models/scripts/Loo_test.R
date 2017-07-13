library(loo)
library(rstan)
install.packages('rstanarm')
library(rstanarm)

fit<-stan_summry(fit.vec.9)
ll<-log_lik(fit)

source('fit_models.R')

log_lik_1<-extract_log_lik(fit.vec.4)
loo_1<- loo(log_lik_1)
print(loo_1)
