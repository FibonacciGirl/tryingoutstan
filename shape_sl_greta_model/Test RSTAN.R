install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
require(rstan)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())




fit <- stan(file = 'STAN_model_hier.stan', data = model_data, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T )



print(fit,par = c('p_intercept_mean'))
plot(fit,par=c('p_slope[4]'))
theta[4]



