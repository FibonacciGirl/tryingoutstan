data {
  int<lower=1> Tr;
  real y[Tr];
}
transformed data {
  real log_unif;
  log_unif = -log(Tr);
}
parameters {
  real mu1;
  real mu2;
}
transformed parameters {
  vector[Tr] lp;
  lp = rep_vector(log_unif, Tr);
  for (s in 1:Tr){
    for (t in 1:Tr){
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? mu1 : mu2 , 10);
    }
  }
}
model {
mu1 ~ normal_lpdf(0,100);
mu2 ~ normal_lpdf(0,100);
target += log_sum_exp(lp);
}
