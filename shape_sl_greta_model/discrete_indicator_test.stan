data {
  int K;
  int I;
  real y[I];
  vector<lower=0>[K] alpha;
  vector[K] mu0;
  real <lower=0> sigma0;
}

parameters {
  simplex[K] pi;
  vector[K] mu;
}
transformed parameters {
  vector[K] log_q_z[I];
  for(i in 1:I){
    log_q_z[i]= log(pi);
    for(k in 1:K){
      log_q_z[i,k]= log_q_z[i,k] +normal_lpdf(y[i]|mu[k],10);
    }
  }
}
model {
    target += dirichlet_lpdf(pi|alpha);
  for(k in 1:K){
    target += normal_lpdf(mu[k]|mu0[k],sigma0);
  }
  for(i in 1:I){
    target += log_sum_exp(log_q_z[i]);
  }
}
