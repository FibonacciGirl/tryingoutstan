data {
  int Tr;
  int J;
  real y[J,Tr];  

}
parameters {
  real mu_mean;
  real mu_sd;
  real mu_offset[J];
}
transformed parameters {
  real mu[J];
  for(j in 1:J){
  mu[j] = mu_mean + mu_offset[j]*mu_sd;
  }
}
model {
  target +=normal_lpdf(mu_mean|0,100);
  target +=gamma_lpdf(mu_sd|100,100);

  for(j in 1:J){
    target += normal_lpdf(mu[j] | mu_mean,mu_sd);
    target += normal_lpdf(y[j,] | mu[j], 10);
 
  }
}
