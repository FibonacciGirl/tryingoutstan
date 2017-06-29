data {
  int<lower =0> J;
  int<lower =0> Tr;
  real y[J,Tr]; 

}
parameters {
  real mu_mean;
  real <lower=0> mu_sd;
  real mu[J]; 
}

model {
  target +=normal_lpdf(mu_mean|0,100);
  target +=gamma_lpdf(mu_sd|100,100);

  for(j in 1:J){
    target += normal_lpdf(mu[j] | mu_mean,mu_sd);
    target += normal_lpdf(y[j,] | mu[j], 10);

  }
}
