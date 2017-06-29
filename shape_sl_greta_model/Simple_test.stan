data {

  real y[100];
}
parameters {
  real mu; 
}
model {
  target += normal_lpdf(mu | 0,100);
  target += normal_lpdf(y | mu, 10);
}
