data {
  int <lower=0> Tr; // number of presentation
  int <lower=0> J; // number of subjects 
  real rtu[J,Tr]; // response time unpredictable
  real rtp[J,Tr]; // response time predictable
}
parameters {
  real <lower=0,upper=1> theta[J]; //
  real <lower=0> sigma_mean; // standard deviation group mean
  real sigma_offset[J]; // the offset of an individual's mean from the group mean
  real <lower=0>sigma_variance;
  real u_intercept_mean; // unpredictable intercept group mean
  real u_intercept_offset[J]; // 
  real <lower=0>u_intercept_variance;
  real p_intercept_mean; //predictable intercept group mean
  real p_intercept_offset[J];
  real <lower=0> p_intercept_variance;
  real p_intercept_mean_2; //predictable intercept group mean
  real p_intercept_offset_2[J];
  real <lower=0> p_intercept_variance_2;
  real p_slope_mean; //predictable slope group mean
  real p_slope_offset[J];
  real <lower=0> p_slope_variance;


 
}
transformed parameters {
  real u_mu[J,Tr];
  real p_mu[J,Tr];
  real <lower = 0> sigma[J];
  real u_intercept[J];
  real p_intercept[J];
  real p_intercept_2[J];
  real p_slope[J];
  
  for(j in 1:J){
    sigma[j] = sigma_mean+sigma_offset[j]*sigma_variance;
    u_intercept[j] = u_intercept_mean+u_intercept_offset[j]*u_intercept_variance;
    p_intercept[j] = p_intercept_mean+p_intercept_offset[j]*p_intercept_variance;
    p_intercept_2[j] = p_intercept_mean_2+p_intercept_offset_2[j]*p_intercept_variance_2;
    p_slope[j] = p_slope_mean+p_slope_offset[j]*p_slope_variance;
    
    for(t in 1:Tr){
    u_mu[j,t] = u_intercept[j];
    p_mu[j,t] = theta[j]*p_intercept[j] + (1-theta[j])*(p_intercept_2[j] + p_slope[j]*t);
    }
  }
}

model {
  //priors
  target += beta_lpdf(theta|1,1);
  target += gamma_lpdf(sigma_mean|9,.5);
  target += normal_lpdf(sigma_offset|0,1);
  target += normal_lpdf(u_intercept_mean|400,10);
  target += normal_lpdf(u_intercept_offset|0,1);
  target += normal_lpdf(p_intercept_mean|-400,10);
  target += normal_lpdf(p_intercept_offset|0,1);
  target += normal_lpdf(p_intercept_mean_2|800,10);
  target += normal_lpdf(p_intercept_offset_2|0,1);
  target += normal_lpdf(p_slope_mean|-15,2.5);
  target += normal_lpdf(p_slope_offset|0,1);
  target += gamma_lpdf(sigma_variance|16,2);
  target += gamma_lpdf(u_intercept_variance|26,.5);
  target += gamma_lpdf(p_intercept_variance|26,.5);
  target += gamma_lpdf(p_intercept_variance_2|26,.5);
  target += gamma_lpdf(p_slope_variance|8,3);
  
  //likelihood
  for(j in 1:J){
    for(t in 1:Tr){
      target += normal_lpdf(rtu[j,t] | u_mu[j,t], sigma[j]);
      target += normal_lpdf(rtp[j,t] | p_mu[j,t], sigma[j]);
    }
  }
}
