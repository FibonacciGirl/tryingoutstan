data {
  int <lower=0> Tr; // number of presentation
  int <lower=0> J; // number of subjects 
  real rtu[J,Tr]; // response time unpredictable
  real rtp[J,Tr]; // response time predictable
}
parameters {
  real <lower=0,upper=1> theta[J]; //
  real <lower=0,upper = 20> sigma_mean; // standard deviation group mean
  real sigma_offset[J]; // the offset of an individual's mean from the group mean
  real u_intercept_mean; // unpredictable intercept group mean
  real u_intercept_offset[J]; // 
  real p_intercept_mean; //predictable intercept group mean
  real p_intercept_offset[J];
  real p_intercept_mean_2; //predictable intercept group mean
  real p_intercept_offset_2[J];
  real p_slope_mean; //predictable slope group mean
  real p_slope_offset[J];

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
    sigma[j] = sigma_mean + sigma_offset[j];
    u_intercept[j] = u_intercept_mean + u_intercept_offset[j];
    p_intercept[j] = p_intercept_mean + p_intercept_offset[j];
    p_intercept_2[j] = p_intercept_mean_2 + p_intercept_offset_2[j];
    p_slope[j] = p_slope_mean + p_slope_offset[j];
    
    for(t in 1:Tr){
    u_mu[j,t] = u_intercept[j];
    p_mu[j,t] = theta[j]*p_intercept[j] + (1-theta[j])*(p_intercept_2[j] + p_slope[j]*t);
    }
  }
}

model {
  //priors
  target += beta_lpdf(theta[J]|1,1);
  target += gamma_lpdf(sigma_mean|8,1);
  target += normal_lpdf(sigma_offset[J]|0,1);
  target += normal_lpdf(u_intercept_mean|100,10);
  target += normal_lpdf(u_intercept_offset[J]|0,1);
  target += normal_lpdf(p_intercept_mean|40,10);
  target += normal_lpdf(p_intercept_offset[J]|0,1);
  target += normal_lpdf(p_intercept_mean_2|-100,10);
  target += normal_lpdf(p_intercept_offset_2[J]|0,1);
  target += normal_lpdf(p_slope_mean|-3,2);
  target += normal_lpdf(p_slope_offset[J]|0,1);
  
  //likelihood
  target += normal_lpdf(rtu[J,Tr] | u_mu[J,Tr], sigma[J]);
  target += normal_lpdf(rtp[J,Tr] | p_mu[J,Tr], sigma[J]);
}
