data{
  //data
  int Tr; //number of presentations per subject
  int J; //number of subjects
  real <lower = 0> rtu[J,Tr]; // unpredictable response time
  real <lower = 0> rtp[J,Tr]; // unpredictable response time
  
  //prior params
  real <lower = 0> alpha;
  real <lower = 0> beta;
  
  real <lower = 0> mu0;
  real <lower = 0> tau0;
  real <lower = 0> mu1;
  real <lower = 0> tau1;
  real <lower = 0> mu2;
  real <lower = 0> tau2;
  real <lower = 0> alpha0;
  real <lower = 0> beta0;
  real <lower = 0> alpha1;
  real <lower = 0> beta1;
  real <lower = 0> alpha2;
  real <lower = 0> beta2;
  real <lower = 0> alpha3;
  real <lower = 0> beta3;
    
  real <lower = 0> a;
  real <lower = 0> b;
  real <lower = 0> a0;
  real <lower = 0> b0;
  real <lower = 0> a1;
  real <lower = 0> b1;
  real <lower = 0> a2;
  real <lower = 0> b2;
  real <lower = 0> a3;
  real <lower = 0> b3;
  real <lower = 0> a5;
  real <lower = 0> b5;
  real <lower = 0> a6;
  real <lower = 0> b6;
}
parameters{

  vector <lower = 0>[J] sigma;
  vector <lower = 0>[J] rtu_intercept;
  vector <lower = 0, upper = 1>[J] rtu_base;
  vector <lower = 0>[J] rtu_rate;
  vector <lower = 0>[J] rlr_intercept;
  vector <lower = 0>[J] rlr_rate;
  vector <lower = 0, upper = 1>[J] rlr_base;

  real <lower = 0> sigma_mean;
  real <lower = 0> sigma_var;
  
  real <lower = 0> rtu_intercept_mean;
  real <lower = 0> rtu_intercept_var;
  
  real <lower = 0, upper = 1> rtu_base_mean;
  real <lower = 0> rtu_base_var;
  
  real <lower = 0> rtu_rate_mean;
  real <lower = 0> rtu_rate_var;
  
  real <lower = 0>  rlr_intercept_mean;
  real <lower = 0>  rlr_intercept_var;
  real <lower = 0>  rlr_rate_mean;
  real <lower = 0>  rlr_rate_var;
  real <lower = 0>  rlr_base_mean;
  real <lower = 0>  rlr_base_var;
}

transformed parameters{
  
  real <lower = 0> rtu_mu[J,Tr];
  real <lower = 0> rtp_mu[J,Tr];
  real <lower = 0> rlr[J,Tr]; 

  for(j in 1:J){
    for(t in 1:Tr){
      rtu_mu[j,t]= rtu_intercept[j]*( 1+ rtu_base[j]*(t^(-rtu_rate[j])-1));
      rlr[j,t]= rlr_intercept[j]*(1 + (rlr_base[j]*((t^(-rlr_rate[j]))-1)));
      rtp_mu[j,t]= rtu_mu[j,t]*rlr[j,t];
    }
  }
}

model{
  //priors
  
  //group prior
  target+= gamma_lpdf(sigma_mean|alpha,beta);
  target+= gamma_lpdf(rtu_intercept_mean|mu0,tau0);
  target+= gamma_lpdf(rtu_base_mean|alpha0,beta0);
  target+= gamma_lpdf(rtu_rate_mean|alpha1,beta1);
  target+= gamma_lpdf(sigma_var|a,b);
  target+= gamma_lpdf(rtu_intercept_var|a0,b0);
  target+= gamma_lpdf(rtu_base_var|a1,b1);
  target+= gamma_lpdf(rtu_rate_var|a2,b2);
  target+= gamma_lpdf(rlr_intercept_mean|mu1,tau1);
  target+= gamma_lpdf(rlr_intercept_var|a3,b3);
  target+= gamma_lpdf(rlr_base_mean|alpha2,beta2);
  target+= gamma_lpdf(rlr_rate_mean|alpha3,beta3);
  target+= gamma_lpdf(rlr_base_var|a5,b5);
  target+= gamma_lpdf(rlr_rate_var|a6,b6);
  
  
  target+= gamma_lpdf(sigma|sigma_mean,sigma_var);
  target+= normal_lpdf(rtu_intercept|rtu_intercept_mean,rtu_intercept_var);
  target+= normal_lpdf(rtu_base|rtu_base_mean,rtu_base_var);
  target+= normal_lpdf(rtu_rate| rtu_rate_mean,rtu_rate_var);
  target += normal_lpdf(rlr_intercept| rlr_intercept_mean, rlr_intercept_var);
  target += normal_lpdf(rlr_base| rlr_base_mean, rlr_base_var);
  target += normal_lpdf(rlr_rate| rlr_rate_mean, rlr_rate_var);


  //likelihood
  for(j in 1:J){
    for(t in 1:Tr){
      target += normal_lpdf(rtu[j,t]|rtu_mu[j,t],sigma[j]);
      target += normal_lpdf(rtp[j,t]|rtp_mu[j,t],sigma[j]);
    }
  }
}

generated quantities{
  
  vector[J*Tr] log_lik;
  
  for(j in 1:J){
    for(t in 1:Tr){
      log_lik[((j-1)*Tr+t)]= normal_lpdf(rtu[j,t]|rtu_mu[j,t],sigma[j]) + normal_lpdf(rtp[j,t]|rtp_mu[j,t],sigma[j]);
    }
  }
}
