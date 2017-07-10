data{
  //data
  int Tr;                       //number of stimuli presentations
  int J;                        //number of subjects
  real <lower = 0,upper=2000> rtu[J,Tr];   // unpredictable response time
  real <lower = 0,upper=2000> rtp[J,Tr];   // predictable response time
  
  //prior parameters 

  real <lower = 0> alpha;       //sigma prior
  real <lower = 0> beta;
  real <lower = 0> a;           
  real <lower = 0> b;
  
  real <lower = 0> mu0;         //rtu intercept prior 
  real <lower = 0> tau0;
  real <lower = 0> alpha0;
  real <lower = 0> beta0;
  
  real <lower = 0> mu1;         //rtu base prior
  real <lower = 0> tau1;
  real <lower = 0> alpha1;
  real <lower = 0> beta1;
  
  real <lower = 0> mu2;         //rtu rate prior
  real <lower = 0> tau2;
  real <lower = 0> alpha2;
  real <lower = 0> beta2;
  
  real <lower = 0> mu3;         //rlr constant intercept prior
  real <lower = 0> tau3;
  real <lower = 0> alpha3;
  real <lower = 0> beta3;
  
  real <lower = 0> mu4;         //rlr constant intercept2 prior
  real <lower = 0> tau4;
  real <lower = 0> alpha4;
  real <lower = 0> beta4;
  
}
transformed data{
  real log_unif;
  log_unif = -log(Tr);
}
parameters{

//individual level
  vector <lower = 0>[J] sigma;
  
  //rtu
  vector <lower = 0, upper = 2000>[J] rtu_intercept;
  vector <lower = 0, upper = 1>[J] rtu_base;
  vector <lower = 0>[J] rtu_rate;
  
  //rlr constant-piecewise
  vector <lower = 0>[J] rlr_constant_intercept;
  vector <lower = 0>[J] rlr_constant_intercept2;

  
//group level
  real <lower = 0> sigma_mean;
  real <lower = 0> sigma_var;
  
  //rtu
  real <lower = 0, upper=2000> rtu_intercept_mean;
  real <lower = 0> rtu_intercept_var;
  
  real <lower = 0, upper = 1> rtu_base_mean;
  real <lower = 0> rtu_base_var;
  
  real <lower = 0> rtu_rate_mean;
  real <lower = 0> rtu_rate_var;
  
  //rlr cosntant
  real <lower = 0, upper = 2000>  rlr_constant_intercept_mean;
  real <lower = 0>  rlr_constant_intercept_var;
  
  //rlr power
  real <lower = 0,upper=2000> rlr_constant_intercept2_mean;
  real <lower = 0> rlr_constant_intercept2_var;
  
}
transformed parameters{
  
  real <lower = 0> rtu_mu[J,Tr];
  real <lower = 0> rtp_mu1[J,Tr];
  real <lower = 0> rtp_mu2[J,Tr];
  vector[Tr] lp_split[J];

  for(j in 1:J){
    for(t in 1:Tr){
      rtu_mu[j,t]= rtu_intercept[j]*( 1+ rtu_base[j]*(t^(-rtu_rate[j])-1));

      rtp_mu1[j,t]= rtu_mu[j,t]*rlr_constant_intercept[j];
      rtp_mu2[j,t]= rtu_mu[j,t]*rlr_constant_intercept2[j];
      for(s in 1:Tr){
        lp_split[j,s] = lp_split[j,s] + normal_lpdf(rtp[j,t] | t < s ? rtp_mu1[j,t] : rtp_mu2[j,t] , sigma[j]);
      }
    }
  }
}
model{
  //priors
  
  //group prior
  target+= gamma_lpdf(sigma_mean|alpha,beta);
  target+= gamma_lpdf(rtu_intercept_mean|mu0,tau0);
  target+= gamma_lpdf(rtu_base_mean|mu1,tau1);
  target+= gamma_lpdf(rtu_rate_mean|mu2,tau2);
  target+= gamma_lpdf(sigma_var|a,b);
  target+= gamma_lpdf(rtu_intercept_var|alpha0,beta0);
  target+= gamma_lpdf(rtu_base_var|alpha1,beta1);
  target+= gamma_lpdf(rtu_rate_var|alpha2,beta2);
  target+= gamma_lpdf(rlr_constant_intercept_mean|mu3,tau3);
  target+= gamma_lpdf(rlr_constant_intercept_var|alpha3,beta3);
  target+= gamma_lpdf(rlr_constant_intercept2_mean|mu4,tau4);
  target+= gamma_lpdf(rlr_constant_intercept2_var|alpha4,beta4);


  
  
  target+= gamma_lpdf(sigma|sigma_mean,sigma_var);
  target+= normal_lpdf(rtu_intercept|rtu_intercept_mean,rtu_intercept_var);
  target+= normal_lpdf(rtu_base|rtu_base_mean,rtu_base_var);
  target+= normal_lpdf(rtu_rate| rtu_rate_mean,rtu_rate_var);
  target += normal_lpdf(rlr_constant_intercept| rlr_constant_intercept_mean, rlr_constant_intercept_var);
  target += normal_lpdf(rlr_constant_intercept2| rlr_constant_intercept2_mean, rlr_constant_intercept2_var);


  //likelihood
  for(j in 1:J){
    for(t in 1:Tr){
      target += normal_lpdf(rtu[j,t]|rtu_mu[j,t],sigma[j]);
      target += lp_split[t,j];
    }
  }
}
