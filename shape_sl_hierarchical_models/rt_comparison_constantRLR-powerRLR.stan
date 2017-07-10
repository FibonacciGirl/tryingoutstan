data{
  //data
  int Tr;                       //number of stimuli presentations
  int J;                        //number of subjects
  real <lower = 0,upper=2000> rtu[J,Tr];   // unpredictable response time
  real <lower = 0,upper=2000> rtp[J,Tr];   // predictable response time
  
  //prior parameters 
  real<lower = 0> p1;                 //model selection prior
  real<lower = 0> p2;
  
  
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
  
  real <lower = 0> mu4;         //rlr power intercept prior
  real <lower = 0> tau4;
  real <lower = 0> alpha4;
  real <lower = 0> beta4;
  
  real <lower = 0> mu5;         //rlr power base prior
  real <lower = 0> tau5;
  real <lower = 0> alpha5;
  real <lower = 0> beta5;
  
  real <lower = 0> mu6;         //rlr power rate prior
  real <lower = 0> tau6;
  real <lower = 0> alpha6;
  real <lower = 0> beta6;
  
}
parameters{
//model selection parameter
  real <lower = 0, upper = 1> pi0;
  
//individual level
  vector <lower = 0>[J] sigma;
  
  //rtu
  vector <lower = 0>[J] rtu_intercept;
  vector <lower = 0, upper = 1>[J] rtu_base;
  vector <lower = 0>[J] rtu_rate;
  
  //rlr constant
  vector <lower = 0>[J] rlr_constant_intercept;
  
  //rlr power
  vector <lower = 0>[J] rlr_power_intercept;
  vector <lower = 0, upper = 1>[J] rlr_power_base;
  vector <lower = 0>[J] rlr_power_rate;
  
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
  real <lower = 0,upper=2000> rlr_power_intercept_mean;
  real <lower = 0> rlr_power_intercept_var;
  
  real <lower = 0, upper = 1> rlr_power_base_mean;
  real <lower = 0> rlr_power_base_var;
  
  real <lower = 0> rlr_power_rate_mean;
  real <lower = 0> rlr_power_rate_var;
  
}
transformed parameters{
  
  real <lower = 0, upper = 2000> rtu_mu[J,Tr];
  real <lower = 0, upper = 2000> rtp_mu_constant[J,Tr];
  real <lower = 0, upper = 2000> rtp_mu_power[J,Tr];  
  real log_q_z1[J,Tr];
  real log_q_z2[J,Tr];


  for(j in 1:J){
    for(t in 1:Tr){
      
      rtu_mu[j,t]= rtu_intercept[j]*( 1+ rtu_base[j]*((t)^(-rtu_rate[j])-1));
      rtp_mu_constant[j,t]= rtu_mu[j,t]*rlr_constant_intercept[j];
      rtp_mu_power[j,t]= rtu_mu[j,t]*rlr_power_intercept[j]*( 1+ rlr_power_base[j]*((t)^(-rlr_power_rate[j])-1));
    

      log_q_z1[j,t]= log(pi0) + normal_lpdf(rtp[j,t]|rtp_mu_constant[j,t],sigma[j])
      + normal_lpdf(rlr_constant_intercept[j]| rlr_constant_intercept_mean, rlr_constant_intercept_var);

      
      log_q_z2[j,t]= log(1-pi0) + normal_lpdf(rtp[j,t]|rtp_mu_power[j,t],sigma[j])
      + normal_lpdf(rlr_power_intercept[j]|rlr_power_intercept_mean,rlr_power_intercept_var)
      + normal_lpdf(rlr_power_base[j]|rlr_power_base_mean,rlr_power_base_var)
      + normal_lpdf(rlr_power_rate[j]| rlr_power_rate_mean,rlr_power_rate_var);
      
      

    }
    
  }
}
model{

//model selection prior
  target += beta_lpdf(pi0|p1,p2);
  
//group priors
  
  target+= gamma_lpdf(sigma_mean|alpha,beta);
  target+= gamma_lpdf(sigma_var|a,b);
  
  //rtu
  target+= gamma_lpdf(rtu_intercept_mean|mu0,tau0);
  target+= gamma_lpdf(rtu_base_mean|mu1,tau1);
  target+= gamma_lpdf(rtu_rate_mean|mu2,tau2);
  

  target+= gamma_lpdf(rtu_intercept_var|alpha0, beta0);
  target+= gamma_lpdf(rtu_base_var|alpha1,beta1);
  target+= gamma_lpdf(rtu_rate_var|alpha2,beta2);

  
//indidvidual priors
  target+= gamma_lpdf(sigma|sigma_mean,sigma_var);
  
  //rtu
  target+= normal_lpdf(rtu_intercept|rtu_intercept_mean,rtu_intercept_var);
  target+= normal_lpdf(rtu_base|rtu_base_mean,rtu_base_var);
  target+= normal_lpdf(rtu_rate| rtu_rate_mean,rtu_rate_var);
  
  
  target += gamma_lpdf(rlr_power_intercept_mean|mu4,tau4);
  target += gamma_lpdf(rlr_power_base_mean|mu5,tau5);
  target += gamma_lpdf(rlr_power_rate_mean|mu6,tau6);
  target += gamma_lpdf(rlr_power_intercept_var|alpha4, beta4);
  target += gamma_lpdf(rlr_power_base_var|alpha5,beta5);
  target += gamma_lpdf(rlr_power_rate_var|alpha6,beta6);

  
  target += gamma_lpdf(rlr_constant_intercept_mean|mu3,tau3);
  target += gamma_lpdf(rlr_constant_intercept_var|alpha3,beta3);
  

  //likelihood
  for(j in 1:J){
    for(t in 1:Tr){

      
      target += normal_lpdf(rtu[j,t]|rtu_mu[j,t],sigma[j]);

      target += log_sum_exp(log_q_z1[j,t],log_q_z2[j,t]);
      
    }

  }
}

generated quantities{
  
  vector[J*Tr] log_lik;
  
  for(j in 1:J){
    for(t in 1:Tr){
      log_lik[((j-1)*Tr+t)]= normal_lpdf(rtu[j,t]|rtu_mu[j,t],sigma[j])+log_sum_exp(log_q_z1[j,t],log_q_z2[j,t]);
    }
  }
}

