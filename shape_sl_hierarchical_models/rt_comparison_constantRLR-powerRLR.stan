data{
  //data
  int Tr;                       //number of stimuli presentations
  int J;                        //number of subjects
  int P;                        //predictable indicator
  int N;                        // number of data points
  
  int <lower = 0,upper=Tr> tt[N];   // response time t index
  int <lower = 0,upper=J> jj[N];   // response time j index
  int <lower = 0,upper=P> pp[N];   // response time p index

  
  real <lower = 0,upper=2000> rt[N];   // response time

  
  //prior parameters 
  simplex[2] p_prior;                 //model selection prior
  
  
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
  simplex[2] p;
  simplex[2] pi0[J];

  
  //individual level

  vector <lower = 0,upper = 300>[J] sigma;
  
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
  real <lower = 0, upper = 300> sigma_mean;
  real <lower = 0> sigma_var;
  
  //rtu
  real <lower = 0, upper=2000> rtu_intercept_mean;
  real <lower = 0> rtu_intercept_var;
  
  real <lower = 0, upper = 1> rtu_base_mean;
  real <lower = 0> rtu_base_var;
  
  real <lower = 0> rtu_rate_mean;
  real <lower = 0> rtu_rate_var;
  
  //rlr cosntant
  real <lower = 0, upper = 5>  rlr_constant_intercept_mean;
  real <lower = 0>  rlr_constant_intercept_var;
  
  //rlr power
  real <lower = 0,upper=5> rlr_power_intercept_mean;
  real <lower = 0> rlr_power_intercept_var;
  
  real <lower = 0, upper = 1> rlr_power_base_mean;
  real <lower = 0> rlr_power_base_var;
  
  real <lower = 0> rlr_power_rate_mean;
  real <lower = 0> rlr_power_rate_var;
  
}
transformed parameters{
  
  real <lower = 0, upper = 2000> rt_mu[N];
  real <lower = 0, upper = 2000> rt1_mu[N];
  real <lower = 0, upper = 2000> rt2_mu[N];
  //real log_q_z1[N];
  //real log_q_z2[N];
  
  vector[J] log_q_z1;
  vector[J] log_q_z2;
  
  
    for(j in 1:J){
      log_q_z1[j] = log(pi0[j,1])
                    + normal_lpdf(rlr_constant_intercept[j]| rlr_constant_intercept_mean, rlr_constant_intercept_var);
      log_q_z2[j] = log(pi0[j,2])
                    + normal_lpdf(rlr_power_intercept[j]|rlr_power_intercept_mean,rlr_power_intercept_var)
                    + normal_lpdf(rlr_power_base[j]|rlr_power_base_mean,rlr_power_base_var)
                    + normal_lpdf(rlr_power_rate[j]| rlr_power_rate_mean,rlr_power_rate_var);
  
    }

  
  for(n in 1:N){
 
          
          
          rt_mu[n]= rtu_intercept[jj[n]]*( 1+ rtu_base[jj[n]]*((tt[n])^(-rtu_rate[jj[n]])-1));   


          rt1_mu[n] = rt_mu[n]*rlr_constant_intercept[jj[n]];
          rt2_mu[n] = rt_mu[n]*(rlr_power_intercept[jj[n]]*( 1+ rlr_power_base[jj[n]]*((tt[n])^(- rlr_power_rate[jj[n]])-1)));
              
            if(pp[n]==1){

              log_q_z1[jj[n]] = log_q_z1[jj[n]] + normal_lpdf(rt[n]|rt1_mu[n],sigma[jj[n]]);

              
              log_q_z2[jj[n]] = log_q_z2[jj[n]] + normal_lpdf(rt[n]|rt2_mu[n],sigma[jj[n]]);

            }
 
  }
  
}
model{
  
  //model selection prior
  target += dirichlet_lpdf(p|p_prior);
  
  //group priors
  
  target+= gamma_lpdf(sigma_mean|alpha,beta);
  target+= gamma_lpdf(sigma_var|a,b);
  
  //rtu
  target+= gamma_lpdf(rtu_intercept_mean|mu0,tau0);
  target+= uniform_lpdf(rtu_base_mean|mu1,tau1);
  target+= gamma_lpdf(rtu_rate_mean|mu2,tau2);
  
  
  target+= gamma_lpdf(rtu_intercept_var|alpha0, beta0);
  target+= gamma_lpdf(rtu_base_var|alpha1,beta1);
  target+= gamma_lpdf(rtu_rate_var|alpha2,beta2);
  
  
  //indidvidual priors
  target+= normal_lpdf(sigma|sigma_mean,sigma_var);
  
  //rtu

  
  

  target += gamma_lpdf(rlr_power_intercept_mean|mu4,tau4);
  target += uniform_lpdf(rlr_power_base_mean|mu5,tau5);
  target += gamma_lpdf(rlr_power_rate_mean|mu6,tau6);
  target += gamma_lpdf(rlr_power_intercept_var|alpha4, beta4);
  target += uniform_lpdf(rlr_power_base_var|alpha5,beta5);
  target += gamma_lpdf(rlr_power_rate_var|alpha6,beta6);
  
  
  target += gamma_lpdf(rlr_constant_intercept_mean|mu3,tau3);
  target += gamma_lpdf(rlr_constant_intercept_var|alpha3,beta3);
  
  
  //likelihood
  for(n in 1:N){
    if(pp[n]==0){
      target += normal_lpdf(rt[n]|rt_mu[n],sigma[jj[n]]);
    }
  }
  for(j in 1:J){
      target += dirichlet_lpdf(pi0[j]|p);
      
      target+= normal_lpdf(rtu_intercept[j]| rtu_intercept_mean,rtu_intercept_var);
      target+= normal_lpdf(rtu_base[j]| rtu_base_mean,rtu_base_var);
      target+= normal_lpdf(rtu_rate[j]| rtu_rate_mean,rtu_rate_var);
    
      target += log_sum_exp(log_q_z1[j],log_q_z2[j]);
      

    }

}


