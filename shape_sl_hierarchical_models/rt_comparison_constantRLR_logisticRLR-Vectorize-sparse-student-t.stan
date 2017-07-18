data{
  
  //SUPPORT FOR POSTERIOR  
  int Tr;                       //total number of stimuli presentations
  int J;                        //number of subjects
  int P;                        //predictable indicator
  int N;                        //number of data points
  
  int <lower = 0,upper=Tr> tt[N];   // trial index
  int <lower = 0,upper=J> jj[N];   // subject index
  int <lower = 0,upper=P> pp[N];   // predictable index
  
  real <lower = 0,upper=2000> rt[N];   // response time



  //HYPER PARAMETERS
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
  
  real <lower = 0> mu4;         //rlr logistic proportion prior
  real <lower = 0> tau4;
  real <lower = 0> alpha4;
  real <lower = 0> beta4;
  
  real <lower = 0> mu5;         //rlr logistic proportion 2 prior
  real <lower = 0> tau5;
  real <lower = 0> alpha5;
  real <lower = 0> beta5;
  
  real <lower = 0> mu6;         //rlr logistic base prior
  real <lower = 0> tau6;
  real <lower = 0> alpha6;
  real <lower = 0> beta6;
  
  real <lower = 0, upper= Tr> mu7;         //rlr logistic split prior
  real <lower = 0> tau7;
  real <lower = 0> alpha7;
  real <lower = 0> beta7;
  
}
parameters{
  real <lower = 0> nu; //degrees of freedom
  
  //model selection parameter
  simplex[2] p;
  simplex[2] pi0[J];

  
//INDIVIDUAL LEVEL PARAMETERS
  vector <lower = 0>[J] sigma;
  
  //rtu
  vector <lower = 0>[J] rtu_intercept;
  vector <lower = 0, upper = 1>[J] rtu_base;
  vector <lower = 0>[J] rtu_rate;
  
  //rlr constant
  vector <lower = 0>[J] rlr_constant_intercept;
  
  //rlr logistic
  vector <lower = 0>[J] rlr_logistic_intercept;
  vector <lower = 0>[J] rlr_logistic_rate;
  vector <lower = 0, upper = 1>[J] rlr_logistic_intercept2;
  vector <lower = 0 , upper = Tr>[J] rlr_logistic_split;
  
//GROUP LEVEL PARAMETERS
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
  real <lower = 0, upper = 5>  rlr_constant_intercept_mean;
  real <lower = 0>  rlr_constant_intercept_var;
  
  //rlr logistic
  real <lower = 0,upper=5> rlr_logistic_intercept_mean;
  real <lower = 0> rlr_logistic_intercept_var;
  
  real <lower = 0, upper = 1> rlr_logistic_rate_mean;
  real <lower = 0> rlr_logistic_rate_var;
  
  real <lower = 0, upper = 1> rlr_logistic_intercept2_mean;
  real <lower = 0> rlr_logistic_intercept2_var;
  
  real <lower = 0> rlr_logistic_split_mean;
  real <lower = 0> rlr_logistic_split_var;
  
}

transformed parameters{
  
  real <lower = 0, upper = 2000> rt_mu[N]; //likelihood rtu mean
  real <lower = 0, upper = 2000> rt1_mu[N]; //likelihood rlr constant mean
  real <lower = 0, upper = 2000> rt2_mu[N]; //likelihood rlr logistic mean
  
  //log likelihood summation over indicator 
  vector[J] log_q_z1; 
  vector[J] log_q_z2;
  
  
  for(j in 1:J){
    
    
    log_q_z1[j] = log(pi0[j,1]) 
                  +  student_t_lpdf(rlr_constant_intercept[j]| nu,rlr_constant_intercept_mean, rlr_constant_intercept_var);
    log_q_z2[j] = log(pi0[j,2]) 
                  + student_t_lpdf(rlr_logistic_intercept[j]|nu,rlr_logistic_intercept_mean,rlr_logistic_intercept_var)
                  + student_t_lpdf(rlr_logistic_rate[j]|nu, rlr_logistic_rate_mean,rlr_logistic_rate_var)
                  + student_t_lpdf(rlr_logistic_intercept2[j]|nu, rlr_logistic_intercept2_mean,rlr_logistic_intercept2_var)
                  + student_t_lpdf(rlr_logistic_split[j]|nu, rlr_logistic_split_mean,rlr_logistic_split_var);
  }

  for(n in 1:N){
    rt_mu[n]= rtu_intercept[jj[n]]*( 1+ rtu_base[jj[n]]*((tt[n])^(-rtu_rate[jj[n]])-1));   
    if(pp[n]==1){
      rt1_mu[n] = rt_mu[n]*rlr_constant_intercept[jj[n]];
      rt2_mu[n] = rt_mu[n]*(rlr_logistic_intercept[jj[n]]+(rlr_logistic_intercept2[jj[n]]-rlr_logistic_intercept[jj[n]])/(1+exp(rlr_logistic_rate[jj[n]]*(tt[n] - rlr_logistic_split[jj[n]]))));

      log_q_z1[jj[n]] = log_q_z1[jj[n]] + normal_lpdf(rt[n]|rt1_mu[n],sigma[jj[n]]);

              
      log_q_z2[jj[n]] = log_q_z2[jj[n]] + normal_lpdf(rt[n]|rt2_mu[n],sigma[jj[n]]);

    }
    else{
      rt1_mu[n] = 0;
      rt2_mu[n] = 0;
    }
  }
}
model{
  //model selection prior
  target += dirichlet_lpdf(p|p_prior);
  
//GROUP PRIOR CONTRIBUTION
  
  //sigma
  target+= gamma_lpdf(sigma_mean|alpha,beta);
  target+= gamma_lpdf(sigma_var|a,b);
  
  //rtu
  target+= gamma_lpdf(rtu_intercept_mean|mu0,tau0);
  target+= uniform_lpdf(rtu_base_mean|mu1,tau1);
  target+= gamma_lpdf(rtu_rate_mean|mu2,tau2);
  

  target+= gamma_lpdf(rtu_intercept_var|alpha0, beta0);
  target+= gamma_lpdf(rtu_base_var|alpha1,beta1);
  target+= gamma_lpdf(rtu_rate_var|alpha2,beta2);
  
  //rlr logistic  
  target += gamma_lpdf(rlr_logistic_intercept_mean|mu4,tau4);
  target += gamma_lpdf(rlr_logistic_rate_mean|mu5,tau5);
  target += beta_lpdf(rlr_logistic_intercept2_mean|mu6,tau6);
  target += uniform_lpdf(rlr_logistic_split_mean|mu7,tau7);  
  
  target += gamma_lpdf(rlr_logistic_intercept_var|alpha4, beta4);
  target += gamma_lpdf(rlr_logistic_rate_var|alpha5,beta5);
  target += gamma_lpdf(rlr_logistic_intercept2_var|alpha6,beta6);
  target += gamma_lpdf(rlr_logistic_split_var|alpha7,beta7); 

  //rlr constant
  target += gamma_lpdf(rlr_constant_intercept_mean|mu3,tau3);
  
  target += gamma_lpdf(rlr_constant_intercept_var|alpha3,beta3);

  
//INDIVIDUAL PRIOR CONTRIBUTION


  //rtu
  

  

  //likelihood
  //likelihood
  for(n in 1:N){
    if(pp[n]==0){
      target += normal_lpdf(rt[n]|rt_mu[n],sigma[jj[n]]);
    }
  }
  for(j in 1:J){
      target += log_sum_exp(log_q_z1[j],log_q_z2[j]);
      
      target += dirichlet_lpdf(pi0[j]|p);
      
      target+= normal_lpdf(rtu_intercept[j]|rtu_intercept_mean,rtu_intercept_var);
      target+= normal_lpdf(rtu_base[j]|rtu_base_mean,rtu_base_var);
      target+= normal_lpdf(rtu_rate[j]| rtu_rate_mean,rtu_rate_var);
  

    }
}

