data{
  //data
  int Tr; //number of presentations per subject
  int J; //number of subjects
  real <lower = 0> rtu[J,Tr]; // unpredictable response time
  
  //prior params
  real <lower = 0> alpha;
  real <lower = 0> beta;
  
  real <lower = 0> mu0;
  real <lower = 0> tau0;
  real <lower = 0> alpha0;
  real <lower = 0> beta0;
  real <lower = 0> alpha1;
  real <lower = 0> beta1;
    
  real <lower = 0> a;
  real <lower = 0> b;
  real <lower = 0> a0;
  real <lower = 0> b0;
  real <lower = 0> a1;
  real <lower = 0> b1;
  real <lower = 0> a2;
  real <lower = 0> b2;
  
}
parameters{

  vector <lower = 0>[J] sigma;
  vector <lower = 0>[J] u_intercept;
  vector <lower = 0, upper = 1>[J] u_base;
  vector <lower = 0>[J] u_rate;

  real <lower = 0> sigma_mean;
  real <lower = 0> sigma_var;
  
  real <lower = 0> u_intercept_mean;
  real <lower = 0> u_intercept_var;
  
  real <lower = 0, upper = 1> u_base_mean;
  real <lower = 0> u_base_var;
  
  real <lower = 0> u_rate_mean;
  real <lower = 0> u_rate_var;
}
transformed parameters{
  
  row_vector <lower = 0>[Tr] mu[J];
  vector <lower = 0>[Tr] trial;
    
  trial=1:Tr;
    
  for(j in 1:J){
    mu[j] = u_intercept[j]*(1 + u_base[j]*(t^(-u_rate[j])-1));
  }
}
model{
  //priors
  
  //group prior
  target+= gamma_lpdf(sigma_mean|alpha,beta);
  target+= gamma_lpdf(u_intercept_mean|mu0,tau0);
  target+= gamma_lpdf(u_base_mean|alpha0,beta0);
  target+= gamma_lpdf(u_rate_mean|alpha1,beta1);
  target+= gamma_lpdf(sigma_var|a,b);
  target+= gamma_lpdf(u_intercept_var|a0,b0);
  target+= gamma_lpdf(u_base_var|a1,b1);
  target+= gamma_lpdf(u_rate_var|a2,b2);
  
  target+= gamma_lpdf(sigma|sigma_mean,sigma_var);
  target+= normal_lpdf(u_intercept|u_intercept_mean,u_intercept_var);
  target+= normal_lpdf(u_base|u_base_mean,u_base_var);
  target+= normal_lpdf(u_rate| u_rate_mean,u_rate_var);

  //likelihood
  for(j in 1:J){
    for(t in 1:Tr){
      target += normal_lpdf(rtu[j,t]|mu[j,t],sigma[j]);
    }
  }
}
