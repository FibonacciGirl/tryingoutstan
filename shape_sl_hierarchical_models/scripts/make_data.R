
source('scripts/model_predict.R')
source('scripts/DBDA2E-utilities.R')

convert.prior<-function(prior, m){
  new.prior<-list()
  
  
  new.prior<-list('intercept' = prior$intercept, 
                  'base' = prior$base,
                  'rate' = prior$rate,
                  'proportion' = prior$proportion,
                  'jump' = prior$jump,
                  'base.1' = prior$base,
                  'rate.1' = prior$rate,
                  'split' = prior$split,
                  'sigma' = prior$sigma)
  return(new.prior)
  
}




model.params<-function(fit,m){
  
  if(m == 'power.constant'){
    params<- c(as.numeric(fit$power.constant$intercept),
               as.numeric(fit$power.constant$base),
               as.numeric(fit$power.constant$rate),
               as.numeric(fit$power.constant$proportion),
               as.numeric(fit$power.constant$sigma))
  }
  
  if(m == 'power.power'){
    params<-c(as.numeric(fit$power.power$intercept),
              as.numeric(fit$power.power$base),
              as.numeric(fit$power.power$rate),
              as.numeric(fit$power.power$proportion),
              as.numeric(fit$power.power$base.1),
              as.numeric(fit$power.power$rate.1),
              as.numeric(fit$power.power$sigma))
  }
  if(m == 'piecewise.power.constant'){
    params<-c(as.numeric(fit$piecewise.power.constant$intercept),
              as.numeric(fit$piecewise.power.constant$base),
              as.numeric(fit$piecewise.power.constant$rate),
              as.numeric(fit$piecewise.power.constant$proportion),
              as.numeric(fit$piecewise.power.constant$jump),
              as.numeric(fit$piecewise.power.constant$split),
              as.numeric(fit$piecewise.power.constant$sigma))
  }
  if(m == 'power.logistic'){
    params<-c(as.numeric(fit$power.logistic$intercept),
              as.numeric(fit$power.logistic$base),
              as.numeric(fit$power.logistic$rate),
              as.numeric(fit$power.logistic$proportion),
              as.numeric(fit$power.logistic$jump),
              as.numeric(fit$power.logistic$rate.1),
              as.numeric(fit$power.logistic$split),
              as.numeric(fit$power.logistic$sigma))
  }
  return(params)
}



r.normal<- function(params){
  mu<-params[1]
  sd<-params[2]
  return(rnorm(1, mean = mu , sd = sd))
}

r.gamma<- function(params){
  a<-gammaShRaFromModeSD(params[1], params[2])$shape
  b<-gammaShRaFromModeSD(params[1], params[2])$rate
  return(rgamma(1, shape = a , rate = b))
}

r.beta<-function (params){
  a<-params[1]
  b<-params[2]
  return(rbeta(1, a, b))
}

r.exp<-function(params){
  return(rexp(1, params[1]))
}

r.unif<-function(params){
  return(runif(1, params[1], params[2]))
}

r.pois<-function(params){
  return(rpois(1, params[1]))
}

r.cat<-function(params){
  return(rmultinom(1, 1, rep(params[1],params[2])))
}


r.prior.fn<- function(fn, params){
  return(fn(params))
}



r.prior<-function(type, priors){
  
  if(type == 'intercept'){
    fn <- priors$intercept$fn
    params<- c(as.numeric(priors$intercept$par1), as.numeric(priors$intercept$par2))
  }
  if(type == 'base'){
    fn <- priors$base$fn
    params<- c(as.numeric(priors$base$par1), as.numeric(priors$base$par2))
  }
  if(type == 'base.1'){

    fn <- priors$base.1$fn
    params<- c(as.numeric(priors$base.1$par1), as.numeric(priors$base.1$par2))

  }
  if(type == 'rate' ){
    fn <- priors$rate$fn
    params<- c(as.numeric(priors$rate$par1), as.numeric(priors$rate$par2))
  }
  if(type == 'rate.1'){
    fn <- priors$rate.1$fn
    params<- c(as.numeric(priors$rate.1$par1), as.numeric(priors$rate.1$par2))
  }
  if(type == 'sigma'){
    fn <- priors$sigma$fn
    params<- c(as.numeric(priors$sigma$par1), as.numeric(priors$sigma$par2))
  }
  if(type == 'split'){
    fn<- priors$split$fn
    params<- c(as.numeric(priors$split$par1), as.numeric(priors$split$par2))
  }
  if(type == 'jump'){
    fn <- priors$jump$fn
    params<- c(as.numeric(priors$jump$par1), as.numeric(priors$jump$par2))
  }
  if(type == 'proportion'){
    fn <- priors$proportion$fn
    params<- c(as.numeric(priors$proportion$par1), as.numeric(priors$proportion$par2))
  }


  


  
  if(fn == 'normal'){
    return(r.prior.fn(r.normal,params))
  }
  if(fn == 'gamma'){
    return(r.prior.fn(r.gamma,params))
  }
  if(fn == 'beta'){
    return(r.prior.fn(r.beta,params))
  }
  if(fn == 'exp'){
    return(r.prior.fn(r.exp,params))
  }
  if(fn == 'unif'){
    return(r.prior.fn(r.unif,params))
  }
  if(fn == 'pois'){
    return(r.prior.fn(r.pois,params))
  }
  if(fn == 'cat'){
    return(r.prior.fn(r.cat,params))
  }
}


generate.parameter.values<-function(priors){
  
  
  recovery.parameter.values<- list(
    power.constant = list(intercept = r.prior('intercept',priors),
                          base= r.prior('base',priors), 
                          rate= - r.prior('rate',priors),
                          proportion = r.prior('proportion',priors),
                          sigma = r.prior('sigma',priors)),
   
    piecewise.power.constant = list(intercept = r.prior('intercept',priors),
                                    base= r.prior('base',priors), 
                                    rate= -r.prior('rate',priors),
                                    proportion = r.prior('proportion',priors),
                                    jump = r.prior('jump',priors),
                                    split = r.prior('split',priors),
                                    sigma=r.prior('sigma',priors)),
    
    power.power = list(intercept = r.prior('intercept',priors),
                       base= r.prior('base',priors), 
                       rate= -r.prior('rate',priors),
                       proportion = r.prior('proportion',priors), 
                       base.1= r.prior('base.1',priors), 
                       rate.1= -r.prior('rate.1',priors),
                       sigma = r.prior('sigma',priors)),
    
    power.logistic = list(intercept = r.prior('intercept',priors), 
                          base= r.prior('base',priors), 
                          rate= - r.prior('rate',priors),
                          proportion = r.prior('proportion',priors), 
                          jump = r.prior('jump',priors), 
                          rate.1= -r.prior('rate.1',priors),
                          split = r.prior('split',priors),
                          sigma=r.prior('sigma',priors)) )
  return(recovery.parameter.values)
}


add.noise<-function(recover,sigma){
  n<-length(recover)
  

  
  recover.plus.noise<-numeric()
  for(i in 1:n){
    noise<-rnorm(1, 0, sigma)
    val<-recover[i]
    while(abs(val)<= abs(noise)){
      noise<-rnorm(1,0,sigma)
    }
    recover.plus.noise[i]<-val+noise
  }
  return(recover.plus.noise)
}


create.recovery.data<-function(par, t, models= c('power.constant',
                                                'power.logistic',
                                                'power.power',
                                                'piecewise.power.constant')){

  
  
  for(m in models){
    
    
    if(m == 'power.constant'){
      params<-model.params(par, 'power.constant')
      rtu<-mapply(u.power.model.predict,t, MoreArgs = list(params[1],params[2],params[3]))
      rl<-mapply(rl.constant.model.predict,t,MoreArgs= list(params[4]))
      rtp<-rl*rtu
      power.constant.recover.u<-add.noise(rtu,params[5])
      power.constant.recover.p<-add.noise(rtp,params[5])
    }

    if(m == 'power.logistic'){

      
      params<-model.params(par, 'power.logistic')
      rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
      rl<-mapply(rl.logistic.model.predict,t,MoreArgs = list(params[4],params[5],params[6],params[7]))
      rtp<-rl*rtu
      power.logistic.recover.u<-add.noise(rtu,params[8])
      power.logistic.recover.p<-add.noise(rtp,params[8])
    }

    if(m == 'power.power'){
      params<-model.params(par, 'power.power')
      rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
      rl<-mapply(rl.power.model.predict,t,MoreArgs = list(params[4],params[5],params[6]))

      rtp<-rl*rtu
      power.power.recover.u<-add.noise(rtu,params[7])
      power.power.recover.p<-add.noise(rtp,params[7])
    }
    if(m == 'piecewise.power.constant'){
      rtp<-numeric(length(t))
      params<-model.params(par, 'piecewise.power.constant')
      t1<-head(t,params[6])
      t2<-tail(t,-params[6])
      
      rtu<-as.numeric(mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3])))
      rl1<-as.numeric(mapply(rl.constant.model.predict,t1,MoreArgs = list(params[4])))
      rl2<-as.numeric(mapply(rl.constant.model.predict,t2,MoreArgs = list(params[5])))

      rtp[t1]<-rl1[t1]*rtu[t1]

      
      rtp<-c(rtp, rl2*rtu[t2])

      piecewise.power.constant.recover.u<-add.noise(rtu,params[7])
      piecewise.power.constant.recover.p<-add.noise(rtp,params[7])
      
   
    }
  }
  
  
  recovery.data<-expand.grid(t=t, model=models)
  
  rtu<-c()
  rtp<-c()
  
  for(m in models){
    if(m == 'power.constant'){
      rtu<-c(rtu,power.constant.recover.u)
      rtp<-c(rtp,power.constant.recover.p)
    }
    
    if(m == 'power.logistic'){
      rtu<-c(rtu,power.logistic.recover.u)
      rtp<-c(rtp,power.logistic.recover.p)
    }
    
    if(m == 'power.power'){
      rtu<-c(rtu,power.power.recover.u)
      rtp<-c(rtp,power.power.recover.p)
    }
    if(m == 'piecewise.power.constant'){
      rtu<-c(rtu,piecewise.power.constant.recover.u)
      rtp<-c(rtp,piecewise.power.constant.recover.p)
    }
  }
  
  recovery.data$predictable<-rtp
  recovery.data$unpredictable<-rtu
  
  return(recovery.data)  
}



fake.subject.data <- function(prior, model, t, n.subjects){
  
  fake.prior<-convert.prior(prior, model)
  
  params = list()
  fake.data = data.frame(subject=numeric(),t=numeric(),model=character(),predictable= numeric(),unpredictable=numeric())
  
  for(i in 1:n.subjects){
    
    params[[i]]<-generate.parameter.values(fake.prior)
    fake.subject <- cbind(subject = rep(i, length(t)),create.recovery.data(params[[i]], t, models = model))
    fake.data <- rbind(fake.data, fake.subject)
  }
  return(list(fake.data = fake.data,params = params))
}

library(tidyr)



get.plot.data<-function(fake.data, prior.mean, prior.var,model){
  
  fake.pred<-data.frame(j= fake.data$subject,t= fake.data$t, rt= fake.data$predictable, p = rep(1,length(fake.data$t)))
  fake.unpred<-data.frame(j = fake.data$subject, t= fake.data$t, rt = fake.data$unpredictable, p = rep(0,length(fake.data$t)))
  fake.data<-rbind(fake.pred,fake.unpred)
  
  
  print(fake.data)
  
  Tr = max(fake.data$t)
  J  = length(unique(fake.data$j))
  P  = length(unique(fake.data$p))
  N  = length(fake.data$rt)
  rt = fake.data$rt
  jj = fake.data$j
  tt = fake.data$t
  pp = fake.data$p
  
  
  
  sigma.mean= prior.params(prior.mean$sigma)
  rtu_intercept_mean= prior.params(prior.mean$intercept)
  rtu_base_mean= prior.params(prior.mean$base)
  rtu_rate_mean= prior.params(prior.mean$rate)
  rlr_constnat_intercept_mean= prior.params(prior.mean$proportion)
  rlr_logistic_intercept_mean= prior.params(prior.mean$proportion)
  rlr_logistic_jump_mean= prior.params(prior.mean$jump)
  rlr_logistic_rate_mean= prior.params(prior.mean$rate)
  rlr_logistic_split_mean= prior.params(prior.mean$split)
  
  
  
  
  alpha= sigma.mean[1]
  beta= sigma.mean[2]
  
  mu0 = rtu_intercept_mean[1]
  tau0= rtu_intercept_mean[2]
  mu1 = rtu_base_mean[1]
  tau1= rtu_base_mean[2]
  mu2 = rtu_rate_mean[1]
  tau2= rtu_rate_mean[2]
  mu3 = rlr_constnat_intercept_mean[1]
  tau3=rlr_constnat_intercept_mean[2]
  mu4 = rlr_logistic_intercept_mean[1]
  tau4= rlr_logistic_intercept_mean[2]
  mu5 = rlr_logistic_jump_mean[1]
  tau5=rlr_logistic_jump_mean[2]
  mu6 = rlr_logistic_rate_mean[1]
  tau6= rlr_logistic_rate_mean[2]
  mu7 = rlr_logistic_split_mean[1]
  tau7= rlr_logistic_split_mean[2]
  
  
  sigma.var= prior.params(prior.var$sigma)
  rtu_intercept_var= prior.params(prior.var$intercept)
  rtu_base_var= prior.params(prior.var$base)
  rtu_rate_var= prior.params(prior.var$rate)
  rlr_constnat_intercept_var= prior.params(prior.var$proportion)
  rlr_logistic_intercept_var= prior.params(prior.var$proportion)
  rlr_logistic_jump_var= prior.params(prior.var$jump)
  rlr_logistic_rate_var= prior.params(prior.var$rate)
  rlr_logistic_split_var= prior.params(prior.var$split)
  
  
  a= sigma.var[1]
  b= sigma.var[2]
  
  alpha0 = rtu_intercept_var[1]
  beta0= rtu_intercept_var[2]
  alpha1 = rtu_base_var[1]
  beta1= rtu_base_var[2]
  alpha2 = rtu_rate_var[1]
  beta2= rtu_rate_var[2]
  alpha3 = rlr_constnat_intercept_var[1]
  beta3=rlr_constnat_intercept_var[2]
  alpha4 = rlr_logistic_intercept_var[1]
  beta4= rlr_logistic_intercept_var[2]
  alpha5 = rlr_logistic_jump_var[1]
  beta5=rlr_logistic_jump_var[2]
  alpha6 = rlr_logistic_rate_var[1]
  beta6= rlr_logistic_rate_var[2]
  alpha7 = rlr_logistic_split_var[1]
  beta7= rlr_logistic_split_var[2]
  
  p1=.5
  p2=.5
  p_prior=c(p1,p2)
  
  
  nchains = 1

  if(model == 'power.logistic'){
    model.data<- list(rep(list(
      Tr = Tr,
      J=J,
      P=P,
      N=N,
      rt = rt,
      jj=jj,
      tt=tt,
      pp=pp,
      p_prior=p_prior,
      alpha = alpha,
      beta =beta,
      a = a,
      b = b,
      mu0 = mu0,
      tau0 = tau0,
      mu1 = mu1,
      tau1 = tau1,
      mu2 = mu2,
      #tau2 = tau2,
      mu3 = mu3,
      tau3 = tau3,
      mu4 = mu4,
      tau4 = tau4,
      mu5 = mu5,
      #tau5 = tau5,
      mu6 = mu6,
      #tau6 = tau6,
      mu7 = mu7,
      tau7 = tau7,
      
      
      alpha0 = alpha0,
      beta0 = beta0,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha2 = alpha2,
      beta2 = beta2,
      alpha3 = alpha3,
      beta3 = beta3,
      alpha4 = alpha4,
      beta4 = beta4,
      alpha5 = alpha5,
      beta5=beta5,
      alpha6 = alpha6,
      beta6=beta6,
      alpha7 = alpha7,
      beta7=beta7
  
    ), nchains))
  }
  if(model == 'power.power'){
    model.data<-list(rep(list(
      Tr = Tr,
      J=J,
      P=P,
      N=N,
      rt = rt,
      jj=jj,
      tt=tt,
      pp=pp,
      p_prior=p_prior,
      alpha = alpha,
      beta = beta,
      a = a,
      b = b,
      mu0 = mu0,
      tau0 = tau0,
      mu1 = mu1,
      tau1 = tau1,
      mu2 = mu2,
      #tau2 = tau2,
      mu3 = mu3,
      tau3 = tau3,
      mu4 = mu4,
      tau4 = tau4,
      mu5 = mu5,
      #tau5 = tau5,
      mu6 = mu6,
      #tau6 = tau6,
  
      
      
      alpha0 = alpha0,
      beta0 = beta0,
      alpha1 = alpha1,
      beta1 = beta1,
      alpha2 = alpha2,
      beta2 = beta2,
      alpha3 = alpha3,
      beta3 = beta3,
      alpha4 = alpha4,
      beta4 = beta4,
      alpha5 = alpha5,
      beta5=beta5,
      alpha6 = alpha6,
      beta6=beta6
  
    ), nchains))
  }
  
  return(model.data)

}

prior.params<-function(prior){

  fn<-prior$fn
  params<-c(prior$par1,prior$par2)
  
  if(fn == 'normal'){
    return(params)
  }
  if(fn == 'gamma'){
    a<-gammaShRaFromModeSD(params[1], params[2])$shape
    b<-gammaShRaFromModeSD(params[1], params[2])$rate
    params<-c(a,b)
    return(params)
  }
  if(fn=='beta'){
    a<-betaABfromModeKappa(params[1], params[2])$a
    b<-betaABfromModeKappa(params[1],params[2])$b
    params<-c(a,b)
    return(params)
  }
  
  if(fn =='unif'){
    return(params)
  }
  
  if(fn == 'exp'){
    return(params)
  }
}



