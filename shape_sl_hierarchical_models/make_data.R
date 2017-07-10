
source('model_predict.R')
source('DBDA2E-utilities.R')


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
              as.numeric(fit$power.logistic$upper.proportion),
              as.numeric(fit$power.logistic$lower.proportion),
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
  if(type=='base.1'){
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
  if(type == 'upper.proportion'){
    fn <- priors$upper.proportion$fn
    params<- c(as.numeric(priors$upper.proportion$par1), as.numeric(priors$upper.proportion$par2))
  }
  if(type == 'lower.proportion'){
    fn <- priors$lower.proportion$fn
    params<- c(as.numeric(priors$lower.proportion$par1), as.numeric(priors$lower.proportion$par2))
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
                          rate= -r.prior('rate',priors),
                          proportion = r.prior('proportion',priors),
                          sigma=r.prior('sigma',priors)),
   
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
                       sigma=r.prior('sigma',priors)),
    
    power.logistic = list(intercept = r.prior('intercept',priors), 
                          base= r.prior('base',priors), 
                          rate= - r.prior('rate',priors),
                          upper.proportion = r.prior('upper.proportion',priors), 
                          lower.proportion = r.prior('lower.proportion',priors), 
                          rate.1= -r.prior('rate',priors),
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
      params<-model.params(par, 'piecewise.power.constant')
      t1<-head(t,params[6])
      t2<-tail(t,-params[6])
      

      rtu<-mapply(u.power.model.predict,t,MoreArgs = list(params[1],params[2],params[3]))
      rl1<-mapply(rl.constant.model.predict,t1,MoreArgs = list(params[4]))
      rl2<-mapply(rl.constant.model.predict,t2,MoreArgs = list(params[5]))
      rtp[t1]<-rl1*rtu[t1]
      rtp[t2]<-rl2*rtu[t2]
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

params<-generate.parameter.values(priors)
data.recover<-create.recovery.data(params)
subset(data.recover,model=='power.power')

plot.data.fits(subject.data=subset(data.recover,model=='power.power'))


rtu
