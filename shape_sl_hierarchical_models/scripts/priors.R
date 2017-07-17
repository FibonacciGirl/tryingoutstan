
library(jsonlite)

prior<-function(par,type){

  if(type == 'intercept'){
    fn <- priors$intercept$fn[[1]]
    params<- c(as.numeric(priors$intercept$par1), as.numeric(priors$intercept$par2))
  }
  if(type == 'proportion'){
    fn <- priors$proportion$fn
    params<- c(as.numeric(priors$proportion$par1), as.numeric(priors$proportion$par2))
  }
  if(type == 'rate'){
    fn <- priors$rate$fn[[1]]
    params<- c(as.numeric(priors$rate$par1), as.numeric(priors$rate$par2))
  }
  if(type == 'sigma'){
    fn <- priors$sigma$fn[[1]]
    params<- c(as.numeric(priors$sigma$par1), as.numeric(priors$sigma$par2))
  }
  if(type == 'jump.proportion'){
    fn <- priors$jump.proportion$fn
    params<- c(as.numeric(priors$jump.proportion$par1), as.numeric(priors$jump.proportion$par2))
  }
  if(type == 'split'){
    fn <- priors$split$fn[[1]]
    params<- c(as.numeric(priors$split$par1), as.numeric(priors$split$par2))
  }
  
  if(fn == 'normal'){
    return(prior.fn(normal,par,params))
  }
  if(fn == 'gamma'){
    return(prior.fn(gamma,par,params))
  }
  if(fn == 'beta'){
    return(prior.fn(beta,par,params))
  }
  if(fn == 'unif'){
    return(prior.fn(unif,par,params))
  }
  if(fn == 'exp'){
    return(prior.fn(exp,par,params))
  }
}


prior.fn<- function(fn, par, params){
  return(fn(par,params))
}

normal<- function(par, params){
  mu<-params[1]
  sd<-params[2]
  return(dnorm(par, mean = mu , sd = sd, log=T))
}

gamma<- function(par, params){
  a<-gammaShRaFromModeSD(params[1], params[2])$shape
  b<-gammaShRaFromModeSD(params[1], params[2])$rate
  return(dgamma(par, shape = a , rate = b,log=T))
}

beta<-function (par ,params){
  a<-betaABfromModeKappa(params[1], params[2])$a
  b<-betaABfromModeKappa(params[1],params[2])$b
  return(dbeta(par, a, b,log=T))
}

unif<-function(par, params){
  return(dunif(par, params[1], params[2], log=T))
}

unif<-function(par, params){
  return(dexp(par, params[1], log=T))
}




