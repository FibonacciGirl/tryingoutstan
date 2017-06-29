#makefake data



linear.params<-array(dim=c(n.subjects,3))
constant.params<-array(dim=c(n.subjects,2))


t<-1:72
n.trials<-length(t)
n.subjects<-100
#theta<-rbinom(n.subjects,1,.2)
theta<-rbeta(n.subjects,.5,.5)
sigma<-rgamma(n.subjects,9,.4)
theta
rtu<-array(dim=c(n.subjects,n.trials))
rtp<-array(dim=c(n.subjects,n.trials))



for(i in 1:n.subjects){
  linear.params[i,]<-c(rnorm(1,400,50), rnorm(1,800,50), rnorm(1,-15,2.5))
  constant.params[i,]<-c(rnorm(1,400,50), rnorm(1,-400,50))
  
  rtu[i,]<-constant.model(t,constant.params[i,])$rtu + rnorm(n.trials,0,sigma[i])
  rtp[i,]<-theta[i]*linear.model(t,linear.params[i,])$rtp + (1-theta[i])*constant.model(t,constant.params[i,])$rtp + rnorm(n.trials,0,sigma[i])
}

gammaShRaFromModeSD(10,2)

J=n.subjects
Tr=n.trials

model_data<-list(J=J,
                 Tr=Tr,
                 rtu=rtu,
                 rtp=rtp)


linear.model<- function(t,params){
  a<-params[1]; b<-params[2]; c<-params[3]
  
  rtu<-mapply(function(x,a){
    rtu=a
    return(rtu)
  },t,MoreArgs= list(a))
  rtp<-mapply(function(x,b,c){
    rl=b + c*x
    return(rl)
  },t,MoreArgs= list(b,c))
  

  
  rt_data<-data.frame(rtu=rtu,rtp=rtp)
  
  return(rt_data)
}

constant.model<- function(t,params){
  a<-params[1]; b<-params[2];
  
  rtu<-mapply(function(x,a){
    rtu=a
    return(rtu)
  },t,MoreArgs= list(a))
  rtp<-mapply(function(x,b){
    rl=b 
    return(rl)
  },t,MoreArgs= list(b))
  

  
  rt_data<-data.frame(rtu=rtu,rtp=rtp)
  
  return(rt_data)
}



