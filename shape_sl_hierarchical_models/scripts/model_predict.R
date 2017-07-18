##Model predit##

u.power.model.predict<-function(t,a,b,c){
  rtu<- a * ( 1 + b*((t^c)-1))
  return(c(rtu=rtu))
}

rl.constant.model.predict<-function(t,d){
  rl<-d
  return(c(rl=rl))
}

rl.logistic.model.predict<-function(t,d,e,f,g){
  rl<-d+((e-d)/(1+exp(f*(t-g))))
  return(c(rl=rl))
}

rl.power.model.predict<-function(t,d,e,f){
  rl<- d*(1 + (e*((t^f)-1)))
  return(c(rl=rl))
}


# 
# source('scripts/Likelihood 2/create recovery data.R')
# priors<- list(proportion = list( fn = 'unif',par1= .75, par2= 1.25),
#               upper.proportion = list( fn = 'unif',par1= .75, par2= 1.25),
#               lower.proportion = list( fn = 'unif',par1= 0, par2= .5),
#               intercept = list(fn = 'gamma',par1 = 1000, par2 = 50),
#               base = list (fn = 'unif',par1 = .3 , par2 = .8 ),
#               rate = list (fn = 'unif',par1 = .05 , par2 =.5 ),
#               base.1 = list (fn = 'unif',par1 = .3 , par2 =.8 ),
#               rate.1 = list (fn = 'unif',par1 = .05 , par2 =.3 ),
#               jump = list (fn = 'unif',par1 = 0 , par2 =.75 ),
#               split = list(fn = 'pois', par1 = 36 , par2 = 0),
#               sigma = list (fn = 'gamma', par1 = 0.01, par2 =200 )
# )
# 
# 
# params<-generate.parameter.values(priors)
# params<-as.list(head(model.params(params,'exponential.power'),-1))
# 
# t<-1:72
# test.p.c<-mapply(power.constant.model.predict,t,MoreArgs = params )
# test.e.c<-mapply(exponential.constant.model.predict,t,MoreArgs = params)
# test.p.p<-mapply(power.power.model.predict,t,MoreArgs = params)
# test.p.e<-mapply(power.exponential.model.predict,t,MoreArgs = params)
# test.e.p<-mapply(exponential.power.model.predict,t,MoreArgs = params)
# test.e.e<-mapply(exponential.exponential.model.predict,t,MoreArgs = params)
# test.pp.c<-mapply(piecewise.power.constant.model.predict,t,MoreArgs = params)
# test.pe.c<-mapply(piecewise.exponential.constant.model.predict,t,MoreArgs = params)
# test.p.l<-mapply(power.logistic.model.predict,t,MoreArgs = params)
# test.e.l<-mapply(exponential.logistic.model.predict,t,MoreArgs = params)
# 
# library(ggplot2)
# library(gridExtra)
# plot.model.predict(t,test.p.p)
# params
# 


plot.model.predict<-function(t,test){
  plot.data<-rbind(data.frame(t=t,rt=test[1,],type = rep('unpredictable',length(t))),data.frame(t=t,rt=test[2,],type=rep('predictable',length(t))))
  rl.data<-data.frame(t=t,rt=test[3,])
  
  p<-ggplot()+
    geom_line(data=plot.data,aes(x=t,y=rt,col=type))+
    ggtitle('Response time vs. trial')+
    theme_minimal()
  
  l<-ggplot()+
    geom_line(data=rl.data,aes(x=t,y=rt))+
    ggtitle('Relative learning rate')+
    theme_minimal()
  
  i<-grid.arrange(grobs=list(p,l))
  return(i)
}
