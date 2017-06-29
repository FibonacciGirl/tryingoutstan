install.packages('tensorflow')
install.packages('devtools')
require(devtools)
tensorflow::install_tensorflow()
devtools::install_github('goldingn/greta')

library(devtools)
options(devtools.install.args = "--no-multiarch")   
require(greta)



linear.model<- function(t,params){
  a<-params[1]; b<-params[2]; c<-params[3]
  rtu<-a
  rl<-(b+c*t)
  rtp<-rtu*rl
  
  return(c(rtu=rtu,rtp=rtp,rl=rl))
}
constant.model<-function(t,params){
  a<-params[1]; b<-params[2]
  rtu<-a
  rl<-b
  rtp<-rtu*rl
  
  return(c(rtu=rtu,rtp=rtp))
}


make.fake.data<-function(t,n.subjects){

  all.subject.data<-list()
  n.trials<-length(t)
  
  for(i in 1:n.subjects){
    mix<-rbinom(1,1,.7)
    constant.params<- c(rnorm(1,100,10),rgamma(1,2,1))
    subject.data.constant<-(mapply(constant.model,t,MoreArgs = list(params=constant.params))  )
    linear.params<- c(rnorm(1,100,10),rgamma(1,2,1),rnorm(1,-15,3))
    subject.data.linear<-(mapply(linear.model,t, MoreArgs = list(params=linear.params))  )

    rtu=(1-mix)*subject.data.linear[1]+mix*subject.data.constant[1]+ rnorm(n.trials,0,10)
    rtp=(1-mix)*subject.data.linear[2]+mix*subject.data.constant[2]+ rnorm(n.trials,0,10)
    
    all.subject.data[[i]]<-list(subject.data=data.frame(rtu=rtu,rtp=rtp), proportion=mix)
  }
  return(all.subject.data)
}



  fake.data<-make.fake.data(c(1:100),50)
rtu<- array(dim=c(length(fake.data),100))
rtp<-array(dim=c(length(fake.data),100))
fake.data
fake.data[[1]]$subject.data

proportion<-c()
for(i in 1: length(fake.data)){
  rtu[i,]<-fake.data[[i]]$subject.data$rtu
  rtp[i,]<-fake.data[[i]]$subject.data$rtp
  proportion<-c(proportion,fake.data[[i]]$proportion)
}
which(proportion==0)

rtu<-rtu
rtp
rtp<-rtp
dim(rtu)
rtp[1,]
plot(1:100,rtp[9,])

model_data<-list(J=50,
                 Tr=100,
                 rtu=rtu,
                 rtp=rtp)


dim(rtu)

# data
t<-greta_array(1:100)
n.trials<-length(t)
n.subjects<-length(rtu[1,])
rtu<-as_data(rtu)
rtp<-as_data(rtp)
rtp[3,]
# variables and priors


mixing.proportions = dirichlet(alpha=rep(1,n.subjects))


sd_mean= gamma(8, 1)

sd_offset = normal(0,1,dim=c(n.subjects))

a_mean = normal(100,10)
b_mean = normal(40,10)
c_mean = normal(-100,10)
d_mean = normal(-3,3)

a_offset = normal(0,1,dim=c(n.subjects))
b_offset = normal(0,1,dim=c(n.subjects))
c_offset = normal(0,1,dim=c(n.subjects))
d_offset = normal(0,1,dim=c(n.subjects))

# operations
a<-t(a_offset*a_mean)
b<-t(b_mean*b_offset)
c<-t(c_mean*c_offset)
d<-t(d_offset*d_mean)

meanu<-1^t%*%a
meanp.1<-1^t%*% (a*b)
meanp.2<- 1^t%*%a * (1^t%*%c + t%*%d)

mu1<-1^t%*%mixing.proportions*meanp.1
mu2<-1^t%*%(1-mixing.proportions)*meanp.2


# likelihood


distribution(rtu)= normal(meanu,sdu)

distribution(rtp)=normal(mu1+mu2,sqrt(sd1^2+sd2^2))
    

# defining the model
m<-model(a_mean,a_offset,b_mean,b_offset,c_mean,c_offset,d_mean,d_offset,sd_mean,sd_offset,mixing.proportions)



# sampling

init=c(a_mean=100,b_mean=40,c_mean= -100,d_mean= -3,
       a_offset= rnorm(n.subjects,0,10),b_offset= rnorm(n.subjects,0,10),c_offset= rnorm(n.subjects,0,10),d_offset= rnorm(n.subjects,0,3),
       sd_mean = rnorm(1,10,1) ,st_offset = rnorm(n.subjects,0,1), mixing.proportions=rbeta(n.subjects,8,1))

init=c(a_mean=1,b_mean=1,c_mean= -1,d_mean= -1,
       a_offset= rnorm(n.subjects,0,1),b_offset= rnorm(n.subjects,0,1),c_offset= rnorm(n.subjects,0,1),d_offset= rnorm(n.subjects,0,1),
       sd_mean = rnorm(1,10,1) ,st_offset = rnorm(n.subjects,0,1), mixing.proportions=rbeta(n.subjects,8,1))


draws <- mcmc(m, n_samples = 10000, warmup=1000,initial_values = init)
    ?mcmc            
c(init)
# plotting
plot(draws)

