library(rstan)

#simple test
y = c(rnorm(100,80,10))
model_data<-list(list(y=y))

model_data
fit <- stan(file = 'Simple_test.stan', data = model_data, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T )


print(fit)
#hierarchical test 

n.subjects=100
n.data =1000

mu_mean0 = 80
mu_meansd= 10

a=16
b=2

mu_mean<-rnorm(1,mu_mean0,mu_meansd)
mu_sd<-rgamma(1,a,b)

mu<-rnorm(n.subjects,mu_mean,mu_sd)
sd<-rgamma(n.subjects, a, b)

y<-array(dim=c(n.subjects,n.data))

for(i in 1:n.subjects){
  y[i,]<-rnorm(n.data,mu[i],sd[i])
}
y

J=n.subjects
Tr=n.data

model_data_hierarchical<-list(list('J' = J, 'Tr' = Tr,
  'y'=y))

fit_hierarchical <- stan(file = 'hierarchical_test.stan', data = model_data_hierarchical, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T )

print(fit_hierarchical,par=c('mu_sd'))
plot(fit_hierarchical,par=c('mu_sd'))

#hierarchical test(reparameterization)



fit_hierarchical_2 <- stan(file = 'hierarchical_test_2.stan', data = model_data_test, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T )

print(fit_hierarchical_2,par=c('mu_sd'))
plot(fit_hierarchical_2,par=c('mu_sd'))

#breakpoint test
Tr=150
y=c(rnorm(100,80,10),rnorm(50,-60,10))

model_data_breakpoint<-list(list(Tr= Tr, y=y ))

fit_breakpoint <- stan(file = 'breakpoint_test.stan', data = model_data_breakpoint, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T )

plot(fit_breakpoint,pars = c('lp'))
print(fit_breakpoint, pars = c('lp[106]'))

#discrete indicator test (binary)

K=2
I=100

mu=c(-100,200)
sigma=c(10,30)
z<-numeric(I)
y=numeric(I)
for(i in 1:I){
  z[i]<- rbinom(1,1,.8)
  y[i]= rnorm(1,mu[(1+z[i])],sigma[(1+z[i])])
}
y
alpha=c(1,1)
mu0=c(0,0)
sigma0=100



model_data_indicator<-list(list(K=K,I=I,y=y,alpha=alpha,mu0=mu0,sigma0=sigma0))

fit_indicator <- stan(file = 'discrete_indicator_test.stan', data = model_data_indicator, iter = 10000,warmup =1000 , 
            chains = 1, verbose = T )
plot(fit_indicator, par=c('log_q_z[40,1]','log_q_z[40,2]'))
plot(fit_indicator, par=c('log_q_z'))
z[40]

print(fit_indicator,par = c(log_q_z[1]))

fit_indicator
pr<-matrix(nrow=100,ncol=2)
j=1
fit<-as.matrix(fit_indicator, par = c(paste0('log_q_z[',j,',1]'),paste0('log_q_z[',j,',2]')))

length(fit[1,])
for(j in 1:100){
  fit<-as.matrix(fit_indicator, par = c(paste0('log_q_z[',j,',1]'),paste0('log_q_z[',j,',2]')))
  probs<-matrix(nrow=9000,ncol=2)
  for(i in 1:9000){
    probs[i,]<-softmax2(fit[i,])
  }

  pr[j,1]<-mean(probs[,1])
  pr[j,2]<-mean(probs[,2])
}


hist(pr[,2])
z.prob<-as.data.frame(pr)
colnames(z.prob)<- c('zero','one')
z.prob
library(ggplot2)
bargraph(z.prob[4,])
rbind(z.prob[,1],z)

a<-as.matrix(z.prob[,1],nrow=100,ncol=1)
b<-as.matrix(z.prob[,2],nrow=100,ncol=1)


data.plot<-data.frame(prob=rbind(a,b))
data.plot
data.plot$val<-rbind(as.matrix(rep(0,100),nrow=100,ncol=1),as.matrix(rep(1,100),nrow=100,ncol=1))
?geom_col
geo
ggplot(data.plot)+
  geom_col(aes(x=val,y=prob),width=.5)+
  theme_bw()



softmax2 <- function(x){ exp(x) / sum(exp(x))}

softmax2(log_q_z)


softmax2(c(2,5))

mcmc <- As.mcmc.list(fit_indicator)


plot(mcmc)
