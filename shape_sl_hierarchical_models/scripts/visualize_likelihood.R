library(jsonlite)
library(ggplot2)
library(gridExtra)

source('scripts/DBDA2E-utilities.R')
source('scripts/priors.R')


from.to<-function(list,n.params){
  n.priors<-length(list)
  
  from<-array(dim=c(n.priors,n.params))
  to<-array(dim=c(n.priors,n.params))
  out<-data.frame(from=numeric(length=n.params),to=numeric(length=n.params))
  
  if(n.priors == 1){
    
    for(j in 1:n.params){
      for(i in 1:n.priors){
        
        
        prior<-list[[i]]
        
        prior<-extract.prior(prior)
        pr<-prior[[j]]
        
        fn<-pr[1]
        par1<-as.numeric(pr[2])
        par2<-as.numeric(pr[3])
        
        
        
        if(fn == 'normal'){
          from[i,j]= par1-2.5*par2
          to[i,j]= par1 + 2.5*par2
          
        }
        if(fn == 'gamma'){
          from[i,j] = 0
          to[i,j] = par1 + 3*par2
        }
        if(fn== 'beta'){
          from[i,j] = 0
          to[i,j] = 1
        }
        if(fn == 'unif'){
          from[i,j] = par1
          to[i,j] = par2
        }
        if(fn == 'exp'){
          from[i,j] = 0
          to[i,j] = par1*3
        }
        
      }
      out$from[j]<- min(from[,j])
      out$to[j]<-max(to[,j])
    }
  }
  else{
    
    for(j in 1:n.params){
      
      prior<-list
      
      prior<-extract.prior(prior)
      pr<-prior[[j]]
      
      fn<-pr[1]
      par1<-as.numeric(pr[2])
      par2<-as.numeric(pr[3])
      
      
      
      if(fn == 'normal'){
        out$from[j]= par1-2.5*par2
        out$to[j]= par1 + 2.5*par2
        
      }
      if(fn == 'gamma'){
        out$from[j] = 0
        out$to[j] = par1 + 3*par2
      }
      if(fn== 'beta'){
        out$from[j] = 0
        out$to[j] = 1
      }
      if(fn == 'unif'){
        out$from[j] = par1
        out$to[j] = par2
      }
      if(fn == 'exp'){
        out$from[j] = 0
        out$to[j] = par1*3
      }
      
    }
  }
  return(out)
}


plot.likelihood<-function(lik,from,to){
  
  fn<-lik[1]
  par1<-as.numeric(lik[2])
  par2<-as.numeric(lik[3])
  title<-lik[4]
  
  by<-(to-from)/1000
  x<-seq(from,to,by)
  
  if(fn == 'normal'){
    
    y<- sapply(x, function(x){dnorm(x,par1,par2)}) 
    
    title.2<- paste('Normal(',toString(round(par1,2)),',',toString(round(par2,2)), ')',sep="")
  }
  if(fn == 'gamma'){
    
    a<-gammaShRaFromModeSD(par1,par2)$shape
    b<-gammaShRaFromModeSD(par1,par2)$rate
    y<- sapply(x, function(x){dgamma(x,a,b)})
    
    title.2<- paste('Gamma(',toString(round(a,2)),',',toString(round(b,2)), ')',sep="")
  }
  if(fn== 'beta'){
    
    a<-betaABfromModeKappa(par1,par2)$a
    b<-betaABfromModeKappa(par1,par2)$b
    
    y<- sapply(x, function(x){dbeta(x,a,b)})  
    
    title.2<- paste('Beta(',toString(round(a,2)),',',toString(round(b,2)), ')',sep="")
  }
  if(fn == 'unif'){
    
    y<- sapply(x, function(x){dunif(x,par1,par2)}) 
    
    title.2<- paste('Uniform(',toString(round(par1,2)),',',toString(round(par2,2)), ')',sep="")
  }
  if(fn == 'exp'){
    
    y<- sapply(x, function(x){dexp(x,par1)}) 
    
    title.2<- paste('Exponential(',toString(round(par1,2)), ')',sep="")
  }
  prior.plot<- data.frame(x=x,y=y,title= rep(title,length(x)),distribution=rep(title.2,length(x)))
  return(prior.plot)
}



graph.prior.multi<-function(plot.data){
  title<-plot.data$title[1]
  
  p<- 
    ggplot(plot.data)+
    geom_line(aes(x=x, y =y,col=distribution))+
    ggtitle(title)+
    labs(y= 'Probability Density', x = 'Parameter')+
    scale_size(range=c(5,20))+
    ylim(c(0, (max(plot.data$y)*2)))+
    theme_minimal()+
    theme(legend.justification=c(.7,.7), legend.position=c(.9,.8), 
          legend.background = element_rect(fill = 'white', colour = 'black', size = .5, linetype = NULL, color = NULL),
          legend.key.size = unit(0.025,'npc'),
          legend.text = element_text(size=7),
          legend.title = element_text(size = 8))
  
  p<-ggplotGrob(p)  
  return(p)
}


visualize.priors.multi<-function(list){
  n.priors<-length(list)
  prior<-list()
  
  if(n.priors>1){
    
    for(i in 1:n.priors){
      prior[[i]]<-extract.prior(list[[i]])
    }
    
    n.params<-length(prior[[1]])
    
    from<-from.to(list,n.params)$from
    to<-from.to(list,n.params)$to
    
    p<-list()
    plot<-list()
    
    
    for(i in 1:n.params){
      plot.data<-numeric()
      for(j in 1:n.priors){
        
        pr<-prior[[j]][[i]]
        plot[[i]]<- plot.prior(pr,from[i],to[i])
        
        plot.data<-rbind(plot.data,plot[[i]])
      }
      p[[i]]<-graph.prior.multi(plot.data)
    }
  }
  
  else{
    prior<-extract.prior(list[[1]])
    n.params<-length(prior)
    
    from<-from.to(list,n.params)$from
    to<-from.to(list,n.params)$to
    
    p<-list()
    plot<-list()
    
    for(i in 1:n.params){
      plot[[i]]<-plot.prior(prior[[i]], from[i],to[i])
      
      p[[i]]<-graph.prior.multi(plot[[i]])
    }
  }
  return(p)
}


grid.arrange(grobs = visualize.priors.multi(list('prior.1.json')) )
