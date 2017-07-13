install.packages('tidyr')

library(dplyr)
library(tidyr)


##extract raw data##
data.all <- read.csv('raw_data/all-data.csv', header=T)
data.test <- subset(data.all, practice==0)

n_subjects <- length(unique(data.test$subject))

data.test$appearanceCount <- rep(as.vector(sapply(1:72, function(x){return(rep(x,12))})), n_subjects)

data.for.model <- subset(data.test, correct==1 & triple_type == 'critical' & triple_position%in%c(1,3), c('subject','cond','appearanceCount','triple_position','rt'))

subject.id <- data.frame(original=unique(data.test$subject))
subject.id$original.numeric <- as.numeric(factor(subject.id$original))

data.for.model$subject <- as.numeric(factor(data.for.model$subject))
data.for.model$cond <- factor(as.character(data.for.model$cond), levels=c('known','unknown','one'))
data.for.model$cond <- as.numeric(data.for.model$cond)


# data.for.model<- data.for.model %>% spread('triple_position','rt')
# colnames(data.for.model) <- c('subject', 'cond', 't', 'unpredictable','predictable')


colnames(data.for.model) <- c('subject', 'cond', 't', 'predictable','rt')

## using differnece between N and W response times as measure for fitting ##

#data.for.model$difference <- -(data.for.model$rt.N - data.for.model$rt.W) #as.vector(mapply(function(x,y){if(!is.na(x)&& !is.na(y)){return(x-y)} else{return(NA)}},data.for.model$rt.N, data.for.model$rt.W))

##organize data by subject, discard NAs##
#data.for.model <- subset(data.for.model,!is.na(data.for.model$difference))

write.csv(data.for.model, file="generated_data/extracted-data-niw.csv", row.names=F)
 
# data.by.subject<-list()
# for(i in 1:n_subjects){
#   data.by.subject[[i]]<-data[which(data$subject==i),]
# }



