#############################################################
#libraries for analysis

library(ggplot2)
library(knitr)
library(lme4)
library(multcomp)
library(reshape2)
library(dplyr)
library(stringr)


###########################################################
#processing the data for analysis

#read in data
data = read.csv("catmedata.csv")

#process the data using Hadley-verse
tmp = melt(data, id.vars=c("teamid", "fyd", "studentid"))

catmeData = tmp %>% mutate(time=as.integer(gsub("[a-z]+", "", variable)),type=gsub("\\d", "", variable),studentid=factor(studentid),teamid=factor(teamid))

#removing NA's
catmeData = na.omit(catmeData)

###############################################################
#visualizing data

############
#plot a subset of student's scores over time
sam = sample(unique(catmeData$studentid),10)
catmeSub = catmeData %>% filter(studentid %in% sam)

ggplot(catmeSub,aes(x=time,y=value,colour=studentid)) + geom_line() + facet_wrap(~type,scales="free_y")
#seems there are several possible time trends - flat, down or up - possible fit as random intercept + random slope

###############
#visualize correlation between score types for a single time point
catmeSub = catmeData %>% filter(time==1)

catmeWide = split(catmeSub, catmeSub$type) %>% lapply(function(x) {colnames(x)[colnames(x) == "value"] <- unique(x$type); x %>% select(-variable, -type)}) %>% Reduce(function(x, y) left_join(x, y, by=c("studentid", "fyd", "time", "teamid")), .)

pairs(catmeWide[,c(4,6:14)])
#score types are VERY correlated - need to adjust for multiplicity in a manner appropraite for high (positive) dependence

##############
#visualize differences between fyd
sam = sample(unique(catmeData$studentid),25)
catmeSub = catmeData %>% filter(studentid %in% sam)

ggplot(catmeSub,aes(x=time,y=value,colour=fyd,group=studentid)) + geom_line() + facet_wrap(~type,scales="free_y")
#hard to tell - the signal will be weak

##############
#plot histograms of student's scores
ggplot(catmeData %>% filter(time==1), aes(x=value)) + geom_histogram() + facet_wrap(~type,scales="free_x")
#scores normal enough for linear mixed effects models 


###########################################################
#Mixed Effects Modeling

############
#first, full model for one type of score at a time & visualize results
metric = "adj"

#look over all time
fit0 = lmer(value~(1|studentid),data=catmeData%>%filter(type==metric))
fit1 = lmer(value~(1+time|studentid),data=catmeData%>%filter(type==metric))
fit2 = lmer(value~(1+time|studentid)+(1|teamid),data=catmeData%>%filter(type==metric))
anova(fit0,fit1,fit2)
#leave out teamID 

fit2 = lmer(value~(1+time|studentid)+time,data=catmeData%>%filter(type==metric))
anova(fit0,fit1,fit2)
#time on average not an effect

fit3 = lmer(value~(1+time|studentid)+fyd,data=catmeData%>%filter(type==metric))
anova(fit0,fit1,fit3)
#fyd is significant

#simulate new data to understand and interpret model
nums = 20
fakedata = data.frame(time = rep(1:4, nums*2), fyd = rep(c(TRUE, FALSE), each=4, times=nums), studentid=rep(1:(nums*2), each=4))
fakeval = simulate(fit3,newdata=fakedata,allow.new.levels=TRUE)
fakedata = cbind(fakedata,fakeval)

ggplot(fakedata,aes(x=time,y=sim_1,colour=fyd))+geom_line(aes(group=studentid)) + stat_smooth(method="lm",se=FALSE,size=2)


#####################
#Model across all 4 time points - fits linear trends over time
#AIC used to determine appropriate terms for model
#Used Holm's FWER and Benjamini-Hochberg FDR method for multiple comparisons

pvals = NULL
for(metric in unique(catmeData$type)){

#random intercept for students
  fit0 = lmer(value~(1|studentid),data=catmeData%>%filter(type==metric))

#random intercept and slope over time for students
  fit1 = lmer(value~(1+time|studentid),data=catmeData%>%filter(type==metric))
  res1 = anova(fit0,fit1)
  if(res1$AIC[2] < res1$AIC[1]){
  
  #random intercept for teams
    fit2 = lmer(value~(1+time|studentid)+(1|teamid),data=catmeData%>%filter(type==metric))
    res2 = anova(fit1,fit2)
    if(res2$AIC[2] < res2$AIC[1]){
    
    #linear time trend - fixed effect
      fit3 = lmer(value~(1+time|studentid)+(1|teamid)+time,data=catmeData%>%filter(type==metric))
      res3 = anova(fit2,fit3)
      if(res3$AIC[2] < res3$AIC[1]){

      #fixed effect for ENGI 120 status
        fit4 = lmer(value~(1+time|studentid)+(1|teamid)+time+fyd,data=catmeData%>%filter(type==metric))
      res4 = anova(fit0,fit1,fit2,fit3,fit4)
      }else{
        fit4 = lmer(value~(1+time|studentid)+(1|teamid)+fyd,data=catmeData%>%filter(type==metric))
        res4 = anova(fit0,fit1,fit2,fit4)
      }
    }else{
      fit3 = lmer(value~(1+time|studentid)+time,data=catmeData%>%filter(type==metric))
      res3 = anova(fit1,fit3)
      if(res3$AIC[2] < res3$AIC[1]){
        fit4 = lmer(value~(1+time|studentid)+time+fyd,data=catmeData%>%filter(type==metric))
        res4 = anova(fit0,fit1,fit3,fit4)
      }else{
        fit4 = lmer(value~(1+time|studentid)+fyd,data=catmeData%>%filter(type==metric))
        res4 = anova(fit0,fit1,fit4)
      }
    }
  }else{
    fit2 = lmer(value~(1|studentid)+(1|teamid),data=catmeData%>%filter(type==metric))
    res2 = anova(fit0,fit2)
    if(res2$AIC[2] < res2$AIC[1]){
      fit3 = lmer(value~(1|studentid)+(1|teamid)+time,data=catmeData%>%filter(type==metric))
      res3 = anova(fit2,fit3)
      if(res3$AIC[2] < res3$AIC[1]){
        fit4 = lmer(value~(1|studentid)+(1|teamid)+time+fyd,data=catmeData%>%filter(type==metric))
        res4 = anova(fit0,fit2,fit3,fit4)
      }else{
        fit4 = lmer(value~(1|studentid)+(1|teamid)+fyd,data=catmeData%>%filter(type==metric))
        res4 = anova(fit0,fit2,fit4)
      }
    }else{
      fit3 = lmer(value~(1|studentid)+time,data=catmeData%>%filter(type==metric))
      res3 = anova(fit0,fit3)
      if(res3$AIC[2] < res3$AIC[1]){
        fit4 = lmer(value~(1|studentid)+time+fyd,data=catmeData%>%filter(type==metric))
        res4 = anova(fit0,fit3,fit4)
      }else{
        fit4 = lmer(value~(1|studentid)+fyd,data=catmeData%>%filter(type==metric))
        res4 = anova(fit0,fit4)
      }
    }
  }
  print(paste("Score Type:",metric,sep=" "))
  print(res4)
  pvals = c(pvals,tail(res4$Pr,n=1))
}

#p-value results
mod1 = data.frame(Type=(unique(catmeData$type)),P.Val=pvals,Adj.P.Val.Holms=p.adjust(pvals,method="holm"),Adj.P.Val.FDR=p.adjust(pvals,method="BH"))


#####################
#mixed effects models - looking at individual time points
#is there an interaction between time and ENGI 120 status?


catmeDataT = catmeData %>% mutate(time=factor(time))
pvals = NULL
for(metric in unique(catmeData$type)){                                        
  fit0 = lmer(value~(1|studentid),data=catmeDataT%>%filter(type==metric))
  fit1 = lmer(value~(1|studentid)+(1|teamid),data=catmeDataT%>%filter(type==metric))
  res1 = anova(fit0,fit1)
  if(res1$AIC[2] < res1$AIC[1]){
    fit2 = lmer(value~(1|studentid)+(1|teamid)+time,data=catmeDataT%>%filter(type==metric))
    fit3 = lmer(value~(1|studentid)+(1|teamid)+time+fyd,data=catmeDataT%>%filter(type==metric))
    fit4 = lmer(value~(1|studentid)+(1|teamid)+time*fyd,data=catmeDataT%>%filter(type==metric))
    res4 = anova(fit0,fit2,fit3,fit4)
    print(paste("Score Type:",metric,sep=" "))
    print(res4)
    if(tail(res4$Pr,n=1)<.05){
      pw.pvals = NULL
      for(i in 1:4){
        f01 = lmer(value~(1|teamid),data=(catmeDataT%>%filter(type==metric))%>%filter(time==i))
        f1 = lmer(value~(1|teamid)+fyd,data=(catmeDataT%>%filter(type==metric))%>%filter(time==i))
        res1 = anova(f01,f1)
        pw.pvals = rbind(pw.pvals,c(paste("Time ",i,sep=""),res1$Pr[2]))
      }
      pwres = data.frame(Test=pw.pvals[,1],P.Vals=as.numeric(pw.pvals[,2]),Adj.P.Vals=p.adjust(as.numeric(pw.pvals[,2]),method="holm"))
      print(paste("Score Type:",metric,"Pair-wise Tests",sep=" "))
      print(pwres)
    }
  }else{
    fit2 = lmer(value~(1|studentid)+time,data=catmeDataT%>%filter(type==metric))
    fit3 = lmer(value~(1|studentid)+time+fyd,data=catmeDataT%>%filter(type==metric))
    fit4 = lmer(value~(1|studentid)+time*fyd,data=catmeDataT%>%filter(type==metric))
    res4 = anova(fit0,fit2,fit3,fit4)
    print(paste("Score Type:",metric,sep=" "))
    print(res4)
    if(tail(res4$Pr,n=1)<.05){
      pw.pvals = NULL
      for(i in 1:4){
        f01 = lmer(value~(1|teamid),data=(catmeDataT%>%filter(type==metric))%>%filter(time==i))
        f1 = lmer(value~(1|teamid)+fyd,data=(catmeDataT%>%filter(type==metric))%>%filter(time==i))
        res1 = anova(f01,f1)
        pw.pvals = rbind(pw.pvals,c(paste("Time ",i,sep=""),res1$Pr[2]))
      }
      pwres = data.frame(Test=pw.pvals[,1],P.Vals=as.numeric(pw.pvals[,2]),Adj.P.Vals=p.adjust(as.numeric(pw.pvals[,2]),method="holm"))
      print(paste("Score Type:",metric,"Pair-wise Tests",sep=" "))
      print(pwres)
    }
  }
  pvals = c(pvals,tail(res4$Pr,n=1))
}
mod2 = data.frame(Type=(unique(catmeData$type)),P.Val=pvals,Adj.P.Val.Holms=p.adjust(pvals,method="holm"),Adj.P.Val.FDR=p.adjust(pvals,method="BH"))

##############################################################
##############################################################


