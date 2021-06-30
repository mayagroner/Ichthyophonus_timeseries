
##Comparison of assigned and actual ages for PWS herring using fishmethods ALkey

rm(list = ls())

library(readr)
require(plyr)
library(tidyr)
library(EnvStats)
require(ggplot2)
library(jtools)
library(dotwhisker)
library(broom)
library(dplyr)
library(ggpmisc)
require(glmmTMB) 
require(lme4)

setwd("C:/Users/mayag/Dropbox/Ichthyophonus/Age-length Data")
#Load surveillance data and merge with severity data
AllHerring <- read_csv("C:/Users/mayag/Dropbox/Ichthyophonus/AK herring samples for R updated 6-4-21 age_QCed.csv", guess_max=5000)

AllHerring<-AllHerring[!is.na(AllHerring$Set),] 
AllHerring<-AllHerring[!is.na(AllHerring$Site),]      
AllHerring$Bin<-floor(AllHerring$Length_mm/10)*10+5
AllInfectedHerring<-AllHerring[AllHerring$IchFinal==1,]
AllInfectedHerring<-AllInfectedHerring[!is.na(AllInfectedHerring$Set),]     

AllHerringSummary<-ddply(AllHerring, .(Year, Site, Bin, Set), summarize, Prevalence=sum(IchFinal, na.rm=TRUE)/length(IchFinal), PrevalenceSample=length(IchFinal))
AllHerringSummary<-AllHerringSummary[AllHerringSummary$Set!="12-8",]
AllHerringSummary<-AllHerringSummary[AllHerringSummary$PrevalenceSample>10,]



#Incorporate age estimates
#download AL keys
Fishmethods_AL_Key_PWS <- read.delim("C:/Users/mayag/Dropbox/Ichthyophonus/Age-length Data/alk_Cordova_fishmethods.txt")

#merge
AllHerring_PWS<-AllHerring[AllHerring$Site=='PWS',]
AllHerringAL<-merge(AllHerring_PWS, Fishmethods_AL_Key_PWS, by=c("Year", "Bin"), all.x=TRUE)# not all samples are there


#Assign an age to each sample
AllHerringAL$Age34<-AllHerringAL$Age3_minus+AllHerringAL$Age4
AllHerringAL$Age345<-AllHerringAL$Age3_minus+AllHerringAL$Age4+AllHerringAL$Age5
AllHerringAL$Age3456<-AllHerringAL$Age3_minus+AllHerringAL$Age4+AllHerringAL$Age5+AllHerringAL$Age6
AllHerringAL$Age34567plus<-AllHerringAL$Age3_minus+AllHerringAL$Age4+AllHerringAL$Age5+AllHerringAL$Age6+AllHerringAL$Age7_plus


AllHerringAL$Random<-runif(length(AllHerringAL$Age3))
AllHerringAL$AssignedAge<-ifelse(is.na(AllHerringAL$Age34567plus),NA, 
                                 ifelse(AllHerringAL$Random<AllHerringAL$Age3, 3,
                                        ifelse(AllHerringAL$Random<AllHerringAL$Age34, 4,
                                               ifelse(AllHerringAL$Random<AllHerringAL$Age345, 5,
                                                      ifelse(AllHerringAL$Random<AllHerringAL$Age3456, 6, 7)))))

check<-ddply(AllHerringAL, .(Year, AssignedAge), summarize, Random=mean(Random))
ggplot(check, aes(x=AssignedAge, y=Random, colour=as.factor(Year)))+geom_point()+geom_line()
mean(AllHerringAL$Random)
############################################
##Compare predicted ages with actual ages in PWS
AgeComp<-AllHerringAL[!is.na(AllHerringAL$Age),]
AgeComp<-AgeComp[AgeComp$Age>2,]
AgeComp<-AgeComp[AgeComp$Age<10,]

AgeComp$Age<-ifelse(AgeComp$Age>7, 7, AgeComp$Age)
AgeComp$AgeDiff<-AgeComp$Age-AgeComp$AssignedAge
mean(AgeComp$AgeDiff, na.rm=TRUE) ##-.75
sd(AgeComp$AgeDiff, na.rm=TRUE)/sqrt(length(AgeComp$AgeDiff))


AgeCompSum<-ddply(AgeComp, .(Year, Age), summarize, meanDiff=mean(AgeDiff, na.rm=TRUE), SampleSize=length(AgeDiff), meanDiffSE<-sd(AgeDiff, na.rm=TRUE)/sqrt(length(AgeDiff)))

ggplot(AgeCompSum, aes(x=Age, y=meanDiff, size=SampleSize))+geom_point()+geom_line()+facet_grid(Year~.)+ geom_hline(yintercept=0, linetype="dashed", color = "red")+ylab("Scale age - Assigned Age")+xlab('Actual Age')+theme_bw()


#Bland Altman method comparison
library(blandr)
library(ggplot2)

comp<-blandr.statistics ( AgeComp$Age, AgeComp$AssignedAge , sig.level=0.95 )

blandr.output.text ( AgeComp$Age, AgeComp$AssignedAge , sig.level=0.95 )
