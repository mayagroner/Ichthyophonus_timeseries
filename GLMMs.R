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


#load severity data and summarize infection data by proportion of ventricle that is infected. 
Areas <- read_csv("Historical measurements 2009-2019.csv")

AreasSummary<-ddply(Areas, .(Tissue, Sample, Set, SetID), summarize, Area=sum(Area_um2))

InfectionSummary<-spread(AreasSummary, Tissue, Area)
colnames(InfectionSummary)[colnames(InfectionSummary)=="Schizont"] <- "SchizArea"
colnames(InfectionSummary)[colnames(InfectionSummary)=="Ventricle"] <- "VentricleArea"


InfectionSummary$PercentInfected<-ifelse(is.na(InfectionSummary$VentricleArea)==TRUE, NA, 
                                         ifelse(is.na(InfectionSummary$SchizArea)==TRUE,0, InfectionSummary$SchizArea/InfectionSummary$VentricleArea))
InfectionSummary<-subset(InfectionSummary,select= -c(Set, SetID))

#Load surveillance data and merge with severity data
Herring <- read_csv("AK herring disease samples for R updated 11-15-20.csv", guess_max=5000)

AllHerring<-merge(InfectionSummary, Herring, by=c("Sample"), all=TRUE) 
AllHerring<-AllHerring[!is.na(AllHerring$Site),]      
AllHerring$LengthBin<-floor(AllHerring$Length_mm/20)*20
AllHerring$Bin<-floor(AllHerring$Length_mm/10)
AllInfectedHerring<-AllHerring[AllHerring$IchFinal==1,]
AllInfectedHerring<-AllInfectedHerring[!is.na(AllInfectedHerring$Site),]     

AllHerringSummary<-ddply(AllHerring, .(Year, Site, LengthBin, Set), summarize, Prevalence=sum(IchFinal, na.rm=TRUE)/length(IchFinal), PrevalenceSample=length(IchFinal), PercSeverity=geoMean(PercentInfected+0.000000001, na.rm=TRUE), SeveritySample=sum(IchFinal))
AllHerringSummary<-AllHerringSummary[AllHerringSummary$Set!="12-8",]
AllHerringSummary<-AllHerringSummary[AllHerringSummary$PrevalenceSample>10,]
AllHerringSummary$LogPercSeverity<-log10(AllHerringSummary$PercSeverity)


#Incorporate age estimates
#download AL keys
Weighted_AL_Key_PWS <- read_delim("alk_Cordova_fishmethods.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Weighted_AL_Key_PWS$Site<-"PWS"
Weighted_AL_Key_PWS<-subset(Weighted_AL_Key_PWS, select=-c(nl, Age1, Age2, Age3))
Weighted_AL_Key_Sitka <- read_delim("alk_Sitka_fishmethods.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
Weighted_AL_Key_Sitka$Site<-"Sitka"
Weighted_AL_Key_Sitka<-subset(Weighted_AL_Key_Sitka, select=-c(nl, Age2, Age3))
Weighted_ALK<-rbind(Weighted_AL_Key_PWS, Weighted_AL_Key_Sitka)

library(tidyr)
Weighted_AL_long<-gather(Weighted_ALK, Age, ProbAge, Age4:Age3_minus)
Weighted_AL_long$Bin<-Weighted_AL_long$Bin/10
Weighted_AL_long$Age<-ifelse(Weighted_AL_long$Age=='Age3_minus', 3, 
                             ifelse(Weighted_AL_long$Age=='Age4', 4, 
                                    ifelse(Weighted_AL_long$Age=='Age5', 5, 
                                           ifelse(Weighted_AL_long$Age=='Age6', 6, 7))))
Weighted_AL_long$Cohort<-Weighted_AL_long$Year-as.numeric(Weighted_AL_long$Age)
Weighted_AL_long_Sitka<-Weighted_AL_long[Weighted_AL_long$Site=='Sitka',]
Weighted_AL_long_PWS<-Weighted_AL_long[Weighted_AL_long$Site=='PWS',]

##Graph Age length key
p<-ggplot(Weighted_AL_long_Sitka, aes(x=Bin, y=ProbAge, fill=as.factor(Cohort)))+geom_col()+facet_grid(Year~Age)+xlim(14, 28)+ scale_y_continuous(breaks=seq(0,1,.5))
p<-p+xlab('Fork Length (cm)')+ylab('Probability')+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ scale_fill_discrete(name = "Year class", l=40, c=35)+ theme(text = element_text(size = 8))+ labs(title = "A: Sitka Sound age-length key")

p1<-ggplot(Weighted_AL_long_PWS, aes(x=Bin, y=ProbAge, fill=as.factor(Cohort)))+geom_col()+facet_grid(Year~Age)+xlim(14, 28)+ scale_y_continuous(breaks=seq(0,1,.5))
p1<-p1+xlab('Fork Length (cm)')+ylab('Probability')+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ scale_fill_discrete(name = "Year class", l=40, c=35)+ theme(text = element_text(size = 8))+labs(title = "B: Prince William Sound age-length key")

require(gridExtra)
grid.arrange(p, p1)


###merge AL Key with disease data
Weighted_AL<-rbind(Weighted_AL_Key_PWS, Weighted_AL_Key_Sitka)
Weighted_AL$Bin<-Weighted_AL$Bin/10-.5


#merge key with additional data
AllHerringAL<-merge(AllHerring, Weighted_AL, by=c("Site", "Year", "Bin"), all.x=TRUE)


#Change all really large fish w/out an assigned age to 'Age 7 plus'
AllHerringAL$Age7_plus<-ifelse(!is.na(AllHerringAL$Age7_plus), AllHerringAL$Age7_plus,ifelse(AllHerringAL$Bin>25, 1, NA))

#Remove fish that are below the 120 mm fork length (210 total fish)
AllHerringAL<-AllHerringAL[(AllHerringAL$Bin>12),]
AllHerringAL<-AllHerringAL[!is.na(AllHerringAL$Bin),]
AllHerringAL<-AllHerringAL[!is.na(AllHerringAL$IchFinal),]

#Change NA to 0 for age3-age6
AllHerringAL$Age3<-ifelse(is.na(AllHerringAL$Age3_minus), 0, AllHerringAL$Age3_minus)
AllHerringAL$Age4<-ifelse(is.na(AllHerringAL$Age4), 0, AllHerringAL$Age4)
AllHerringAL$Age5<-ifelse(is.na(AllHerringAL$Age5), 0, AllHerringAL$Age5)
AllHerringAL$Age6<-ifelse(is.na(AllHerringAL$Age6), 0, AllHerringAL$Age6)

#summarize to create dataset on prevalence
AllHerringPrev<-ddply(AllHerringAL,  .(Year, Site),summarize, Age3Prev=sum(Age3*IchFinal, na.rm=TRUE)/sum(Age3, na.rm=TRUE), Age4Prev=sum(Age4*IchFinal, na.rm=TRUE)/sum(Age4, na.rm=TRUE), Age5Prev=sum(Age5*IchFinal, na.rm=TRUE)/sum(Age5, na.rm=TRUE), Age6Prev=sum(Age6*IchFinal, na.rm=TRUE)/sum(Age6, na.rm=TRUE), Age7Prev=sum(Age7_plus*IchFinal, na.rm=TRUE)/sum(Age7_plus, na.rm=TRUE))


#Wide to Long Form
AllHerringPrevLong<-gather(AllHerringPrev, Age, Prevalence, Age3Prev:Age7Prev, factor_key=TRUE)
AllHerringPrevLong$Age<-ifelse(AllHerringPrevLong$Age=='Age3Prev', 3,
                               ifelse(AllHerringPrevLong$Age=='Age4Prev', 4, 
                                      ifelse(AllHerringPrevLong$Age=='Age5Prev', 5,
                                             ifelse(AllHerringPrevLong$Age=='Age6Prev', 6, 7))))

#Samplesize by age
AllHerringSampSize<-ddply(AllHerringAL,  .(Year, Site),summarize, Age3SS=sum(Age3, na.rm=TRUE), Age4SS=sum(Age4, na.rm=TRUE), Age5SS=sum(Age5, na.rm=TRUE), Age6SS=sum(Age6, na.rm=TRUE), Age7SS=sum(Age7_plus, na.rm=TRUE))#length(x[!is.na(x)]), Age4Prev=sum(Age4*IchFinal, na.rm=TRUE), Age5Prev=sum(IchFinal, na.rm=TRUE)/sum(Age5, na.rm=TRUE), Age6Prev=sum(Age6*IchFinal, na.rm=TRUE)/sum(Age6, na.rm=TRUE), Age7Prev=sum(Age7_plus*IchFinal, na.rm=TRUE)/sum(Age7_plus, na.rm=TRUE))

#Wide to Long Form
AllHerringSSLong<-gather(AllHerringSampSize, Age, SampleSize, Age3SS:Age7SS, factor_key=TRUE)
AllHerringSSLong$Age<-ifelse(AllHerringSSLong$Age=='Age3SS', 3,
                               ifelse(AllHerringSSLong$Age=='Age4SS', 4, 
                                      ifelse(AllHerringSSLong$Age=='Age5SS', 5,
                                             ifelse(AllHerringSSLong$Age=='Age6SS', 6, 7))))



AllHerringPrevLong$Cohort<-as.factor(AllHerringPrevLong$Year-AllHerringPrevLong$Age)
AllHerringPrevLong$Age<-as.factor(AllHerringPrevLong$Age)
AllHerringSSLong$Age<-as.factor(AllHerringSSLong$Age)

AllHerringPrevLong<-merge(AllHerringPrevLong, AllHerringSSLong, by=c("Site", "Year", "Age"))

#exclude data with sample size < 5
AllHerringPrevLong$Prevalence<-ifelse(AllHerringPrevLong$SampleSize<5, NA, AllHerringPrevLong$Prevalence)
AllHerringPrevLong$SampleSize<-signif(AllHerringPrevLong$SampleSize, digits = 3)

AllHerringPrevLong$Location<-ifelse(AllHerringPrevLong$Site=='Sitka', 'Sitka Sound', 
                         ifelse(AllHerringPrevLong$Site=='PWS', "Prince William Sound", AllHerringPrevLong$Site))

#Graph prevalence over time by age

PrevAge<-ggplot(AllHerringPrevLong, aes(x=as.numeric(Year), y=Prevalence, color=Cohort))+geom_point()+geom_smooth(method='lm', se=FALSE)
PrevAge+facet_grid(Location~.)+xlab('Sampling year')+ylab('Infection prevalence')+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                            panel.background = element_blank())+ labs(color = "Year class")


#Assign an age to each sample
AllHerringAL$Age34<-AllHerringAL$Age3+AllHerringAL$Age4
AllHerringAL$Age345<-AllHerringAL$Age3+AllHerringAL$Age4+AllHerringAL$Age5
AllHerringAL$Age3456<-AllHerringAL$Age3+AllHerringAL$Age4+AllHerringAL$Age5+AllHerringAL$Age6
AllHerringAL$Age34567plus<-AllHerringAL$Age3+AllHerringAL$Age4+AllHerringAL$Age5+AllHerringAL$Age6+AllHerringAL$Age7_plus


AllHerringAL$Random<-runif(length(AllHerringAL$Age3))
AllHerringAL$AssignedAge<-ifelse(is.na(AllHerringAL$Age34567plus),NA, 
                              ifelse(AllHerringAL$Random<AllHerringAL$Age3, 3,
                                  ifelse(AllHerringAL$Random<AllHerringAL$Age34, 4,
                                      ifelse(AllHerringAL$Random<AllHerringAL$Age345, 5,
                                          ifelse(AllHerringAL$Random<AllHerringAL$Age3456, 6, 7)))))


AllHerringALSev<-ddply(AllHerringAL, .(Year, Site, AssignedAge), summarize, Severity=mean(PercentInfected, na.rm=TRUE), SeveritySD=sd(PercentInfected, na.rm=TRUE), SampleSize=length(PercentInfected[!is.na(PercentInfected)]))
AllHerringALSev$Cohort<-AllHerringALSev$Year-AllHerringALSev$AssignedAge 
AllHerringALSev$Age<-as.factor(AllHerringALSev$AssignedAge)
AllHerringALSev$Cohort<-as.factor(AllHerringALSev$Cohort)

#ddply(AllHerringALSev, .(Site), summarize, ssTotal=sum(SampleSize, na.rm=TRUE))
AllHerringALSev$Location<-ifelse(AllHerringALSev$Site=='PWS', 'Prince William Sound',
                                 ifelse(AllHerringALSev$Site=='Sitka', 'Sitka Sound', AllHerringALSev$Site))

#Graph severity over time and year class
pSev<-ggplot(AllHerringALSev, aes(x=Year, y=Severity, fill=Cohort))+geom_bar(stat='identity')+facet_grid(Cohort~Site)
pSev+geom_text(aes(label=SampleSize), size=3, vjust=-.25)

pSev4<-ggplot(AllHerringALSev, aes(x=(as.numeric(Year)), y=Severity, color=Cohort))+geom_point(stat='identity')+geom_smooth(method='lm', se=FALSE)+facet_grid(Location~.)
pSev4+theme_bw()+xlab('Sampling year')+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                          panel.background = element_blank())+ labs(color = "Year class")+ylab('Infection intensity') +scale_x_continuous(breaks=seq(2009, 2019, 2))+ xlim(2009, 2019)##+geom_text(aes(label=SampleSize), size=3, vjust=-.25, color='black')

#No strong correlations between severity and prevalence
names(AllHerringALSev)[names(AllHerringALSev=="AssignedAge")]<-'Age'
AllHerringALSevPrev<-merge(AllHerringPrevLong, AllHerringALSev, by=c('Year', 'Site', 'Age', 'Cohort'), all=TRUE)
p<-ggplot(AllHerringALSevPrev, aes(x=Prevalence, y=Severity, color=as.factor(Year)))+geom_point()+geom_smooth(method="lm", se=FALSE)
p+facet_grid(Site~.)+theme_bw()+ theme(legend.title=element_blank()) 


######################################################################################
###Run GLMMS to identify correlates of infection severity and infection prevalence

#Upload environmental covariates

Env1<- read_csv("Summarized PDO NPGO.csv")
Env2 <- read_csv("SummarizedUpwelling.csv")
Env<-merge(Env1, Env2, by='Year')
TempAnomPWS<- read_csv("PWSSC 2m Temp Anomaly.csv")

#Upload Demographic covariates
RecruitmentPWS<-read_csv("PWS_Recruitment.csv")
RecruitmentSitka<-read_csv("Sitka Recruitment.csv")
Recruitment<-rbind(RecruitmentPWS, RecruitmentSitka)

#Quick look at temp anomalies. How much missing data is there?

TempAnomPWS<-TempAnomPWS[TempAnomPWS$Year>2007,]
TempAnomPWS$SSTAnom_C<-as.numeric(TempAnomPWS$SSTAnom_C)
TempAnomPWS$Year<-as.factor(TempAnomPWS$Year)

#plot it
p<-ggplot(TempAnomPWS, aes(x=Month, y=SSTAnom_C, color=Year))+geom_point()+geom_line()
p

#July is the only month with temp anomalies for all years
TempAnomPWSJuly<-as.data.frame(TempAnomPWS[TempAnomPWS$Month==7,])



#merge datasets
AllHerringAL<-merge(AllHerringAL, Env, by=c('Year','Site'))
AllHerringAL<-merge(AllHerringAL, Recruitment, by=c('Year', 'Site'))
AllHerringAL$UniqueSet<-as.factor(AllHerringAL$UniqueSet)

AllHerringAL<-AllHerringAL[!is.na(AllHerringAL$Year),]
AllHerringAL$PercentInfected<-ifelse(AllHerringAL$IchFinal==0,0,AllHerringAL$PercentInfected)

AllHerringAL$Year<-as.numeric(AllHerringAL$Year)
AllHerringAL$Cohort<-AllHerringAL$Year-AllHerringAL$AssignedAge
AllHerringAL$Year<-as.factor(AllHerringAL$Year)
AllHerringAL$Cohort<-as.factor(AllHerringAL$Cohort)

AllHerringALSitka<-AllHerringAL[AllHerringAL$Site=='Sitka',]
AllHerringALSitka<-merge(AllHerringALSitka, RecruitmentSitka, by='Year')
AllHerringALPWS<-AllHerringAL[AllHerringAL$Site=='PWS',]
AllHerringALPWS<-merge(AllHerringALPWS, RecruitmentPWS, by='Year')
AllHerringALInfected<-AllHerringAL[AllHerringAL$IchFinal==1,]



#differentiate between symptomatic and asymptomatic infections (Asymptomatic are culture positive, histo negative, symptomatic are culture positive and histo positive)
AllHerringALInfected$Symptomatic<-ifelse(is.na(AllHerringALInfected$VentricleArea),NA,
                                         ifelse(is.na(AllHerringALInfected$PercentInfected), NA, 
                                                ifelse(AllHerringALInfected$PercentInfected==0,0, 1)))
                                                
AllHerringALInfected<-AllHerringALInfected[!is.na(AllHerringALInfected$Year),]
AllHerringALInfected<-AllHerringALInfected[!is.na(AllHerringALInfected$Symptomatic),]

AllHerringALInfectedSitka<-AllHerringALInfected[AllHerringALInfected$Site=='Sitka',]
AllHerringALInfectedSitka<-merge(AllHerringALInfectedSitka, RecruitmentSitka, by='Year')
AllHerringALInfectedPWS<-AllHerringALInfected[AllHerringALInfected$Site=='PWS',]
AllHerringALInfectedPWS<-merge(AllHerringALInfectedPWS, RecruitmentPWS, by='Year')
AllHerringALInfectedPWS<-merge(AllHerringALInfectedPWS, TempAnomPWSJuly, by="Year")


###Create dataset for Prevalence reconstruction
AllHerringAL$Asymptomatic<-ifelse(AllHerringAL$IchFinal==0, 0, ifelse(is.na(AllHerringAL$VentricleArea),0, ifelse(is.na(AllHerringAL$SchizArea), 1, 0)))
AllHerringAL$UnknownSympt<-ifelse(AllHerringAL$IchFinal==0, 0, ifelse(is.na(AllHerringAL$VentricleArea),1,0))
AllHerringAL$Symptomatic<-ifelse(AllHerringAL$IchFinal==0, 0, ifelse(AllHerringAL$PercentInfected>0,1,0))
AllHerringAL$Symptomatic<-ifelse(AllHerringAL$UnknownSympt==1, NA, AllHerringAL$Symptomatic)
AllHerringAL$Asymptomatic<-ifelse(AllHerringAL$UnknownSympt==1, NA, AllHerringAL$Asymptomatic)


AllHerringSummary2<-ddply(AllHerringAL, .(Site, Year, AssignedAge, Cohort), summarize, SampleSize=length(Symptomatic), Prevalence=sum(IchFinal, na.rm=TRUE)/length(IchFinal), SPrevalence=(sum(Symptomatic, na.rm=TRUE)/length(Symptomatic)), APrevalence=sum(Asymptomatic, na.rm=TRUE)/length(Asymptomatic), ASUnknown=sum(UnknownSympt, na.rm=TRUE)/length(UnknownSympt))
names(AllHerringSummary2)[names(AllHerringSummary2)=="AssignedAge"]<-"Age"
#write.csv(AllHerringSummary2,"C:\\Users\\mayag\\Dropbox\\Ichthyophonus\\Prevalence_age_site_v3.csv", row.names = FALSE)



####################################################################
####################################################################
#Model selection for Zero-inflated models of infection patterns. 
AllHerringALInfected$PercentInfected<-ifelse(is.na(AllHerringALInfected$PercentInfected), 0, AllHerringALInfected$PercentInfected)

mod1.0 <-glmmTMB(PercentInfected~(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.0)

mod1.1 <-glmmTMB(PercentInfected~Site + (1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.1)

mod1.2<-glmmTMB(PercentInfected~Site+PDO_Oct_Mar+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.2)

mod1.3<-glmmTMB(PercentInfected~Site+PDO_April_Sept_lagged+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.3)

mod1.4<-glmmTMB(PercentInfected~Site+NPGO_April_Sept_lagged+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.4)

mod1.5<-glmmTMB(PercentInfected~Site+NPGO_Oct_Mar+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.5)

mod1.6<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.6)

mod1.7<-glmmTMB(PercentInfected~Cohort+Site+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.7)

mod1.8<-glmmTMB(PercentInfected~Site+Cohort+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.8)

mod1.9<-glmmTMB(PercentInfected~Site+MatureMillions+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula=~1+(1|UniqueSet), family=beta_family())
summary(mod1.9)

mod1.10<-glmmTMB(PercentInfected~Site+RecruitmentMillions+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.10)

mod1.11<-glmmTMB(PercentInfected~OffshoreUpwelling_Oct_Mar+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.11)

mod1.12<-glmmTMB(PercentInfected~OffshoreUpwelling_April_Sept_lagged+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.12)

mod1.13<-glmmTMB(PercentInfected~AssignedAge+OffshoreUpwelling_April_Sept_lagged+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.13)

mod1.14<-glmmTMB(PercentInfected~Site+AssignedAge+OffshoreUpwelling_April_Sept_lagged+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.14)

mod1.14<-glmmTMB(PercentInfected~Site+AssignedAge+OffshoreUpwelling_April_Sept_lagged+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.14)

mod1.15<-glmmTMB(PercentInfected~Site+AssignedAge+Cohort+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.15)

mod1.16<-glmmTMB(PercentInfected~Site+AssignedAge+Year+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~1+(1|UniqueSet), family=beta_family())
summary(mod1.16)


#
require(AICcmodavg)
AICc(mod1.0)
AICc(mod1.1)
AICc(mod1.2)
AICc(mod1.3)
AICc(mod1.4)
AICc(mod1.5)
AICc(mod1.6)## BEST
AICc(mod1.7)
AICc(mod1.8)
AICc(mod1.8)
AICc(mod1.9)
AICc(mod1.10)
AICc(mod1.11)
AICc(mod1.12)
AICc(mod1.13)
AICc(mod1.14)
AICc(mod1.15)
AICc(mod1.16)

##Now take best conditional model and evaluate zero-inflated models

mod1.6a<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~Site+(1|UniqueSet), family=beta_family())
summary(mod1.6a)

mod1.6b<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~AssignedAge+(1|UniqueSet), family=beta_family())
summary(mod1.6b)

mod1.6c<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~MatureMillions+(1|UniqueSet), family=beta_family())
summary(mod1.6c)

mod1.6d<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~RecruitmentMillions+(1|UniqueSet), family=beta_family())
summary(mod1.6d)

mod1.6e<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~PDO_Oct_Mar+(1|UniqueSet), family=beta_family())
summary(mod1.6e)

mod1.6f<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~PDO_April_Sept_lagged+(1|UniqueSet), family=beta_family())
summary(mod1.6e)

mod1.6g<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~NPGO_Oct_Mar+(1|UniqueSet), family=beta_family())
summary(mod1.6g)

mod1.6h<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~NPGO_April_Sept_lagged+(1|UniqueSet), family=beta_family())
summary(mod1.6h)

mod1.6i<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~OffshoreUpwelling_Oct_Mar+(1|UniqueSet), family=beta_family())
summary(mod1.6i)

mod1.6j<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~OffshoreUpwelling_April_Sept_lagged+(1|UniqueSet), family=beta_family())
summary(mod1.6j)

mod1.6k<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~Cohort+(1|UniqueSet), family=beta_family())
summary(mod1.6k)

mod1.6l<-glmmTMB(PercentInfected~Site+AssignedAge+(1|UniqueSet), data=AllHerringALInfected, na.action=na.omit, ziformula= ~Year+(1|UniqueSet), family=beta_family())
summary(mod1.6l)



AICc(mod1.6) ##Best model
AICc(mod1.6a)
AICc(mod1.6b)
AICc(mod1.6c)
AICc(mod1.6d)
AICc(mod1.6e)
AICc(mod1.6f)
AICc(mod1.6g)
AICc(mod1.6h)
AICc(mod1.6i)
AICc(mod1.6j)
AICc(mod1.6k)
AICc(mod1.6l)


############################################################################################
############################################################################################
#Model selection for prevalence data


mod2.0<-glmer(IchFinal~(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.0)

mod2.1<-glmer(IchFinal~Site+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.1)

mod2.2<-glmer(IchFinal~Site+AssignedAge+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.2)

mod2.3<-glmer(IchFinal~Site+AssignedAge+RecruitmentMillions+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.3)

mod2.4<-glmer(IchFinal~Site+AssignedAge+MatureMillions+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.4)

mod2.5<-glmer(IchFinal~Site+AssignedAge+PDO_Oct_Mar+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.5)

mod2.6<-glmer(IchFinal~Site+AssignedAge+PDO_April_Sept_lagged+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.6)

mod2.7<-glmer(IchFinal~Site+AssignedAge+NPGO_Oct_Mar+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.7)

mod2.8<-glmer(IchFinal~Site+AssignedAge+NPGO_April_Sept_lagged+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.8)

mod2.9<-glmer(IchFinal~Site+AssignedAge+OffshoreUpwelling_Oct_Mar+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.9)

mod2.10<-glmer(IchFinal~Site+AssignedAge+OffshoreUpwelling_April_Sept_lagged+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.10)

mod2.11<-glmer(IchFinal~Site+AssignedAge+Cohort+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.11)

mod2.12<-glmer(IchFinal~Site+AssignedAge+Year+(1|UniqueSet), data=AllHerringAL, family=binomial)
summary(mod2.12)


AICc(mod2.0)
AICc(mod2.1)
AICc(mod2.2) 
AICc(mod2.3)
AICc(mod2.4)
AICc(mod2.5)
AICc(mod2.6)
AICc(mod2.7)
AICc(mod2.8)
AICc(mod2.9) 
AICc(mod2.10)
AICc(mod2.11)
AICc(mod2.12)


#Predictions for best model

newdata <- with(AllHerringAL, expand.grid(AssignedAge=seq(3,7,1), PDO_Oct_Mar=seq(-1.5, 2, .3), Site=c('Sitka', 'PWS')))
predicted <- predict(mod2.5, newdata, type='response', re.form=NA)    # new data, level-0
newdata<-cbind(newdata, predicted)
newdata$Location<-ifelse(newdata$Site=='PWS', 'Prince William Sound', 'Sitka Sound')

require(plyr)
newdata2<-ddply(newdata, .(Site), summarize, meanPrev=mean(predicted))


p<-ggplot(newdata, aes(x=AssignedAge, y=predicted, color=as.factor(PDO_Oct_Mar)))+geom_line(size=2.5)+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p<-p+xlab("Age")+ylab("Predicted prevalence")+ scale_color_discrete(l=45, name = "PDO Index \n(October - March)", guide = guide_legend(reverse = TRUE))
p+facet_grid(.~Location) #+ scale_color_discrete(low = "red", high = "blue", na.value = NA)

