

require(readr)
require(lubridate)
#Read in the data, available from ADFG, contact Sherri Dressel (sherri.dressel@alaska.gov)

#SitkaAgeLength <- read_csv("Sitka ALEX cast net for R.csv")#Aged samples over time
SitkaAgeLength$SAMPLE_DATE<-as.Date(SitkaAgeLength$SAMPLE_DATE, "%m/%d/%y")
SitkaAgeLength$Month<-month(SitkaAgeLength$SAMPLE_DATE)

summary(SitkaAgeLength$Month) #looks like all sampling is in March and April. This is good!

SitkaAgeLength<-SitkaAgeLength[!is.na(SitkaAgeLength$LENGTH_MM),]
SitkaAgeLength<-SitkaAgeLength[!is.na(SitkaAgeLength$AGE),]
SitkaAgeLength<-SitkaAgeLength[SitkaAgeLength$LENGTH_MM!=0,]

#Calculate Fork lengths (data are in Standard lengths)
SitkaAgeLength$Fork_Length<--2.348+1.060*SitkaAgeLength$LENGTH_MM

#Set aside random 10% of dataset for validation
SitkaAgeLength$Random<-runif(n=length(SitkaAgeLength$Month))
SitkaAgeLengthSubset<-SitkaAgeLength[SitkaAgeLength$Random>.9,]
SitkaAgeLength<-SitkaAgeLength[SitkaAgeLength$Random<.9,]

#get mean and sd of length for each age by year
require(plyr)
AgeLengthSummary<-ddply(SitkaAgeLength, .(YEAR, AGE), summarize, LengthMean=mean(Fork_Length), LengthSD=sqrt(var(Fork_Length)),n=length(Fork_Length))
AgeLengthSummary<-AgeLengthSummary[AgeLengthSummary$AGE>2,] #limit ages to between 2 and 17
AgeLengthSummary<-AgeLengthSummary[AgeLengthSummary$AGE<17,]
AgeLengthSummary<-AgeLengthSummary[AgeLengthSummary$n>10,] #need a sample size greater than 10 to get calculate a mean and sd
AgeLengthSummary$AGE<-as.factor(AgeLengthSummary$AGE)
AgeLengthSummary$YEAR<-as.factor(AgeLengthSummary$YEAR)


require(ggplot2)
require(RColorBrewer)
require(colorspace)

AgeLengthSummary$AGE<-as.factor(AgeLengthSummary$AGE)

p0<-ggplot(AgeLengthSummary, aes(YEAR,LengthMean, fill=AGE), col="black")+geom_col(position="dodge", color="black")+ scale_fill_brewer()+theme_dark()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))
p0 



#Remove fish with an age greater than 17 (Error codes are 18 and 19)
SitkaAgeLength<-SitkaAgeLength[SitkaAgeLength$AGE<17,]
SitkaAgeLength<-SitkaAgeLength[SitkaAgeLength$Fork_Length<400,]

#create a plus group
SitkaAgeLength$AGE<-ifelse(SitkaAgeLength$AGE>7, 7, SitkaAgeLength$AGE)
require(fishmethods)


#calculate alk for each year
x<-SitkaAgeLength[SitkaAgeLength$YEAR==2007,]
alk_Sitka<-with(x,alk(age=round(AGE,0),size=Fork_Length, binsize=10))
alk_Sitka$year<-2007


for (i in 2008:2019){
  x<-SitkaAgeLength[SitkaAgeLength$YEAR==i,]
  y<-with(x,alk(age=round(AGE,0),size=Fork_Length,binsize=10))
  y$year<-i
  alk_Sitka<-rbind.fill(alk_Sitka, y)
}

alk_Sitka$A2_prob<-alk_Sitka$A2/alk_Sitka$nl
alk_Sitka$A3_prob<-alk_Sitka$A3/alk_Sitka$nl
alk_Sitka$A4_prob<-alk_Sitka$A4/alk_Sitka$nl
alk_Sitka$A5_prob<-alk_Sitka$A5/alk_Sitka$nl
alk_Sitka$A6_prob<-alk_Sitka$A6/alk_Sitka$nl
alk_Sitka$A7_prob<-alk_Sitka$A7/alk_Sitka$nl

age_length_key_Sitka_Fishmethods<-subset(alk_Sitka, select = -c(A2, A3, A4, A5, A6, A7) )
age_length_key_Sitka_Fishmethods$A2_prob<-ifelse(is.na(age_length_key_Sitka_Fishmethods$A2_prob), 0, age_length_key_Sitka_Fishmethods$A2_prob)
write.table(age_length_key_Sitka_Fishmethods, "c:/Users/mayag/Dropbox/Ichthyophonus/Age-length Data/alk_Sitka_fishmethods.txt", sep="\t")


#######################################################################
#######################################################################
#Check table against 10% chosen at random
names(alk_Sitka)[names(alk_Sitka) == "len"] <- "Fork_Length_Bin"
names(alk_Sitka)[names(alk_Sitka) == "year"] <- "YEAR"

SitkaAgeLengthSubset$Fork_Length_Bin<-floor(SitkaAgeLengthSubset$Fork_Length/10)*10+5
SitkaAgeLengthSubset<-merge(SitkaAgeLengthSubset, alk_Sitka, by=c("YEAR", "Fork_Length_Bin"), all.x=TRUE)

#Create minus group for fish 3 and younger (not many under 3 in the dataset)
SitkaAgeLengthSubset$A3_prob_minus<-ifelse(is.na(SitkaAgeLengthSubset$A2_prob),SitkaAgeLengthSubset$A3_prob, SitkaAgeLengthSubset$A2_prob+SitkaAgeLengthSubset$A3_prob)


#Assign an age to each sample
SitkaAgeLengthSubset$Age34<-SitkaAgeLengthSubset$A3_prob_minus+SitkaAgeLengthSubset$A4_prob
SitkaAgeLengthSubset$Age345<-SitkaAgeLengthSubset$A3_prob_minus+SitkaAgeLengthSubset$A4_prob+SitkaAgeLengthSubset$A5_prob
SitkaAgeLengthSubset$Age3456<-SitkaAgeLengthSubset$A3_prob_minus+SitkaAgeLengthSubset$A4_prob+SitkaAgeLengthSubset$A5_prob+SitkaAgeLengthSubset$A6_prob
SitkaAgeLengthSubset$Age34567plus<-SitkaAgeLengthSubset$A3_prob_minus+SitkaAgeLengthSubset$A4_prob+SitkaAgeLengthSubset$A5_prob+SitkaAgeLengthSubset$A6_prob+SitkaAgeLengthSubset$A7_prob

SitkaAgeLengthSubset$Random<-runif(length(SitkaAgeLengthSubset$A3_prob_minus))
SitkaAgeLengthSubset$AssignedAge<-ifelse(is.na(SitkaAgeLengthSubset$Age34567plus),NA, 
                                           ifelse(SitkaAgeLengthSubset$Random<SitkaAgeLengthSubset$A3_prob_minus, 3,
                                                  ifelse(SitkaAgeLengthSubset$Random<SitkaAgeLengthSubset$Age34, 4,
                                                         ifelse(SitkaAgeLengthSubset$Random<SitkaAgeLengthSubset$Age345, 5,
                                                                ifelse(SitkaAgeLengthSubset$Random<SitkaAgeLengthSubset$Age3456, 6, 7)))))

check<-ddply(SitkaAgeLengthSubset, .(YEAR, AssignedAge), summarize, Random=mean(Random))
ggplot(check, aes(x=AssignedAge, y=Random, colour=as.factor(YEAR)))+geom_point()+geom_line()
mean(SitkaAgeLengthSubset$Random)


############################################
##Compare predicted ages with actual ages in PWS
AgeComp<-SitkaAgeLengthSubset[!is.na(SitkaAgeLengthSubset$AGE),]
AgeComp<-AgeComp[AgeComp$AGE>2,]
AgeComp<-AgeComp[AgeComp$AGE<10,]

AgeComp$AGE<-ifelse(AgeComp$AGE>7, 7, AgeComp$AGE)
AgeComp$AgeDiff<-AgeComp$AGE-AgeComp$AssignedAge
mean(AgeComp$AgeDiff, na.rm=TRUE) ##-1.053
sd(AgeComp$AgeDiff, na.rm=TRUE)/sqrt(length(AgeComp$AgeDiff))


AgeCompSum<-ddply(AgeComp, .(YEAR, AGE), summarize, meanDiff=mean(AgeDiff, na.rm=TRUE), SampleSize=length(AgeDiff), meanDiffSE<-sd(AgeDiff, na.rm=TRUE)/sqrt(length(AgeDiff)))

ggplot(AgeCompSum, aes(x=AGE, y=meanDiff, color=YEAR, size=SampleSize))+geom_point()+geom_line()+facet_grid(YEAR~.)+ geom_hline(yintercept=0, linetype="dashed", color = "red")+ylab("Scale age - ALK Age")+xlab('Scale Age')



##bland altman comparison of methods
library(blandr)
library(ggplot2)

comp<-blandr.statistics ( AgeComp$AGE, AgeComp$AssignedAge , sig.level=0.95 )

blandr.output.text ( AgeComp$AGE, AgeComp$AssignedAge , sig.level=0.95 )

wright.stats <- blandr.statistics( AgeComp$AGE, AgeComp$AssignedAge )
wright.plot <- blandr.plot.limits (wright.stats )
plot( x=wright.stats$means , y=wright.stats$differences )


