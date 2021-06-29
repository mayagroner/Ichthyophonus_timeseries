require(readr)

#Import and combine datasets for 2007-2019
#data available here: https://researchworkspace.com/file/2532988/df35b.276.1-pwsHerrASWL.csv.csv
early_Cordova <- read_csv("C:/Users/mayag/Dropbox/Ichthyophonus/Age-length Data/1973-2014 Cordova Formatted to merge.csv")#aged herring from PWS (data from Stormy)
early_Cordova<-as.data.frame(early_Cordova)
early_Cordova$Year<-as.numeric(early_Cordova$Year)
early_Cordova<-early_Cordova[early_Cordova$Year>2006,]
early_Cordova<-early_Cordova[early_Cordova$Year<2014,]

#data available here: https://researchworkspace.com/file/6327437/PWS_HAWL_Data_2014-2019.csv
late_Cordova<- read_csv("C:/Users/mayag/Dropbox/Ichthyophonus/Age-length Data/PWS_HAWL_Data_2014-2019.csv")

colnames(early_Cordova)
colnames(late_Cordova)
CordovaAgeLength<-rbind(early_Cordova, late_Cordova)

summary(CordovaAgeLength$Month)
summary(CordovaAgeLength$Year)

#restrict data to spawning time (months 3 and 4)
CordovaAgeLength<-CordovaAgeLength[CordovaAgeLength$Month<5,]
CordovaAgeLength<-CordovaAgeLength[CordovaAgeLength$Month>2,]

#remove data where age or length is NA
CordovaAgeLength<-CordovaAgeLength[!is.na(CordovaAgeLength$LENGTH),]
CordovaAgeLength<-CordovaAgeLength[!is.na(CordovaAgeLength$AGE),]
CordovaAgeLength<-CordovaAgeLength[CordovaAgeLength$LENGTH!=0,]


#Calculate Fork lengths (data are in Standard lengths)- Kristen has a dataset, which I used to calculate the conversion
CordovaAgeLength$Fork_Length<- -2.348+1.060*CordovaAgeLength$LENGTH                # -2.348+1.060*CordovaAgeLength$LENGTH

CordovaAgeLength$Random<-runif(n=length(CordovaAgeLength$Year))
CordovaAgeLengthSubset<-CordovaAgeLength[CordovaAgeLength$Random>.9,]
CordovaAgeLength<-CordovaAgeLength[CordovaAgeLength$Random<.9,]

#get mean and sd of length for each age by year
require(plyr)
AgeLengthSummary<-ddply(CordovaAgeLength, .(Year, SEX_code, AGE), summarize, LengthMean=mean(Fork_Length), LengthSD=sqrt(var(Fork_Length)),n=length(LENGTH))
AgeLengthSummary<-AgeLengthSummary[AgeLengthSummary$AGE>2,] #limit ages to between 2 and 11
AgeLengthSummary<-AgeLengthSummary[AgeLengthSummary$AGE<11,]
AgeLengthSummary<-AgeLengthSummary[AgeLengthSummary$n>10,] #need a sample size greater than 10 to get calculate a mean and sd
AgeLengthSummary$AGE<-as.factor(AgeLengthSummary$AGE)
AgeLengthSummary$Year<-as.factor(AgeLengthSummary$Year)


require(ggplot2)
require(RColorBrewer)
require(colorspace)

AgeLengthSummary$AGE<-as.factor(AgeLengthSummary$AGE)

p0<-ggplot(AgeLengthSummary, aes(Year,LengthMean, fill=AGE), col="black")+geom_col(position="dodge", color="black")+ scale_fill_brewer()+theme_dark()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p0+facet_grid(SEX_code~.) 

p1<-ggplot(AgeLengthSummary, aes(AGE,LengthMean, fill=Year))+geom_col(position="dodge", color="black")+scale_fill_brewer()+theme_dark()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1+facet_grid(SEX_code~.)

AgeLengthSummary$AGE<-as.numeric(AgeLengthSummary$AGE)
p2<-ggplot(AgeLengthSummary, aes(AGE,LengthMean, color=Year))+geom_point()+geom_line(aes(lty=as.factor(SEX_code)))+theme_dark()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"))                                                                                        
p2

#Length of SEX_code '2' (I think these are females) is > than 1 Ideally, should have different age length keys for each sex code. 

#Need to convert standard length to fork length

#Fork length = -2.348+1.06*Std Length: Valid for fish with standard lengths between 120 and 242 mm



#Remove fish with an age greater than 15 (Suspicious outliers)
CordovaAgeLength<-CordovaAgeLength[CordovaAgeLength$AGE<15,]
CordovaAgeLength<-CordovaAgeLength[CordovaAgeLength$Fork_Length<400,]

#create a plus group
CordovaAgeLength$AGE<-ifelse(CordovaAgeLength$AGE>7, 7, CordovaAgeLength$AGE)
require(fishmethods)

x<-CordovaAgeLength[CordovaAgeLength$Year==2007,]
alk_Cordova<-with(x,alk(age=round(AGE,0),size=Fork_Length,binsize=10))
alk_Cordova$year<-2007


for (i in 2008:2019){
  x<-CordovaAgeLength[CordovaAgeLength$Year==i,]
  y<-with(x,alk(age=round(AGE,0),size=Fork_Length,binsize=10))
  y$year<-i
  alk_Cordova<-rbind.fill(alk_Cordova, y)
}


alk_Cordova$A1_prob<-alk_Cordova$A1/alk_Cordova$nl
alk_Cordova$A2_prob<-alk_Cordova$A2/alk_Cordova$nl
alk_Cordova$A3_prob<-alk_Cordova$A3/alk_Cordova$nl
alk_Cordova$A4_prob<-alk_Cordova$A4/alk_Cordova$nl
alk_Cordova$A5_prob<-alk_Cordova$A5/alk_Cordova$nl
alk_Cordova$A6_prob<-alk_Cordova$A6/alk_Cordova$nl
alk_Cordova$A7_prob<-alk_Cordova$A7/alk_Cordova$nl

age_length_key_Cordova_Fishmethods<-subset(alk_Cordova, select = -c(A1, A2, A3, A4, A5, A6, A7) )
#write.table(age_length_key_Cordova_Fishmethods, "c:/Users/mayag/Dropbox/Ichthyophonus/Age-length Data/alk_Cordova_fishmethods.txt", sep="\t")


#######################################################################
#######################################################################
#Check table against 10% chosen at random
names(alk_Cordova)[names(alk_Cordova) == "len"] <- "Fork_Length_Bin"
names(alk_Cordova)[names(alk_Cordova) == "year"] <- "Year"

CordovaAgeLengthSubset$Fork_Length_Bin<-floor(CordovaAgeLengthSubset$Fork_Length/10)*10+5
CordovaAgeLengthSubset<-merge(CordovaAgeLengthSubset, alk_Cordova, by=c("Year", "Fork_Length_Bin"), all.x=TRUE)

CordovaAgeLengthSubset$A1_prob<-ifelse(is.na(CordovaAgeLengthSubset$A1_prob), 0, CordovaAgeLengthSubset$A1_prob)
CordovaAgeLengthSubset$A3_prob_minus<-CordovaAgeLengthSubset$A1_prob+CordovaAgeLengthSubset$A2_prob+CordovaAgeLengthSubset$A3_prob

#Assign an age to each sample
CordovaAgeLengthSubset$Age34<-CordovaAgeLengthSubset$A3_prob_minus+CordovaAgeLengthSubset$A4_prob
CordovaAgeLengthSubset$Age345<-CordovaAgeLengthSubset$A3_prob_minus+CordovaAgeLengthSubset$A4_prob+CordovaAgeLengthSubset$A5_prob
CordovaAgeLengthSubset$Age3456<-CordovaAgeLengthSubset$A3_prob_minus+CordovaAgeLengthSubset$A4_prob+CordovaAgeLengthSubset$A5_prob+CordovaAgeLengthSubset$A6_prob
CordovaAgeLengthSubset$Age34567plus<-CordovaAgeLengthSubset$A3_prob_minus+CordovaAgeLengthSubset$A4_prob+CordovaAgeLengthSubset$A5_prob+CordovaAgeLengthSubset$A6_prob+CordovaAgeLengthSubset$A7_prob

CordovaAgeLengthSubset$Random<-runif(length(CordovaAgeLengthSubset$A3_prob_minus))
CordovaAgeLengthSubset$AssignedAge<-ifelse(is.na(CordovaAgeLengthSubset$Age34567plus),NA, 
                                 ifelse(CordovaAgeLengthSubset$Random<CordovaAgeLengthSubset$A3_prob_minus, 3,
                                        ifelse(CordovaAgeLengthSubset$Random<CordovaAgeLengthSubset$Age34, 4,
                                               ifelse(CordovaAgeLengthSubset$Random<CordovaAgeLengthSubset$Age345, 5,
                                                      ifelse(CordovaAgeLengthSubset$Random<CordovaAgeLengthSubset$Age3456, 6, 7)))))

check<-ddply(CordovaAgeLengthSubset, .(Year, AssignedAge), summarize, Random=mean(Random))
ggplot(check, aes(x=AssignedAge, y=Random, colour=as.factor(Year)))+geom_point()+geom_line()
mean(CordovaAgeLengthSubset$Random)
############################################
##Compare predicted ages with actual ages in PWS
AgeComp<-CordovaAgeLengthSubset[!is.na(CordovaAgeLengthSubset$AGE),]
AgeComp<-AgeComp[AgeComp$AGE>2,]
AgeComp<-AgeComp[AgeComp$AGE<10,]

AgeComp$AGE<-ifelse(AgeComp$AGE>7, 7, AgeComp$AGE)
AgeComp$AgeDiff<-AgeComp$AGE-AgeComp$AssignedAge
mean(AgeComp$AgeDiff, na.rm=TRUE) ##-1.053
sd(AgeComp$AgeDiff, na.rm=TRUE)/sqrt(length(AgeComp$AgeDiff))


AgeCompSum<-ddply(AgeComp, .(Year, AGE), summarize, meanDiff=mean(AgeDiff, na.rm=TRUE), SampleSize=length(AgeDiff), meanDiffSE<-sd(AgeDiff, na.rm=TRUE)/sqrt(length(AgeDiff)))

ggplot(AgeCompSum, aes(x=AGE, y=meanDiff, color=Year, size=SampleSize))+geom_point()+geom_line()+facet_grid(Year~.)+ geom_hline(yintercept=0, linetype="dashed", color = "red")+ylab("Scale age - ALK Age")+xlab('Scale Age')



##bland altman comparison
library(blandr)
library(ggplot2)

comp<-blandr.statistics ( AgeComp$AGE, AgeComp$AssignedAge , sig.level=0.95 )

blandr.output.text ( AgeComp$AGE, AgeComp$AssignedAge , sig.level=0.95 )

wright.stats <- blandr.statistics( AgeComp$AGE, AgeComp$AssignedAge )
wright.plot <- blandr.plot.limits (wright.stats )
plot( x=wright.stats$means , y=wright.stats$differences )



