#Calculate yearly prevalence for each site
require(readr)
Prev_at_age <- read_csv("Prevalence_age_site.csv")
Mature_at_age <- read_csv("Mature_prop_at_age_PWS_Sitka.csv")

PrevMat<-merge(Prev_at_age, Mature_at_age, by=c('Site', 'Year', 'Age'))

PrevMat$PrevalenceSD<-sqrt(PrevMat$Prevalence*(1-PrevMat$Prevalence)/PrevMat$SampleSize)
PrevMat$PropPrev<-PrevMat$Prevalence*PrevMat$Prop_at_age
PrevMat$PropPrevSD<-sqrt(PrevMat$PropPrev*(1-PrevMat$PropPrev)/PrevMat$SampleSize)
PrevMat$PropPrevCI<-sqrt(PrevMat$PropPrev*(1-PrevMat$PropPrev)/PrevMat$SampleSize)*1.96
PrevMat$PropPrevCI_sq<-PrevMat$PropPrevCI*PrevMat$PropPrevCI
PrevMat$Site<-ifelse(PrevMat$Site=='PWS', 'Prince William Sound', 'Sitka Sound')

require(plyr)
require(ggplot2)


##Graph Age-adjusted prevalence
PrevMat_Sum<-ddply(PrevMat, .(Site, Year), summarize, Prevalence=sum(PropPrev, na.rm=TRUE), PrevCI=(sqrt(sum(PropPrevCI_sq, na.rm=TRUE))))

PrevMat_Sum$PrevSE<-PrevMat_Sum$PrevCI/1.96

p<-ggplot(PrevMat_Sum, aes(x=Year, y=Prevalence, color=Site))+geom_point()+geom_line(size=1)+ylim(0, 0.65)+
  geom_ribbon(aes(ymin=Prevalence-PrevSE, ymax=Prevalence+PrevSE, fill=Site), colour=NA, alpha=0.1)

p<-p+theme_bw()+ylab("Age-based Ichthyophonus sp. prevalence \nestimate for mature Pacific herring")+ theme(text = element_text(size = 8)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_smooth(method='lm', se=FALSE, linetype=2)+ theme(legend.position='none' )+ labs(tag = "B")+
  theme(plot.tag = element_text(),
        plot.tag.position = c(0.1, 0.94))

p

##Graph sample prevalence

PrevMat_Sum_Sample<-ddply(PrevMat, .(Site, Year), summarize, PrevalenceSamp=sum(Prevalence*Age, na.rm=TRUE)/sum(Age), SS=sum(SampleSize))
PrevMat_Sum_Sample$PrevalenceSampSE<-sqrt((PrevMat_Sum_Sample$PrevalenceSamp*(1-PrevMat_Sum_Sample$PrevalenceSamp))/PrevMat_Sum_Sample$SS)


p3<-ggplot(PrevMat_Sum_Sample, aes(x=Year, y=PrevalenceSamp, color=Site))+geom_point()+geom_line(size=1)+
  ylim(0, 0.65)+geom_ribbon(aes(ymin=PrevalenceSamp-PrevalenceSampSE, ymax=PrevalenceSamp+PrevalenceSampSE, fill=Site), colour=NA, alpha=0.1)

p3<-p3+theme_bw()+ylab("Ichthyophonus sp. prevalence in mature \nPacific herring sample ")+ theme(text = element_text(size = 8)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_smooth(method='lm', se=FALSE, linetype=2)+
  theme(legend.position = 'none')+ labs(tag = "A")+
  theme(plot.tag = element_text(),
        plot.tag.position = c(0.1, 0.94))

p3
#Graph proportion of Prevalent that are symptomatic

y_title2 <- expression(paste("Proportion of Pacific herring infected with ", italic("Ichthyophonus sp.")," that are symptomatic"))
PrevMat$PropSymptomatic<-ifelse(PrevMat$SPrevalence==0, 0, PrevMat$SPrevalence/(PrevMat$Prevalence-PrevMat$ASUnknown)*PrevMat$Prop_at_age) 

PrevMatS_Sum<-ddply(PrevMat, .(Site, Year), summarize, PrevalenceS=sum(PropSymptomatic, na.rm=TRUE), PrevCI=(sqrt(sum(PropPrevCI_sq))))
PrevMatS_Sum$PrevalenceS<-ifelse(PrevMatS_Sum$PrevalenceS==0, NA, PrevMatS_Sum$PrevalenceS)#remove 0s, these are cases where symptoms were not quantified
PrevMatS_Sum$PrevSE<-PrevMatS_Sum$PrevCI/1.96
PrevMatS_Sum$LowerSE<-ifelse(PrevMatS_Sum$PrevalenceS-PrevMatS_Sum$PrevSE==0, 0, PrevMatS_Sum$PrevalenceS-PrevMatS_Sum$PrevSE)

PrevMatS_Sum<-PrevMatS_Sum[!is.na(PrevMatS_Sum$PrevalenceS),]
p1<-ggplot(PrevMatS_Sum, aes(x=Year, y=PrevalenceS, color=Site))+geom_line(size=1)+geom_ribbon(aes(ymin=LowerSE, ymax=PrevalenceS+PrevSE, fill=Site), colour=NA, alpha=0.1)

p1<-p1+theme_bw()+ylab("Proportion of mature infected Pacific herring \nwith high Ichthyophonus sp. infection") + theme(text=element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_point()+geom_smooth(method='lm', se=FALSE, linetype=2)+ theme(legend.position='none')+ scale_x_continuous(breaks=seq(2009,2019,2))+ labs(tag = "D")+
  theme(plot.tag = element_text(),
        plot.tag.position = c(0.1, 0.94))+ scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))
p1
#Calculate mean age for each year
Mature_at_age$WeightedAge<-Mature_at_age$Age*Mature_at_age$Prop_at_age

Mature_at_age_Sum<-ddply(Mature_at_age, .(Year, Site), summarize, MeanAge=sum(WeightedAge))

Mature_at_age_Sum$Site<-ifelse(Mature_at_age_Sum$Site=='PWS', 'Prince William Sound', 'Sitka Sound')

PrevMat_Sum<-merge(PrevMat_Sum, Mature_at_age_Sum, by=c('Year', 'Site'))
Mature_at_age_Sum<-Mature_at_age_Sum[!is.na(Mature_at_age_Sum$Year),]

p2<-ggplot(PrevMat_Sum, aes(x=MeanAge, y=Prevalence, color=Site))+geom_point(size=3)+geom_smooth(method='lm', aes(fill=Site), alpha=.1, linetype=2)+xlab('Mean age of mature population')+ylab("Age-based Ichthyophonus sp. prevalence \nin mature Pacific herring")+ theme(legend.position='none')
p2<-p2+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position='none')
p2<-p2+geom_errorbar(aes(ymin=ifelse(Prevalence-PrevSE<0, 0, Prevalence-PrevSE), ymax=Prevalence+PrevSE), width=.1,
                     position=position_dodge(0.05))+ labs(tag = "C")+
  theme(plot.tag = element_text(),
        plot.tag.position = c(0.1, 0.94), text=element_text(size = 8),)
p2

require(gridExtra)
grid.arrange(p3, p, p2, p1, ncol = 1)
