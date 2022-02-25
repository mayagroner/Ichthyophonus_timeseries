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

y_title <- expression(paste(italic("Ichthyophonus sp."), " prevalence in mature Pacific herring"))

PrevMat_Sum<-ddply(PrevMat, .(Site, Year), summarize, Prevalence=sum(PropPrev), PrevCI=(sqrt(sum(PropPrevCI_sq))))

PrevMat_Sum$PrevSE<-PrevMat_Sum$PrevCI/sqrt(length(PrevMat_Sum$PrevCI))
PrevMat_Sum$LowerSE<-ifelse(PrevMat_Sum$Prevalence-PrevMat_Sum$PrevSE<0, 0, PrevMat_Sum$Prevalence-PrevMat_Sum$PrevSE)

p<-ggplot(PrevMat_Sum, aes(x=Year, y=Prevalence, color=Site))+geom_point()+geom_line(size=1)+ylim(0, 0.65)+geom_ribbon(aes(ymin=LowerSE, ymax=Prevalence+PrevSE, fill=Site), linetype=2, alpha=0.1)

p<-p+theme_bw()+ylab("Ichthyophonus sp. prevalence in mature \nPacific herring")+ theme(text = element_text(size = 8)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_smooth(method='lm', se=FALSE, linetype=2)+ theme(legend.position = c(0.8, 0.8))+ theme(legend.title = element_blank())+ labs(tag = "A")+
  theme(plot.tag = element_text(),
        plot.tag.position = c(0.1, 0.94))

ddply(PrevMat_Sum, .(Site), summarize, Prevalence=mean(Prevalence))

#Graph proportion of Prevalent that are symptomatic

y_title2 <- expression(paste("Proportion of Pacific herring infected with ", italic("Ichthyophonus sp.")," that are symptomatic"))
PrevMat$PropSymptomatic<-ifelse(PrevMat$SPrevalence==0, 0, PrevMat$SPrevalence/(PrevMat$Prevalence-PrevMat$ASUnknown)*PrevMat$Prop_at_age) 

PrevMatS_Sum<-ddply(PrevMat, .(Site, Year), summarize, PrevalenceS=sum(PropSymptomatic, na.rm=TRUE), PrevalenceU=sum(ASUnknown, na.rm=TRUE), PrevCI=(sqrt(sum(PropPrevCI_sq))))
PrevMatS_Sum$PrevalenceS<-ifelse(PrevMatS_Sum$PrevalenceS==0, NA, PrevMatS_Sum$PrevalenceS)#remove 0s, these are cases where symptoms were not quantified
PrevMatS_Sum$PrevSE<-PrevMatS_Sum$PrevCI/1.96
PrevMatS_Sum$LowerSE<-ifelse(PrevMatS_Sum$PrevalenceS-PrevMatS_Sum$PrevSE==0, 0, PrevMatS_Sum$PrevalenceS-PrevMatS_Sum$PrevSE)

PrevMatS_Sum<-PrevMatS_Sum[!is.na(PrevMatS_Sum$PrevalenceS),]
p1<-ggplot(PrevMatS_Sum, aes(x=Year, y=PrevalenceS, color=Site))+geom_line(size=1)+geom_ribbon(aes(ymin=LowerSE, ymax=PrevalenceS+PrevSE, fill=Site), linetype=2, alpha=0.1)

p1<-p1+theme_bw()+ylab("Proportion of mature infected Pacific herring with \nclinical signs of Ichthyophonus sp.") + theme(text=element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_point()+geom_smooth(method='lm', se=FALSE, linetype=2)+ theme(legend.position='none')+ scale_x_continuous(breaks=seq(2009,2019,2))+ labs(tag = "C")+
  theme(plot.tag = element_text(),
        plot.tag.position = c(0.1, 0.94))+ scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2))
p1
#Calculate mean age for each year
Mature_at_age$WeightedAge<-Mature_at_age$Age*Mature_at_age$Prop_at_age

Mature_at_age_Sum<-ddply(Mature_at_age, .(Year, Site), summarize, MeanAge=sum(WeightedAge))

Mature_at_age_Sum$Site<-ifelse(Mature_at_age_Sum$Site=='PWS', 'Prince William Sound', 'Sitka Sound')
  
PrevMat_Sum<-merge(PrevMat_Sum, Mature_at_age_Sum, by=c('Year', 'Site'))
Mature_at_age_Sum<-Mature_at_age_Sum[!is.na(Mature_at_age_Sum$Year),]

p2<-ggplot(PrevMat_Sum, aes(x=MeanAge, y=Prevalence, color=Site))+geom_point(size=3)+geom_smooth(method='lm', aes(fill=Site), alpha=.1, linetype=2)+xlab('Mean age of mature population')+ylab("Ichthyophonus sp. prevalence in mature \nPacific herring")+ theme(legend.position='none')
p2<-p2+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ylim(0, 0.61)+ theme(legend.position='none')
p2<-p2+geom_errorbar(aes(ymin=LowerSE, ymax=Prevalence+PrevSE), width=.1,
                position=position_dodge(0.05))+ labs(tag = "B")+
      theme(plot.tag = element_text(),
      plot.tag.position = c(0.1, 0.94), text=element_text(size = 8),)
p2

require(gridExtra)
grid.arrange(p, p2, p1, ncol = 1)


##Weighted linear regressions of age-prevalence relationships for each year class for each site
PrevMat$Weight<-1/(PrevMat$PrevalenceSD^2)
PrevMat<-PrevMat[PrevMat$Weight!=Inf,]
PrevMat<-PrevMat[!is.na(PrevMat$Site),]
PrevMat$Cohort<-as.factor(PrevMat$Cohort)
model<-lm(Prevalence~Site*Age*Cohort, weight=Weight, data=PrevMat)
summary(model)

PrevMatSitka<-PrevMat[PrevMat$Site=='Sitka Sound',]
PrevMatPWS<-PrevMat[PrevMat$Site=='Prince William Sound',]


library(dplyr)
library(forestplot)

results <- lapply(unique(PrevMatSitka$Cohort), function(pos) {
  fit = filter(PrevMatSitka, Cohort == pos) %>% 
    lm(data = ., Prevalence ~ Age, weight=Weight)
  Age_lm_index = 2 # the second term in lm() output is "stm_white"
  coefficient = coef(fit)[Age_lm_index]
  lb = confint(fit)[Age_lm_index,1] # lower bound confidence 
  ub = confint(fit)[Age_lm_index,2] # upper bound confidence
  output = data.frame(Cohort = pos, coefficient, lb, ub)
  return(output)
}) %>% bind_rows() # bind_rows() combines output from each model in the list

resultsSitka<-results[!is.na(results$lb),]
resultsSitka$Site<-'Sitka Sound'

results <- lapply(unique(PrevMatPWS$Cohort), function(pos) {
  fit = filter(PrevMatPWS, Cohort == pos) %>% 
    lm(data = ., Prevalence ~ Age, weight=Weight)
  Age_lm_index = 2 # the second term in lm() output is "stm_white"
  coefficient = coef(fit)[Age_lm_index]
  lb = confint(fit)[Age_lm_index,1] # lower bound confidence 
  ub = confint(fit)[Age_lm_index,2] # upper bound confidence
  output = data.frame(Cohort = pos, coefficient, lb, ub)
  return(output)
}) %>% bind_rows() # bind_rows() combines output from each model in the list

resultsPWS<-results[!is.na(results$lb),]
resultsPWS$Site<-'Prince William Sound'


results <- lapply(unique(PrevMat$Site), function(pos) {
  fit = filter(PrevMat, Site == pos) %>% 
    lm(data = ., Prevalence ~ Age)
  Age_lm_index = 2 # the second term in lm() output is "stm_white"
  coefficient = coef(fit)[Age_lm_index]
  lb = confint(fit)[Age_lm_index,1] # lower bound confidence 
  ub = confint(fit)[Age_lm_index,2] # upper bound confidence
  output = data.frame(Site = pos, coefficient, lb, ub)
  return(output)
}) %>% bind_rows() # bind_rows() combines output from each model in the list

resultsSite<-results[!is.na(results$lb),]
resultsSite$Cohort<-ifelse(resultsSite$Site=='Sitka Sound', 'Sitka Sound', 'Prince William Sound')




results<-rbind(resultsPWS, resultsSitka)
results<-rbind(results, resultsSite)


results$ub<-((results$ub-results$coefficient)/1.96)+results$coefficient
results$lb<-((results$lb-results$coefficient)/1.96)+results$coefficient


results$ub<-ifelse(results$ub>.2, .2, results$ub)
results$lb<-ifelse(results$lb< -.2, -.2, results$lb)

library('ggplot2') 

results<-results[results$Cohort!='Sitka Sound',]
results<-results[results$Cohort!='Prince William Sound',]


p <- ggplot(results, aes(x=Cohort, y=coefficient, ymin=lb, ymax=ub,col=Site,fill=Site)) + 
  #specify position here
  geom_linerange(size=4,position=position_dodge(width = .8), alpha=.4) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = .8)) +
  scale_x_discrete(name="Predictor") +
  scale_y_continuous(name="Coefficient (effect of age on prevalence)", limits = c(-.2,.2)) +
  coord_flip() +
  theme_minimal()
  
p +facet_grid(Site~.)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ theme(legend.position= "none")

