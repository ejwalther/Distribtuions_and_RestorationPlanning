library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)

MULTITEST_ULO<-function(dataset,cutpoint,species,prediction,filename){
 speciescheck<-species
 summaryforalldata<-data.frame(matrix(NA,nrow =length(cutpoint),ncol = length(speciescheck)))
 colnames(summaryforalldata)<-speciescheck
 rownames(summaryforalldata)<-cutpoint
 for(s in speciescheck){
   storedat<-data.frame()
   sumoflength<-NULL
   cuts<-cutpoint
   for( i in 1:length(cuts)){
     storedat<-LocationofULOs_SPcheck(dataset,cuts[i],s,prediction)
     write.csv(storedat,paste0(filename,cuts[i],"_wUpstreamMn_merged.csv"))
   sumoflength[i]<-(sum(subset(storedat,storedat$Range=="Within")[,"Shape_Leng"]))/1000
   }
   summaryforalldata[,s]<-sumoflength
 }
 return(summaryforalldata)
 }

filenameCo<-"Predicted_CohoULOs_glmm"
filenameSthd<-"Predicted_SthdULOs_glmm"
filenameCH<-"Predicted_ChumULOs_glmm"
cutpointUM<-seq(0.1,0.9,by=0.05)
checks<-c("pCo","p.CoUp","p.CoL")
checksSD<-c("pSD","p.SDUp","p.SDL")
checksCH<-c("pCH","p.CHUp","p.CHL")
Co_CI<-MULTITEST_ULO(orderedNHD_allCo,cutpointUM,checks,"PCo_ULO",filenameCo)
SD_CI<-MULTITEST_ULO(orderedNHD_allSD,cutpointUM,checksSD,"pSD_ULO",filenameSthd)
CH_CI<-MULTITEST_ULO(orderedNHD_allCH,cutpointUM,checksCH,"pCH_ULO",filenameCH)

CH_all<-cbind(CH_CI,cutpointUM)
colnames(CH_all)<-c("Mean","Upper","Lower","Cut")
SD_all<-cbind(SD_CI,cutpointUM)
colnames(SD_all)<-c("Mean","Upper","Lower","Cut")
Co_all<-cbind(Co_CI,cutpointUM)
colnames(Co_all)<-c("Mean","Upper","Lower","Cut")

Co.predReg<-lm(Mean~Cut,data = Co_all)
summary(SD.predReg<-lm(Mean~Cut,data = SD_all))
summary(CH.predReg<-lm(Mean~Cut,data = CH_all))

PerChange<-function(x,y){
((x-y)/abs(x))*100
  }


Co_all<-subset(allsppMulticut_CI,allsppMulticut_CI$Species=="Coho")
PerChange(max(Co_all$Mean),min(Co_all$Mean))
(max(Co_all$Mean)-min(Co_all$Mean))


SD_all<-subset(allsppMulticut_CI,allsppMulticut_CI$Species=="Steelhead")
PerChange(max(SD_all$Mean),min(SD_all$Mean))
(max(SD_all$Mean)-min(SD_all$Mean))

PerChange(max(CH_all$Mean),min(CH_all$Mean))
PerChange(2514.4,830.3)
(max(CH_all$Mean)-min(CH_all$Mean))



allsppMulticut_CI<-rbind(Co_all,SD_all,CH_all)
allsppMulticut_CI$Species<-c(rep("Coho",17),rep("Steelhead",17),rep("Chum",17))
write.csv(allsppMulticut_CI,"PDT_values_Range_allspp.csv")

allspp_CI<-rbind(Co_all[c(4,9,14),],SD_all[c(4,9,14),],CH_all[c(4,9,14),])
allspp_CI$Species<-c(rep("Coho",3),rep("Steelhead",3),rep("Chum",3))
Species<-c(rep("Coho",3),rep("Steelhead",3),rep("Chum",3))
allspp_CI$Cutpoint<-rep(c("0.25","0.50","0.75"),3)
cbind(round(allspp_CI[,1:3],1),Species)

###----------------------------###
###       Figure 2             ###
###----------------------------###

Multispp_nulticutUM<-ggplot()+
  geom_point(data=Co_all,aes(x=Cut,y=Mean),color="salmon", size = 2)+
  geom_errorbar(data=Co_all,aes(x=Cut,ymin = Lower, ymax = Upper),size = .5, width = .05, position=position_dodge(width=0.1),color="salmon") + 
  geom_point(data=SD_all,aes(x=Cut,y=Mean),color="grey53", size = 2)+
  geom_errorbar(data=SD_all,aes(x=Cut,ymin = Lower, ymax = Upper),size = .5, width = .05, position=position_dodge(width=0.1),color="grey53") +
  geom_point(data=CH_all,aes(x=Cut,y=Mean),color="purple", size = 2)+
  geom_errorbar(data=CH_all,aes(x=Cut,ymin = Lower, ymax = Upper),size = .5, width = .05, position=position_dodge(width=0.1),color="purple") +
  #scale_y_continuous(limits=c(693,910),breaks = seq(2500,7500,by=2500))+
  theme_classic()+
  xlab("Probability Decision Threshold") + ylab("Range of Occurrence (km)")+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))


allsppMulticut_CI$Species<-factor(allsppMulticut_CI$Species, levels = c("Coho","Steelhead","Chum"),labels = c("Coho","Steelhead","Chum"))

Multispp_nulticutWithLegend<-ggplot(allsppMulticut_CI,aes(x=Cut,y=Mean,colour=Species))+
  geom_point(size = 2)+
  geom_errorbar(aes(x=Cut,ymin = Lower, ymax = Upper),size = .5, width = .05) +
  scale_color_manual(values = c("salmon","grey53","purple"))+
  theme_classic()+
  xlab("Probability Decision Threshold") + ylab("Range of Occurrence (km)")+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=18))+
  theme(legend.position = c(0.9, 0.8))

ppi <- 600 
tiff("Multispp_nulticut_wLegend.tiff", width=10, height=8, units = 'in', res = 600) #save as fi
Multispp_nulticutWithLegend
dev.off()


###------------------------------------------------------------###
##    Objective 2: Comparing predicted range of occurrence
##                 to currently described distribution  
##-------------------------------------------------------------###

#coho summmary----------------
Coho_SWIFD<-read.csv("Coho_SWIFD_NHDsegs_NoGH.csv")
Coho_SWIFD_desc<-read.csv("Coho_SWIFD_Desc.csv")
Coho_SWIFD_NHD<-unique(Coho_SWIFD$NHD_FID)



Co.5within<-subset(CohoULOs0.5UM,CohoULOs0.5UM$Range=="Within")
Co.25within<-subset(CohoULOs0.25UM,CohoULOs0.25UM$Range=="Within")
Co.75within<-subset(CohoULOs0.75UM,CohoULOs0.75UM$Range=="Within")

Co.5within_NHD<-unique(Co.5within$NHD_FID)
Co.25within_NHD<-unique(Co.25within$NHD_FID)
Co.75within_NHD<-unique(Co.75within$NHD_FID)



Co_SWIFD_length<-sum(Coho_SWIFD[,"Shape_Leng"])/1000

Co.25within_length<-sum(Co.25within[,"Shape_Leng"])/1000
Co.5within_length<-sum(Co.5within[,"Shape_Leng"])/1000
Co.75within_length<-sum(Co.75within[,"Shape_Leng"])/1000

#shared coho
Sharedwithin_COpred0.5<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%in%Co.5within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_COpred0.25<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%in%Co.25within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_COpred0.75<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%in%Co.75within_NHD)[,"Shape_Leng"]))/1000

#only swifd coho
onlySWIFD_co0.5<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.5within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_co0.25<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.25within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_co0.75<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.75within_NHD)[,"Shape_Leng"]))/1000

#swifd sum 0.5
onlySWIFD_co0.5dat<-subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.5within_NHD)
onlySWIFD_co0.5datNHD<-unique(onlySWIFD_co0.5dat$NHD_FID)

DesconlySWIFD_co0.5<-subset(Coho_SWIFD_desc,Coho_SWIFD_desc$NHD_FID%in%onlySWIFD_co0.5datNHD)
DesconlySWIFD_co0.5%>%
  group_by(DISTTYPE_D) %>%
  summarise(totlngth=sum(Shape_Leng)/1000)

#swifd sum 0.75
onlySWIFD_co0.75dat<-subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.75within_NHD)
onlySWIFD_co0.75datNHD<-unique(onlySWIFD_co0.75dat$NHD_FID)

DesconlySWIFD_co0.75<-subset(Coho_SWIFD_desc,Coho_SWIFD_desc$NHD_FID%in%onlySWIFD_co0.75datNHD)
DesconlySWIFD_co0.75%>%
  group_by(DISTTYPE_D) %>%
  summarise(totlngth=sum(Shape_Leng)/1000)


#export summary data
write.csv(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.5within_NHD),"coho_onlySWIFD.csv")
write.csv(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.25within_NHD),"coho_onlySWIFD_25.csv")
write.csv(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.75within_NHD),"coho_onlySWIFD_75.csv")

#only pred coho
onlyPred_co0.5<-(sum(subset(Co.5within,Co.5within$NHD_FID%nin%Coho_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_co0.25<-(sum(subset(Co.25within,Co.25within$NHD_FID%nin%Coho_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_co0.75<-(sum(subset(Co.75within,Co.75within$NHD_FID%nin%Coho_SWIFD_NHD)[,"Shape_Leng"]))/1000

#check values 0.5
(Sharedwithin_COpred0.5+onlyPred_co0.5)
((sum(subset(CohoULOs0.5UM,CohoULOs0.5UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_COpred0.5+onlySWIFD_co0.5)
sum(Coho_SWIFD$Shape_Leng/1000)

#check values 0.25
(Sharedwithin_COpred0.25+onlyPred_co0.25)
((sum(subset(CohoULOs0.25UM,CohoULOs0.25UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_COpred0.25+onlySWIFD_co0.25)
sum(Coho_SWIFD$Shape_Leng/1000)

#check values 0.75
(Sharedwithin_COpred0.75+onlyPred_co0.75)
((sum(subset(CohoULOs0.75UM,CohoULOs0.75UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_COpred0.75+onlySWIFD_co0.75)
sum(Coho_SWIFD$Shape_Leng/1000)



#identify where variance occures
SWIFDonlydatCO<-subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.5within_NHD)
PredOnlydat_co<-subset(Co.5within,Co.5within$NHD_FID%nin%Coho_SWIFD_NHD)
  
HUCSumSWIFDCo<-as.data.frame(SWIFDonlydatCO%>%
  group_by(HUC_10)%>%
  summarise(streamlng=sum(Shape_Leng)/1000))

HUCSumPredCo<-as.data.frame(PredOnlydat_co%>%
  group_by(HUC_10)%>%
  summarise(streamlng=sum(Shape_Leng)/1000))

comparedifCO<-cbind(HUCSumSWIFDCo,HUCSumPredCo[,2])
colnames(comparedifCO)<-c("HUC_10","SWIFD","Prediction")
comparedifCO$HUC_10<-as.factor(comparedifCO$HUC_10)

long_dat<-comparedifCO%>% gather(Type, Dist, 2:3)


ggplot(long_dat,aes(x=HUC_10,y=Dist))+
  geom_col(aes(fill=Type, group=Type),position = "dodge")


#check values
(Sharedwithin_COpred+onlyPred_co)
((sum(subset(CohoULOs0.5UM,CohoULOs0.5UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_COpred+onlySWIFD_co)
sum(Coho_SWIFD$Shape_Leng/1000)

#steelhead summmary--------------
SD_SWIFD<-read.csv("Sthd_SWIFD_NHDsegs_NoGH.csv")
SD_SWIFD_desc<-read.csv("Sthd_SWIFD_Desc.csv")
SD_SWIFD_NHD<-unique(SD_SWIFD$NHD_FID)
SD.5within<-subset(SteelheadULOs0.5UM,SteelheadULOs0.5UM$Range=="Within")
SD.25within<-subset(SteelheadULOs0.25UM,SteelheadULOs0.25UM$Range=="Within")
SD.75within<-subset(SteelheadULOs0.75UM,SteelheadULOs0.75UM$Range=="Within")

SD.5within_NHD<-unique(SD.5within$NHD_FID)
SD.25within_NHD<-unique(SD.25within$NHD_FID)
SD.75within_NHD<-unique(SD.75within$NHD_FID)

SD_SWIFD_length<-sum(SD_SWIFD[,"Shape_Leng"])/1000

#shared steelhead
Sharedwithin_SDpred0.5<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%in%SD.5within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_SDpred0.25<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%in%SD.25within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_SDpred0.75<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%in%SD.75within_NHD)[,"Shape_Leng"]))/1000

#only Swifd steelhead
onlySWIFD_SD0.5<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%nin%SD.5within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_SD0.25<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%nin%SD.25within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_SD0.75<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%nin%SD.75within_NHD)[,"Shape_Leng"]))/1000

#only preds sthd
onlyPred_SD0.5<-(sum(subset(SD.5within,SD.5within$NHD_FID%nin%SD_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_SD0.25<-(sum(subset(SD.25within,SD.25within$NHD_FID%nin%SD_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_SD0.75<-(sum(subset(SD.75within,SD.75within$NHD_FID%nin%SD_SWIFD_NHD)[,"Shape_Leng"]))/1000

#check values 0.5
(Sharedwithin_SDpred0.5+onlyPred_SD0.5)
((sum(subset(SteelheadULOs0.5UM,SteelheadULOs0.5UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_SDpred0.5+onlySWIFD_SD0.5)
sum(SD_SWIFD$Shape_Leng/1000)

#check values 0.25
(Sharedwithin_SDpred0.25+onlyPred_SD0.25)
((sum(subset(SteelheadULOs0.25UM,SteelheadULOs0.25UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_SDpred0.25+onlySWIFD_SD0.25)
sum(SD_SWIFD$Shape_Leng/1000)

#check values 0.75
(Sharedwithin_SDpred0.75+onlyPred_SD0.75)
((sum(subset(SteelheadULOs0.75UM,SteelheadULOs0.75UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_SDpred0.75+onlySWIFD_SD0.75)
sum(SD_SWIFD$Shape_Leng/1000)

#chum summmary--------------
Chum_SWIFD<-read.csv("Chum_SWIFD_NHDsegs_NoGH.csv")
Chum_SWIFD_desc<-read.csv("Chum_SWIFD_Desc.csv")
CH_SWIFD_NHD<-unique(Chum_SWIFD$NHD_FID)
CH.5within<-subset(ChumULOs0.5UM,ChumULOs0.5UM$Range=="Within")
CH.75within<-subset(ChumULOs0.75UM,ChumULOs0.75UM$Range=="Within")
CH.25within<-subset(ChumULOs0.25UM,ChumULOs0.25UM$Range=="Within")

CH.5within_NHD<-unique(CH.5within$NHD_FID)
CH.25within_NHD<-unique(CH.25within$NHD_FID)
CH.75within_NHD<-unique(CH.75within$NHD_FID)

CH_SWIFD_length<-sum(Chum_SWIFD[,"Shape_Leng"])/1000

#Shared amount 
Sharedwithin_CHpred<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%in%CH.5within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_CHpred.25<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%in%CH.25within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_CHpred.75<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%in%CH.75within_NHD)[,"Shape_Leng"]))/1000

Sharedwithin_CHpred.75/CH_SWIFD_length
Sharedwithin_CHpred.25/CH_SWIFD_length

#only swifd
onlySWIFD_CH<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%nin%CH.5within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_CH.25<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%nin%CH.25within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_CH.75<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%nin%CH.75within_NHD)[,"Shape_Leng"]))/1000

#only preds chum
onlyPred_CH<-(sum(subset(CH.5within,CH.5within$NHD_FID%nin%CH_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_CH.25<-(sum(subset(CH.25within,CH.25within$NHD_FID%nin%CH_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_CH.75<-(sum(subset(CH.75within,CH.75within$NHD_FID%nin%CH_SWIFD_NHD)[,"Shape_Leng"]))/1000


#check values chum  0.5
(Sharedwithin_CHpred+onlyPred_CH)
((sum(subset(ChumULOs0.5UM,ChumULOs0.5UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_CHpred+onlySWIFD_CH)
sum(Chum_SWIFD$Shape_Leng/1000)

#check values chum  0.25
(Sharedwithin_CHpred.25+onlyPred_CH.25)
((sum(subset(ChumULOs0.25UM,ChumULOs0.25UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_CHpred.25+onlySWIFD_CH.25)
sum(Chum_SWIFD$Shape_Leng/1000)

#check values chum  0.75
(Sharedwithin_CHpred.25+onlyPred_CH.25)
((sum(subset(ChumULOs0.25UM,ChumULOs0.25UM$Range=="Within")[,"Shape_Leng"]))/1000)

(Sharedwithin_CHpred.25+onlySWIFD_CH.25)
sum(Chum_SWIFD$Shape_Leng/1000)


###----------------------###
###     Figure 3
###----------------------###
allspp_CI_melt<-melt(allspp_CI,id.vars = c("Species","Cutpoint"),measure.vars = c("Upper", "Lower","Mean"))

PCohocompV<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Coho")[,"value"]
CoPercDif<-(PCohocompV-Co_SWIFD_length)/((PCohocompV+Co_SWIFD_length)/2)*100
CoPerChng<-((PCohocompV-Co_SWIFD_length)/(Co_SWIFD_length))*100
CohoPreds<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Coho")
CohoPreds$PerDif<-CoPercDif
CohoPreds$PerChng<-CoPerChng
CohoPreds
CohoPreds$KmDiff<-CohoPreds$value-Co_SWIFD_length


PSTHDcompV<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Steelhead")[,"value"]
SthdPercDif<-(PSTHDcompV-SD_SWIFD_length)/((PSTHDcompV+SD_SWIFD_length)/2)*100
SthdoPerChng<-((PSTHDcompV-SD_SWIFD_length)/(SD_SWIFD_length))*100
SthdPreds<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Steelhead")
SthdPreds$PerDif<-SthdPercDif
SthdPreds$PerChng<-SthdoPerChng
SthdPreds
SthdPreds$KmDiff<-SthdPreds$value-SD_SWIFD_length



PChumcompV<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Chum")[,"value"]
ChumPercDif<-(PChumcompV-CH_SWIFD_length)/((PChumcompV+CH_SWIFD_length)/2)*100
ChumPerChng<-((PChumcompV-CH_SWIFD_length)/(CH_SWIFD_length))*100
ChumPreds<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Chum")
ChumPreds$PerDif<-ChumPercDif
ChumPreds$PerChng<-ChumPerChng
ChumPreds
ChumPreds$kmDiff<-ChumPreds$value-CH_SWIFD_length

PerDif.Table<-rbind(CohoPreds,SthdPreds,ChumPreds)
write.csv(PerDif.Table,"PerDifSWIFD_Table.csv") #export summary table

ChumPreds.mean<-subset(ChumPreds,ChumPreds$variable=="Mean")
ChumPreds.mean<-ChumPreds.mean[,-c(4,5)]
spread.ChumPreds.mean<-ChumPreds.mean %>% spread(variable,PerChng)
ChumPreds.CI<-subset(ChumPreds,ChumPreds$variable!="Mean")
ChumPreds.CI<-ChumPreds.CI[,-c(4,5)]
spread.ChumPreds.CI<-ChumPreds.CI %>% spread(variable,PerChng)
chumPerDif<-cbind(spread.ChumPreds.CI,spread.ChumPreds.mean[,"Mean"])
colnames(chumPerDif)<-c("Species","Cutpoint","Upper","Lower","Mean")

SthdPreds.mean<-subset(SthdPreds,SthdPreds$variable=="Mean")
SthdPreds.mean<-SthdPreds.mean[,-c(4,5)]
spread.SthdPreds.mean<-SthdPreds.mean %>% spread(variable,PerChng)
SthdPreds.CI<-subset(SthdPreds,SthdPreds$variable!="Mean")
SthdPreds.CI<-SthdPreds.CI[,-c(4,5)]
spread.SthdPreds.CI<-SthdPreds.CI %>% spread(variable,PerChng)
sthdPerDif<-cbind(spread.SthdPreds.CI,spread.SthdPreds.mean[,"Mean"])
colnames(sthdPerDif)<-c("Species","Cutpoint","Upper","Lower","Mean")

CohoPreds.mean<-subset(CohoPreds,CohoPreds$variable=="Mean")
CohoPreds.mean<-CohoPreds.mean[,-c(4,5)]
spread.CohoPreds.mean<-CohoPreds.mean %>% spread(variable,PerChng)
CohoPreds.CI<-subset(CohoPreds,CohoPreds$variable!="Mean")
CohoPreds.CI<-CohoPreds.CI[,-c(4,5)]
spread.CohoPreds.CI<-CohoPreds.CI %>% spread(variable,PerChng)
cohoPerDif<-cbind(spread.CohoPreds.CI,spread.CohoPreds.mean[,"Mean"])
colnames(cohoPerDif)<-c("Species","Cutpoint","Upper","Lower","Mean")

PerDifAllspp<-rbind(cohoPerDif,sthdPerDif,chumPerDif)

PerDifAllspp$Species<-factor(PerDifAllspp$Species, levels = c("Coho","Steelhead","Chum"),labels = c("Coho","Steelhead","Chum"))

Multispp_PerDif<-ggplot(PerDifAllspp,aes(x=Mean, y=Cutpoint,color=Species,group=Species))+
  geom_point(aes(shape=Species),size = 4.5,position=position_dodge(width=0.3))+
  geom_errorbar(aes(xmin = Lower, xmax = Upper,y=Cutpoint), size = .5, width = .2, position=position_dodge(width=0.3))+
  theme_classic()+
  geom_vline(xintercept = 0, lty=3)+
  scale_x_continuous(limits=c(-50,200),breaks = seq(-50,200,by=50))+
  xlab("Percent Change") + ylab("Probabiltiy Decision Threshold")+
  scale_color_manual(values=c("salmon","grey53","purple"))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  #theme(legend.position="none")
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.position = c(.9,.8))

ppi <- 600 
tiff("Multispp_PerChng_Fig3.tiff", width=8, height=6, units = 'in', res = 600) #save as fi
Multispp_PerDif
dev.off()

