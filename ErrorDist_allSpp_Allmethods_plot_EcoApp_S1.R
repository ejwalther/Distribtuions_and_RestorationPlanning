library(MuMIn)
library(merTools)
library(lme4)
library(dplyr)
library(ResourceSelection)
library(pROC)
library(lmtest)
library(ggplot2)

setwd(...) #set working directory


#-----------
#   COHO
#-----------
#---------------------------------------------------------
## adjusted error distance function to incorportate neighborhood means for predicting ULO locations
CoEDplotDat_mn<-read.csv("Coho error distance data for all species plot_500_Mn.csv")
CoEDplotDat_mn$Method<-rep("NP",length.out=length(CoEDplotDat_mn[,"ErrorDist"]))
CoEDplotDat_mn<-CoEDplotDat_mn[,-1]
head(CoEDplotDat_mn)
#---------------------------------------------------------
## adjusted error distance function to incorportate upstream means for predicting ULO locations
CoEDplotDat_UP<-read.csv("Coho error distance data for all species plot_500_SPcheck.csv")
CoEDplotDat_UP$Method<-rep("UP",length.out=length(CoEDplotDat_UP[,"ErrorDist"]))
CoEDplotDat_UP<-CoEDplotDat_UP[,-1]
head(CoEDplotDat_UP)
#---------------------------------------------------------
##calculate error distance based on last consequtive reach
CoEDplotDatSingle<-read.csv("Coho error distance data for all species plot_500_sinlge.csv")
CoEDplotDatSingle$Method<-rep("P",length.out=length(CoEDplotDatSingle[,"ErrorDist"]))

head(CoEDplotDat_mn)
head(CoEDplotDat_UP)
head(CoEDplotDatSingle)

#sum stats
NB_EDmean<-mean(CoEDplotDat_mn$ErrorDist)
NB_EDmed<-median(CoEDplotDat_mn$ErrorDist)
NB_EDsd<-sd(CoEDplotDat_mn$ErrorDist)
NB_EDrange<-quantile(CoEDplotDat_mn$ErrorDist)[4]-quantile(CoEDplotDat_mn$ErrorDist)[2]
(NBallsumstat<-cbind(NB_EDmean,NB_EDmed,NB_EDsd,NB_EDrange))

UP_EDmean<-mean(CoEDplotDat_UP$ErrorDist)
UP_EDmed<-median(CoEDplotDat_UP$ErrorDist)
UP_EDsd<-sd(CoEDplotDat_UP$ErrorDist)
UP_EDrange<-quantile(CoEDplotDat_UP$ErrorDist)[4]-quantile(CoEDplotDat_UP$ErrorDist)[2]
(UPallsumstat<-cbind(UP_EDmean,UP_EDmed,UP_EDsd,UP_EDrange))

single_EDmean<-mean(CoEDplotDatSingle$ErrorDist)
single_EDmed<-median(CoEDplotDatSingle$ErrorDist)
single_EDsd<-sd(CoEDplotDatSingle$ErrorDist)
single_EDrange<-quantile(CoEDplotDatSingle$ErrorDist)[4]-quantile(CoEDplotDatSingle$ErrorDist)[2]
(Sinlgeallsumstat<-cbind(single_EDmean,single_EDmed,single_EDsd,single_EDrange))

coEDsumdat<-rbind(Sinlgeallsumstat,NBallsumstat,UPallsumstat)
rownames(coEDsumdat)<-c("P","NB","UP")
colnames(coEDsumdat)<-c("Mean","Median","SD","75perQuantileRange")
coEDsumdat

CoEDallmethods<-rbind(CoEDplotDat_mn,CoEDplotDat_UP,CoEDplotDatSingle)

ggplot(data=CoEDallmethods,aes(x=Species, y=ErrorDist))+
  geom_boxplot(aes(fill=Method))+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-2300,6000))+
  theme_classic()

#-----------
#   Steelhead
#-----------
#---------------------------------------------------------
## adjusted error distance function to incorportate neighborhood means for predicting ULO locations
SDEDplotDat_mn<-read.csv("Steelhead error distance data for all species plot_500_Mn.csv")
SDEDplotDat_mn$Method<-rep("NP",length.out=length(SDEDplotDat_mn[,"ErrorDist"]))
SDEDplotDat_mn<-SDEDplotDat_mn[,-1]
head(SDEDplotDat_mn)
tail(SDEDplotDat_mn)
#---------------------------------------------------------
## adjusted error distance function to incorportate upstream means for predicting ULO locations
SDEDplotDat_UP<-read.csv("Steelhead error distance data for all species plot_500_Upstream.csv")
SDEDplotDat_UP$Method<-rep("UP",length.out=length(SDEDplotDat_UP[,"ErrorDist"]))
SDEDplotDat_UP<-SDEDplotDat_UP[,-1]
head(SDEDplotDat_UP)
tail(SDEDplotDat_UP)
#---------------------------------------------------------
##calculate error distance based on last consequtive reach
SDsingle_EDmean<-read.csv("Steelhead error distance data for all species plot_500_sinlge.csv")
SDsingle_EDmean<-SDsingle_EDmean[,-1]
SDsingle_EDmean$Method<-rep("P",length.out=length(SDsingle_EDmean[,"ErrorDist"]))
head(SDsingle_EDmean)

#sum stats
sthd_NB_EDmean<-mean(SDEDplotDat_mn$ErrorDist)
sthd_NB_EDmed<-median(SDEDplotDat_mn$ErrorDist)
sthd_NB_EDsd<-sd(SDEDplotDat_mn$ErrorDist)
sthd_NB_EDrange<-quantile(SDEDplotDat_mn$ErrorDist)[4]-quantile(SDEDplotDat_mn$ErrorDist)[2]
(sthd_NBallsumstat<-cbind(sthd_NB_EDmean,sthd_NB_EDmed,sthd_NB_EDsd,sthd_NB_EDrange))

sthd_UP_EDmean<-mean(SDEDplotDat_UP$ErrorDist)
sthd_UP_EDmed<-median(SDEDplotDat_UP$ErrorDist)
sthd_UP_EDsd<-sd(SDEDplotDat_UP$ErrorDist)
sthd_UP_EDrange<-quantile(SDEDplotDat_UP$ErrorDist)[4]-quantile(SDEDplotDat_UP$ErrorDist)[2]
(sthd_UPallsumstat<-cbind(sthd_UP_EDmean,sthd_UP_EDmed,sthd_UP_EDsd,sthd_UP_EDrange))

sthd_single_EDmean<-mean(SDsingle_EDmean$ErrorDist)
sthd_single_EDmed<-median(SDsingle_EDmean$ErrorDist)
sthd_single_EDsd<-sd(SDsingle_EDmean$ErrorDist)
sthd_single_EDrange<-quantile(SDsingle_EDmean$ErrorDist)[4]-quantile(SDsingle_EDmean$ErrorDist)[2]
(sthd_Sinlgeallsumstat<-cbind(sthd_single_EDmean,sthd_single_EDmed,sthd_single_EDsd,sthd_single_EDrange))

sthdEDsumdat<-rbind(sthd_Sinlgeallsumstat,sthd_NBallsumstat,sthd_UPallsumstat)
rownames(sthdEDsumdat)<-c("P","NB","UP")
colnames(sthdEDsumdat)<-c("Mean","Median","SD","75perQuantileRange")
sthdEDsumdat


SthdEDallmethods<-rbind(SDEDplotDat_mn,SDEDplotDat_UP,SDsingle_EDmean)

ggplot(data=SthdEDallmethods,aes(x=Species, y=ErrorDist))+
  geom_boxplot(aes(fill=Method))+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-2300,6000))+
  theme_classic()

#-----------
#   Chum
#-----------
#---------------------------------------------------------
## adjusted error distance function to incorportate neighborhood means for predicting ULO locations
CHEDplotDat_mn<-read.csv("Chum error distance data for all species plot_500_Mn.csv")
CHEDplotDat_mn$Method<-rep("NP",length.out=length(CHEDplotDat_mn[,"ErrorDist"]))
CHEDplotDat_mn<-CHEDplotDat_mn[,-1]
head(CHEDplotDat_mn)
#---------------------------------------------------------
## adjusted error distance function to incorportate upstream means for predicting ULO locations
CHEDplotDat_UP<-read.csv("Chum error distance data for all species plot_500_Upstream.csv")
CHEDplotDat_UP$Method<-rep("UP",length.out=length(CHEDplotDat_UP[,"ErrorDist"]))
CHEDplotDat_UP<-CHEDplotDat_UP[,-1]

#---------------------------------------------------------
##calculate error distance based on last consequtive reach
CHEDplotDat_Single<-read.csv("Chum error distance data for all species plot_500_sinlge.csv")
CHEDplotDat_Single$Method<-rep("P",length.out=length(CHEDplotDat_Single[,"ErrorDist"]))
CHEDplotDat_Single<-CHEDplotDat_Single[,-1]
head(CHEDplotDat_Single)

#sum stats
chum_NB_EDmean<-mean(CHEDplotDat_mn$ErrorDist)
chum_NB_EDmed<-median(CHEDplotDat_mn$ErrorDist)
chum_NB_EDsd<-sd(CHEDplotDat_mn$ErrorDist)
chum_NB_EDrange<-quantile(CHEDplotDat_mn$ErrorDist)[4]-quantile(CHEDplotDat_mn$ErrorDist)[2]
(chum_NBallsumstat<-cbind(chum_NB_EDmean,chum_NB_EDmed,chum_NB_EDsd,chum_NB_EDrange))

chum_UP_EDmean<-mean(CHEDplotDat_UP$ErrorDist)
chum_UP_EDmed<-median(CHEDplotDat_UP$ErrorDist)
chum_UP_EDsd<-sd(CHEDplotDat_UP$ErrorDist)
chum_UP_EDrange<-quantile(CHEDplotDat_UP$ErrorDist)[4]-quantile(CHEDplotDat_UP$ErrorDist)[2]
(chum_UPallsumstat<-cbind(chum_UP_EDmean,chum_UP_EDmed,chum_UP_EDsd,chum_UP_EDrange))

chum_single_EDmean<-mean(CHEDplotDat_Single$ErrorDist)
chum_single_EDmed<-median(CHEDplotDat_Single$ErrorDist)
chum_single_EDsd<-sd(CHEDplotDat_Single$ErrorDist)
chum_single_EDrange<-quantile(CHEDplotDat_Single$ErrorDist)[4]-quantile(CHEDplotDat_Single$ErrorDist)[2]
(chum_Sinlgeallsumstat<-cbind(chum_single_EDmean,chum_single_EDmed,chum_single_EDsd,chum_single_EDrange))

chumEDsumdat<-rbind(chum_Sinlgeallsumstat,chum_NBallsumstat,chum_UPallsumstat)
rownames(chumEDsumdat)<-c("P","NP","UP")
colnames(chumEDsumdat)<-c("Mean","Median","SD","75perQuantileRange")
chumEDsumdat


ChumEDallmethods<-rbind(CHEDplotDat_mn,CHEDplotDat_UP,CHEDplotDat_Single)


#combine coho, steelhead, and chum species dataframe and plot
allspp_EDplot<-rbind(CoEDallmethods,SthdEDallmethods,ChumEDallmethods)
allspp_EDplot$Method<-factor(allspp_EDplot$Method,levels=c("P","NP","UP"),labels=c("P","NP","UP"))
allspp_EDplot$Species<-factor(allspp_EDplot$Species,levels=c("Coho","Steelhead","Chum"),labels=c("Coho","Steelhead","Chum"))

ErrorDist_allmethodAllsppplot<-ggplot(data=allspp_EDplot,aes(x=Species, y=ErrorDist))+
  geom_boxplot(aes(fill=Method))+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-10000,15000))+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

ppi <- 600 
png("ErrorDist_allSpp_allmehods_plot_Figure_S1.png", width=8, height=6, units = 'in', res = 600) #save as fi
ErrorDist_allmethodAllsppplot
dev.off()
