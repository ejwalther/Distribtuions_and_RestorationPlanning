###------------------------------------------------####
#   Objective 3: Code for generating summary data.
#                Datasets used in this code were
#                exported from ArcMap.    

#               Author: Eric J. Walther
#
###------------------------------------------------####

setwd(...) #set working directory

library(dplyr)
library(tidyr)
library(Hmisc)

AllBarries<-read.csv("Dist_AllCulvert_allspp_FINAL.csv")
Lowerculvertsummary<-read.csv("Dist_Lowestculvert_all.csv")
names(Lowerculvertsummary)

lowerculvertsforanalysis<-subset(AllBarries,AllBarries$SiteId%in%unique(Lowerculvertsummary$SiteId))

Lowerculvertsummary[duplicated(Lowerculvertsummary$SiteId),]
Lowerculvertsummary[Lowerculvertsummary$SiteId%nin%unique(lowerculvertsforanalysis$SiteId)==TRUE,]


CohoHab0.5_blocked<-sum(lowerculvertsforanalysis[,"Lng_50Co"])/1000
CohoHab0.25_blocked<-sum(lowerculvertsforanalysis[,"Lng_25Co"])/1000
CohoHab0.75_blocked<-sum(lowerculvertsforanalysis[,"Lng_75Co"])/1000
SDHab0.5_blocked<-sum(lowerculvertsforanalysis[,"Lng_50SD"])/1000
SDHab0.25_blocked<-sum(lowerculvertsforanalysis[,"Lng_25SD"])/1000
SDHab.075_blocked<-sum(lowerculvertsforanalysis[,"Lng_75SD"])/1000
ChumHab0.5_blocked<-sum(lowerculvertsforanalysis[,"Lng_50CH"])/1000
ChumHab0.25_blocked<-sum(lowerculvertsforanalysis[,"Lng_25CH"])/1000
ChumHab.75_blocked<-sum(lowerculvertsforanalysis[,"Lng_75CH"])/1000

cohoDistB<-rbind(CohoHab0.5_blocked,CohoHab0.25_blocked,CohoHab0.75_blocked)[,1]
SteetheadDistB<-rbind(SDHab0.5_blocked,SDHab0.25_blocked,SDHab.075_blocked)[,1]
ChumDistB<-rbind(ChumHab0.5_blocked,ChumHab0.25_blocked,ChumHab.75_blocked)[,1]

SumDistBlocked<-round(rbind(cohoDistB,SteetheadDistB,ChumDistB),1)
rownames(SumDistBlocked)<-c("Coho","Steelhead","Chum")
SumDistBlocked2<-cbind(SumDistBlocked[,2],SumDistBlocked[,1],SumDistBlocked[,3])
colnames(SumDistBlocked2)<-c("0.25","0.5","0.75")

write.csv(SumDistBlocked2,"SumDistBlocked.csv")#export summary table for distance blocked

##-------------------------------
##  All culvert summary
##-------------------------------

AllBarries[duplicated(AllBarries$SiteId),"SiteId"]

#--------------------------------------------#
#  All barriers summary using a 0.50 PDT
#--------------------------------------------#
AllBarries_Co_No0<-subset(AllBarries,AllBarries$Lng_50Co!=0)
Co_allculsum<-AllBarries_Co_No0%>%
  summarise(Med.dist=median(Lng_50Co),
            Mn.dist = mean(Lng_50Co),
            min.dist = min(Lng_50Co),
            max.dist = max(Lng_50Co),
            number =length(Lng_50Co))

AllBarries_SD_No0<-subset(AllBarries,AllBarries$Lng_50SD!=0)
SD_allculsum<-AllBarries_SD_No0%>%
  summarise(Med.dist=median(Lng_50SD),
            Mn.dist = mean(Lng_50SD),
            min.dist = min(Lng_50SD),
            max.dist = max(Lng_50SD),
            number =length(Lng_50SD))

AllBarries_CH_No0<-subset(AllBarries,AllBarries$Lng_50CH!=0)
CH_allculsum<-as.data.frame(AllBarries_CH_No0%>%
  summarise(Med.dist=median(Lng_50CH),
            Mn.dist = mean(Lng_50CH),
            min.dist = min(Lng_50CH),
            max.dist = max(Lng_50CH),
            number =length(Lng_50CH)))

sumall50<-round(rbind(Co_allculsum,SD_allculsum,CH_allculsum),1)
sumall50[,1:4]<-sumall50[,1:4]/1000
sumall50$Species<-c("Coho","Steelhead","Chum")
sumall50$PerCul<-c((sumall50[1,"number"]/1737)*100,(sumall50[2,"number"]/1737)*100,(sumall50[3,"number"]/1737)*100)
sumall50$PDT<-rep("50",3)
sumall50

#--------------------------------------------#
#  All barriers summary using a 0.25 PDT
#--------------------------------------------#
AllBarries_Co25_No0<-subset(AllBarries,AllBarries$Lng_25Co!=0)
Co25_allculsum<-AllBarries_Co25_No0%>%
  summarise(Med.dist=median(Lng_25Co),
            Mn.dist = mean(Lng_25Co),
            min.dist = min(Lng_25Co),
            max.dist = max(Lng_25Co),
            number =length(Lng_25Co))

AllBarries_SD25_No0<-subset(AllBarries,AllBarries$Lng_25SD!=0)
SD25_allculsum<-AllBarries_SD25_No0%>%
  summarise(Med.dist=median(Lng_25SD),
            Mn.dist = mean(Lng_25SD),
            min.dist = min(Lng_25SD),
            max.dist = max(Lng_25SD),
            number =length(Lng_25SD))

AllBarries_CH25_No0<-subset(AllBarries,AllBarries$Lng_25CH!=0)
CH25_allculsum<-AllBarries_CH25_No0%>%
  summarise(Med.dist=median(Lng_25CH),
            Mn.dist = mean(Lng_25CH),
            min.dist = min(Lng_25CH),
            max.dist = max(Lng_25CH),
            number =length(Lng_25CH))
# 

sumall25<-round(rbind(Co25_allculsum,SD25_allculsum,CH25_allculsum),1)
sumall25[,1:4]<-sumall25[,1:4]/1000
sumall25$Species<-c("Coho","Steelhead","Chum")
sumall25$PerCul<-c((sumall25[1,"number"]/1737)*100,(sumall25[2,"number"]/1737)*100,(sumall25[3,"number"]/1737)*100)
sumall25$PDT<-rep("25",3)
sumall25


#--------------------------------------------#
#  All barriers summary using a 0.75 PDT
#--------------------------------------------#
AllBarries_Co75_No0<-subset(AllBarries,AllBarries$Lng_75Co!=0)
Co75_allculsum<-AllBarries_Co75_No0%>%
  summarise(Med.dist=median(Lng_75Co),
            Mn.dist = mean(Lng_75Co),
            min.dist = min(Lng_75Co),
            max.dist = max(Lng_75Co),
            number =length(Lng_75Co))

AllBarries_SD75_No0<-subset(AllBarries,AllBarries$Lng_75SD!=0)
SD75_allculsum<-AllBarries_SD75_No0%>%
  summarise(Med.dist=median(Lng_75SD),
            Mn.dist = mean(Lng_75SD),
            min.dist = min(Lng_75SD),
            max.dist = max(Lng_75SD),
            number =length(Lng_75SD))

AllBarries_CH75_No0<-subset(AllBarries,AllBarries$Lng_75CH!=0)
CH75_allculsum<-AllBarries_CH75_No0%>%
  summarise(Med.dist=median(Lng_75CH),
            Mn.dist = mean(Lng_75CH),
            min.dist = min(Lng_75CH),
            max.dist = max(Lng_75CH),
            number =length(Lng_75CH))
# 


sumall75<-round(rbind(Co75_allculsum,SD75_allculsum,CH75_allculsum),1)
sumall75[,1:4]<-sumall75[,1:4]/1000
sumall75$Species<-c("Coho","Steelhead","Chum")
sumall75$PerCul<-c((sumall75[1,"number"]/1737)*100,(sumall75[2,"number"]/1737)*100,(sumall75[3,"number"]/1737)*100)
sumall75$PDT<-rep("75",3)
sumall75



allculsumdat<-rbind(sumall25,sumall50,sumall75)



SppCulComp<-read.csv("culvert_sppComp.csv")
head(SppCulComp)


shared<-subset(SppCulComp,SppCulComp$Ovr==3)
blockcul<-shared[,18:20]
minblocked<-NULL
for(i in 1:length(blockcul$LngCo)){
minblocked[i]<-min(blockcul[i,])
}
sum(minblocked)/1000

#sum of total blocked
sum50<-Lowerculvertsummary[,c(2,5,9:11)]
for(i in 1:length(sum50$POINTID)){
sum50$maxlength[i]<-max(sum50[i,c(3:5)])
}
(sum(sum50$maxlength)/1000)

sum25<-Lowerculvertsummary[,c(2,5,13,15,17)]
for(i in 1:length(sum25$POINTID)){
  sum25$maxlength[i]<-max(sum25[i,c(3:5)])
}
sum(sum25$maxlength)/1000

sum75<-Lowerculvertsummary[,c(2,5,12,14,16)]
for(i in 1:length(sum75$POINTID)){
  sum75$maxlength[i]<-max(sum75[i,c(3:5)])
}
(sum(sum75$maxlength)/1000)

rangehabblockedallspp<-round(c((sum(sum75$maxlength)/1000),(sum(sum50$maxlength)/1000),(sum(sum25$maxlength)/1000)),1)


