###------------------------------------------------####
#   Objective 3: Code for generating summary data.
#                Datasets used in this code were
#                exported from ESRI ArcMap GIS.    

#               Author: Eric J. Walther
#
###------------------------------------------------####

setwd(...) #set working directory

library(dplyr)
library(tidyr)
library(Hmisc)

AllBarries<-read.csv("Dist_AllCulvert_allspp_FINAL.csv")
Lowerculvertsummary<-read.csv("Dist_Lowestculvert_all.csv") #these culverts that are lowest down in the river network were identified using a geoprocessing model developed in ArcMap

lowerculvertsforanalysis<-subset(AllBarries,AllBarries$SiteId%in%unique(Lowerculvertsummary$SiteId))

#lowest culvert summary (total amount of habitat blocked without replication)
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

sumall75<-round(rbind(Co75_allculsum,SD75_allculsum,CH75_allculsum),1)
sumall75[,1:4]<-sumall75[,1:4]/1000
sumall75$Species<-c("Coho","Steelhead","Chum")
sumall75$PerCul<-c((sumall75[1,"number"]/1737)*100,(sumall75[2,"number"]/1737)*100,(sumall75[3,"number"]/1737)*100)
sumall75$PDT<-rep("75",3)
sumall75



allculsumdat<-rbind(sumall25,sumall50,sumall75)

#Identify number of culverts that impact more that one species-------------------

#function developed to sumarize the number of unique culverts for each species
summarizeUniqueCulvert<-function(testculvert){
  for (i in 1:ncol(testculvert)){
    colID<-i
    for(r in 1:nrow(testculvert)){
      if((testculvert[r,colID]>0)=="TRUE"){
        testculvert[r,colID]<-1}
      else{testculvert[r,colID]<-0}}
  }
  sumrow<-function(x, output){
    A<-x[1]
    B<-x[2]
    C<-x[3]
    # return sum
    return(A+B+C)
  }
  testculvert$over<-(apply(testculvert,1,sumrow))
  unique<-subset(testculvert,testculvert$over==1)
  countUniqueCulvert<-NULL
  for (i in 1:ncol(unique[,-1])){
    numCul<-sum(unique[,i])
    countUniqueCulvert<-c(countUniqueCulvert,numCul)
  }
  return(countUniqueCulvert)
}

#Summarize Unique Barriers using a 0.50 PDT
AllBarriers_50<-AllBarries[,33:35]
UniqueBar_50<-summarizeUniqueCulvert(AllBarriers_50)

#Summarize Unique Barriers using a 0.75 PDT
AllBarriers_75<-AllBarries[,c(28,30,32)]
UniqueBar_75<-summarizeUniqueCulvert(AllBarriers_75)

#Summarize Unique Barriers using a 0.25 PDT
AllBarriers_25<-AllBarries[,c(27,29,31)]
UniqueBar_25<-summarizeUniqueCulvert(AllBarriers_25)

#summary Table
FinalSumTable<-rbind(cbind(rbind(UniqueBar_25[1],UniqueBar_50[1],UniqueBar_75[1]),rep("Coho",3),c("0.25","0.50,","0.75")),
                     cbind(rbind(UniqueBar_25[3],UniqueBar_50[3],UniqueBar_75[3]),rep("Steelhead",3),c("0.25","0.50,","0.75")),
                     cbind(rbind(UniqueBar_25[2],UniqueBar_50[2],UniqueBar_75[2]),rep("Chum",3),c("0.25","0.50,","0.75")))
colnames(FinalSumTable)<-c("count","Species","PDT")
FinalSumTable
