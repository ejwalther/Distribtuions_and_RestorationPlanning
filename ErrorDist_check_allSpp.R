rm(list=ls(all=TRUE)) #clear workspace

library(MuMIn)
library(merTools)
library(lme4)
library(dplyr)
library(ResourceSelection)
library(pROC)
library(lmtest)

setwd("E:\\Z167574\\Chehalis upper extent\\04_Anaylsis\\GIS")


zscore<-function(x,y,z){
  (x-y)/z
}


#error didstance among ULO prediciton methods

#-----------
#   COHO
#-----------

#LOAD AND FORMAT DATA
cohodatNBall<-read.csv("CoAllYears_01s_NHD_noBarriers_HUC_R.csv") #no barrier
cohodatNBall<-na.omit(cohodatNBall)#remove stream reaches where there are missing covariate data
row.names(cohodatNBall)<-1:length(cohodatNBall[,1]) #change row names for continous sequence
head(cohodatNBall)
cohodatNBall$Geolgoy<-as.factor(cohodatNBall$Geolgoy)
cohodatNBall$HUC_10<-as.factor(cohodatNBall$HUC_10)
cohodatNBall$SiteID<-as.factor(cohodatNBall$SiteID)
#cohodat$Wetland_Type<-as.factor(cohodat$Wetland_Type) #this covariate is not included in the anaylsis so do not need to transform to categorical variable
cohodatNBall$Wetland_P<-as.factor(cohodatNBall$Wetland_P)
cohodatNBall$Arealog<-log10(cohodatNBall$Area) #log transform drainage area variable

NHD_all_forCoPred<-read.csv("NHD_forCoPreds_noBar_R.csv") #all basin data for approproate standardization

NHD_all_forCoPred<-na.omit(NHD_all_forCoPred)#remove stream reaches where there are missing covariate data
row.names(NHD_all_forCoPred)<-1:length(NHD_all_forCoPred[,1]) #change row names for continous sequence
head(NHD_all_forCoPred)
NHD_all_forCoPred$Geolgoy<-as.factor(NHD_all_forCoPred$Geolgoy)
NHD_all_forCoPred$HUC_10<-as.factor(NHD_all_forCoPred$HUC_10)
NHD_all_forCoPred$LLID<-as.factor(NHD_all_forCoPred$LLID)
NHD_all_forCoPred$Wetland_P<-as.factor(NHD_all_forCoPred$Wetland_P)
NHD_all_forCoPred$Arealog<-log10(NHD_all_forCoPred$Area) #log transform drainage area variable



cohodatNBall2<-cohodatNBall
cohodatNBall2[,c(2:7,15,18)]
cohodatNBall2$ArealogZ<-zscore(cohodatNBall$Arealog,mean(NHD_all_forCoPred$Arealog),sd(NHD_all_forCoPred$Arealog))
cohodatNBall2$MAPZ<-zscore(cohodatNBall$mn_precip,mean(NHD_all_forCoPred$mn_precip),sd(NHD_all_forCoPred$mn_precip))
cohodatNBall2$ElevZ<-zscore(cohodatNBall$Elevation_m,mean(NHD_all_forCoPred$Elevation_m),sd(NHD_all_forCoPred$Elevation_m))
names(cohodatNBall2)


COpred_withMNs<-read.csv("CO_predictions_neighborhood_values.csv")

#read in file that has predicted neighborhood means
for(i in 1:length(cohodatNBall2[,1])){
  if(any(any(COpred_withMNs[,"NHD_FID"]==cohodatNBall2[i,"NHD_FID"]))){
  cohodatNBall2[i,"Mn_pOC"]<-COpred_withMNs[as.numeric(rownames(COpred_withMNs[which(COpred_withMNs[,"NHD_FID"]==cohodatNBall2[i,"NHD_FID"]),])),"RASTERVALU"]}
  else{cohodatNBall2[i,"Mn_pOC"]<-NA}
}

length(cohodatNBall2$NHD_FID)-(length(which(is.na(cohodatNBall2$Mn_pOC)))) #wont be equal because some segments are above barriers that were removed from NHD layer for prediction analysis

cohodatNBall2_mn<-cohodatNBall2[!is.na(cohodatNBall2$Mn_pOC),]

source("E:\\Z167574\\R_data\\Chehalis upper extent\\Stopping Criteria\\functionsForStoppingCriteria.R")

#-------------------------------------------------------------
##    adjusted error distance function to incorportate
##    neighborhood means for predicting ULO locations
##------------------------------------------------------------
CoErrorDistCalc_Mn<-function(dataset,Iter,per,cut){
  cutpointV<-cut
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modCO<-glmer(OBSCO_V_2~ArealogZ+MAPZ+ElevZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],MAPZ=DatForPred[,"MAPZ"],ElevZ=DatForPred[,"ElevZ"],Wetland_P=DatForPred[,"Wetland_P"],Geolgoy=DatForPred[,"Geolgoy"],HUC_10=DatForPred[,"HUC_10"])
    #make predictions
    p <- predictInterval(merMod = Best.modCO, newdata = newdat.k, 
                         level = 0.95, n.sims = 1000,
                         stat = "median", which = "all", 
                         include.resid.var = FALSE)
    p.CV<-as.vector(invlogit(p[,2][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect
    p.CVupper<-as.vector(invlogit(p[,3][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
    p.CVlower<-as.vector(invlogit(p[,4][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSCO_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    Mn_pOC<-as.vector(DatForPred[,"Mn_pOC"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength,Mn_pOC)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      StreamProfileForPred<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(StreamProfileForPred[,"p.CV"])#number of rows (i.e stream reaches) being evaluated 
      if(anyNA(StreamProfileForPred[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (j in 1:nC) compR[j]<- {
        if(any(StreamProfileForPred[1:j,"p.CV"]<cutpointV)=="FALSE"){compR[j]<-1}
        else{
          if((any(StreamProfileForPred[j:nC,"p.CV"]>cutpointV)=="TRUE")&(StreamProfileForPred[j,"Mn_pOC"]>cutpointV)){compR[j]<-1}else{compR[j]<-0
          break}}}# value t
      segCount<-sum(compR)
      sumdat<-summarizefostream(StreamProfileForPred,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}

Coerroreval_500_Mn<-CoErrorDistCalc_Mn(cohodatNBall2_mn,500,.2,.5)
write.csv(CHerroreval_500_Mn,"CoErrorDistCV_500_Mn.csv")
Coerroreval_500_Mn<-read.csv("CoErrorDistCV_500_Mn.csv")

Co.speciesID<-rep("Coho",length.out=length(Coerroreval_500_Mn[,"errorDist"]))
Co.methodID<-rep("Neighborhood",length.out=length(Coerroreval_500_Mn[,"errorDist"]))
CoEDplotDat<-NULL
CoEDplotDat<-Coerroreval_500_Mn$errorDist
CoEDplotDat<-as.data.frame(CoEDplotDat)
CoEDplotDat$Species<-Co.speciesID
CoEDplotDat$Method<-Co.methodID
colnames(CoEDplotDat)<-c("ErrorDist","Species","Method")
write.csv(CoEDplotDat,"Coho error distance data for all species plot_500_Mn.csv")

colnames(CoEDplotDat)<-c("ErrorDist","Species")

ggplot(data=CoEDplotDat,aes(x=ErrorDist))+
  geom_histogram()+
  xlab("Error Distance")+ylab("Number of Streams")+
  coord_cartesian(xlim = c(-2300,6000))

ggplot(data=CoEDplotDat,aes(x=Species, y=ErrorDist))+
  geom_boxplot()+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-2300,6000))

mean(CoEDplotDat$ErrorDist)

median(CoEDplotDat$ErrorDist)
sd(CoEDplotDat$ErrorDist)
#-------------------------------------------------------------
##    adjusted error distance function to incorportate
##    upstream means for predicting ULO locations
##------------------------------------------------------------
CoErrorDistCalc_SPcheck<-function(dataset,Iter,per,cut){
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modCO<-glmer(OBSCO_V_2~ArealogZ+MAPZ+ElevZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],MAPZ=DatForPred[,"MAPZ"],ElevZ=DatForPred[,"ElevZ"],Wetland_P=DatForPred[,"Wetland_P"],Geolgoy=DatForPred[,"Geolgoy"],HUC_10=DatForPred[,"HUC_10"])
    #make predictions
    p <- predictInterval(merMod = Best.modCO, newdata = newdat.k, 
                         level = 0.95, n.sims = 1000,
                         stat = "median", which = "all", 
                         include.resid.var = FALSE)
    p.CV<-as.vector(invlogit(p[,2][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect
    p.CVupper<-as.vector(invlogit(p[,3][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
    p.CVlower<-as.vector(invlogit(p[,4][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSCO_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      singlestreamdat<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(singlestreamdat$p.CV) #number of rows (i.e stream reaches) being evaluated 
      if(anyNA(singlestreamdat[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (j in 1:nC) compR[j]<- {
        if(any(singlestreamdat[1:j,"p.CV"]<cut)=="FALSE"){compR[j]<-1}
        else{
          if((singlestreamdat[j,"order"]+5)>max(singlestreamdat[,"order"])){
            if(singlestreamdat[j,"order"]!=max(singlestreamdat[,"order"])){
              if(mean(singlestreamdat[(j+1):max(singlestreamdat[,"order"]),"p.CV"])>cut){compR[j]<-1}else{compR[j]<-0
              break}}
            else{if(singlestreamdat[j,"p.CV"]>cut){compR[j]<-1}else{compR[j]<-0
            break}}}
          else{if(mean(singlestreamdat[(j+1):(j+5),"p.CV"])>cut){compR[j]<-1}else{compR[j]<-0
          break}}
        }}# value t
      segCount<-sum(compR)
      sumdat<-summarizefostream(singlestreamdat,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}


CHerroreval_500_SPcheck<-CoErrorDistCalc_SPcheck(cohodatNBall2_mn,500,.2,.5)
write.csv(CHerroreval_500_SPcheck,"CoErrorDistCV_500_SPcheck.csv")
CHerroreval_500_SPcheck<-read.csv("CoErrorDistCV_500_SPcheck.csv")

Co.speciesIDSP<-rep("Coho",length.out=length(CHerroreval_500_SPcheck[,"errorDist"]))
Co.methodIDSP<-rep("Upstream",length.out=length(CHerroreval_500_SPcheck[,"errorDist"]))
CoEDplotDatSP<-NULL
CoEDplotDatSP<-CHerroreval_500_SPcheck$errorDist
CoEDplotDatSP<-as.data.frame(CoEDplotDatSP)
CoEDplotDatSP$Species<-Co.speciesIDSP
CoEDplotDatSP$Method<-Co.methodIDSP
colnames(CoEDplotDatSP)<-c("ErrorDist","Species","Method")
write.csv(CoEDplotDatSP,"Coho error distance data for all species plot_500_SPcheck.csv")

colnames(CoEDplotDatSP)<-c("ErrorDist","Species")

ggplot(data=CoEDplotDatSP,aes(x=ErrorDist))+
  geom_histogram()+
  xlab("Error Distance")+ylab("Number of Streams")+
  coord_cartesian(xlim = c(-2300,6000))

ggplot(data=CoEDplotDatSP,aes(x=Species, y=ErrorDist))+
  geom_boxplot()+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-2300,6000))


mean(CoEDplotDatSP$ErrorDist)
median(CoEDplotDatSP$ErrorDist)
sd(CoEDplotDatSP$ErrorDist)

#------
##calculate error distance based on last consequtive reach
#------

CoErrorDistCalc_single<-function(dataset,Iter,per,cut){
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modCO<-glmer(OBSCO_V_2~ArealogZ+MAPZ+ElevZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],MAPZ=DatForPred[,"MAPZ"],ElevZ=DatForPred[,"ElevZ"],Wetland_P=DatForPred[,"Wetland_P"],Geolgoy=DatForPred[,"Geolgoy"],HUC_10=DatForPred[,"HUC_10"])
    #make predictions
    p <- predictInterval(merMod = Best.modCO, newdata = newdat.k, 
                         level = 0.95, n.sims = 1000,
                         stat = "median", which = "all", 
                         include.resid.var = FALSE)
    p.CV<-as.vector(invlogit(p[,2][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect
    p.CVupper<-as.vector(invlogit(p[,3][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
    p.CVlower<-as.vector(invlogit(p[,4][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSCO_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop  will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      singlestreamdat<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(singlestreamdat$p.CV) #number of rows (i.e stream reaches) being evaluated 
      if(anyNA(singlestreamdat[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (i in 1:nC) compR[i]<- {
        if(any(singlestreamdat[1:i,"p.CV"]<cut)=="FALSE"){ # value that here is .5 will be generic vairiable "cutpoint" in function
          compR[i]<-1
        }
        else{
          compR[i]<-0
        } 
      }
      segCount<-sum(compR)
      sumdat<-summarizefostream(singlestreamdat,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}

Coerroreval_500_single<-CoErrorDistCalc_single(cohodatNBall2_mn,500,.2,.5)


Co.speciesID_single<-rep("Coho",length.out=length(Coerroreval_500_single[,"errorDist"]))
Co.methodID_single<-rep("Single",length.out=length(Coerroreval_500_single[,"errorDist"]))
CoEDplotDatSingle<-NULL
CoEDplotDatSingle<-Coerroreval_500_single$errorDist
CoEDplotDatSingle<-as.data.frame(CoEDplotDatSingle)
CoEDplotDatSingle$Species<-Co.speciesID_single
CoEDplotDatSingle$Method<-Co.methodID_single
colnames(CoEDplotDatSingle)<-c("ErrorDist","Species","Method")
write.csv(CoEDplotDatSingle,"Coho error distance data for all species plot_500_sinlge.csv")


mean(CoEDplotDatSingle$ErrorDist)
median(CoEDplotDatSingle$ErrorDist)
sd(CoEDplotDatSingle$ErrorDist)



head(CoEDplotDat)
head(CoEDplotDatSP)
head(CoEDplotDatSingle)

CoEDallmethods<-rbind(CoEDplotDat,CoEDplotDatSP,CoEDplotDatSingle)

ggplot(data=CoEDallmethods,aes(x=Species, y=ErrorDist))+
  geom_boxplot(aes(fill=Method))+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-2300,6000))+
  theme_classic()


##------------
#  Steelhead
##------------

#LOAD AND FORMAT DATA
STHDdatall<-read.csv("SDAllYears_01s_NHD_HUC_R.csv") #all streams dataset
STHDdatall<-na.omit(STHDdatall)#remove stream reaches where there are missing covariate data
row.names(STHDdatall)<-1:length(STHDdatall[,1]) #change row names for continous sequence
head(STHDdatall)
STHDdatall$Geolgoy<-as.factor(STHDdatall$Geolgoy)
STHDdatall$HUC_10<-as.factor(STHDdatall$HUC_10)
STHDdatall$Wetland_P<-as.factor(STHDdatall$Wetland_P)
STHDdatall$Arealog<-log10(STHDdatall$Area) #log transform drainage area variable
STHDdatall$SiteID<-as.factor(STHDdatall$SiteID)
STHDdatall$WORK<-STHDdatall$Dist*STHDdatall$Elevation_m


NHD_all_forSDPred<-read.csv("NHD_forSDPreds_noBar_R.csv")
anyNA(NHD_all_forSDPred)

NHD_all_forSDPred<-na.omit(NHD_all_forSDPred)#remove stream reaches where there are missing covariate data
row.names(NHD_all_forSDPred)<-1:length(NHD_all_forSDPred[,1]) #change row names for continous sequence
head(NHD_all_forSDPred)
NHD_all_forSDPred$Geolgoy<-as.factor(NHD_all_forSDPred$Geolgoy)
NHD_all_forSDPred$HUC_10<-as.factor(NHD_all_forSDPred$HUC_10)
NHD_all_forSDPred$LLID<-as.factor(NHD_all_forSDPred$LLID)
NHD_all_forSDPred$Wetland_P<-as.factor(NHD_all_forSDPred$Wetland_P)
NHD_all_forSDPred$Arealog<-log10(NHD_all_forSDPred$Area) #log transform drainage area variable


STHDdatall2<-STHDdatall
names(STHDdatall2)
STHDdatall2$ArealogZ<-zscore(STHDdatall$Arealog,mean(NHD_all_forSDPred$Arealog),sd(NHD_all_forSDPred$Arealog))
STHDdatall2$ElevZ<-zscore(STHDdatall$Elevation_m,mean(NHD_all_forSDPred$Elevation_m),sd(NHD_all_forSDPred$Elevation_m))
STHDdatall2$SlopeZ<-zscore(STHDdatall$slope,mean(NHD_all_forSDPred$slope),sd(NHD_all_forSDPred$slope))
head(STHDdatall2)


SDpred_withMNs<-read.csv("SD_predictions_neighborhood_values.csv")

#read in file that has predicted neighborhood means
for(i in 1:length(STHDdatall2[,1])){
  if(any(any(SDpred_withMNs[,"NHD_FID"]==STHDdatall2[i,"NHD_FID"]))){
    STHDdatall2[i,"Mn_pOC"]<-SDpred_withMNs[as.numeric(rownames(SDpred_withMNs[which(SDpred_withMNs[,"NHD_FID"]==STHDdatall2[i,"NHD_FID"]),])),"Mn_pOC"]}
  else{STHDdatall2[i,"Mn_pOC"]<-NA}
}

length(STHDdatall2$NHD_FID)-(length(which(is.na(STHDdatall2$Mn_pOC)))) #wont be equal because some segments are above barriers that were removed from NHD layer for prediction analysis

STHDdatall2_mn<-STHDdatall2[!is.na(STHDdatall2$Mn_pOC),]


#-------------------------------------------------------------
##    adjusted error distance function to incorportate
##    neighborhood means for predicting ULO locations
##------------------------------------------------------------
SDErrorDistCalc_Mn<-function(dataset,Iter,per,cut){
  cutpointV<-cut
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modSD<-glmer(OBSSD_V_2~ArealogZ+SlopeZ+ElevZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],SlopeZ=DatForPred[,"SlopeZ"],ElevZ=DatForPred[,"ElevZ"],Wetland_P=DatForPred[,"Wetland_P"],Geolgoy=DatForPred[,"Geolgoy"],HUC_10=DatForPred[,"HUC_10"])
    #make predictions
    p <- predictInterval(merMod = Best.modSD, newdata = newdat.k, 
                         level = 0.95, n.sims = 1000,
                         stat = "median", which = "all", 
                         include.resid.var = FALSE)
    p.CV<-as.vector(invlogit(p[,2][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect
    p.CVupper<-as.vector(invlogit(p[,3][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
    p.CVlower<-as.vector(invlogit(p[,4][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSSD_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    Mn_pOC<-as.vector(DatForPred[,"Mn_pOC"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength,Mn_pOC)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      StreamProfileForPred<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(StreamProfileForPred[,"p.CV"])#number of rows (i.e stream reaches) being evaluated 
      if(anyNA(StreamProfileForPred[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (j in 1:nC) compR[j]<- {
        if(any(StreamProfileForPred[1:j,"p.CV"]<cutpointV)=="FALSE"){compR[j]<-1}
        else{
          if((any(StreamProfileForPred[j:nC,"p.CV"]>cutpointV)=="TRUE")&(StreamProfileForPred[j,"Mn_pOC"]>cutpointV)){compR[j]<-1}else{compR[j]<-0
          break}}}# value t
      segCount<-sum(compR)
      sumdat<-summarizefostream(StreamProfileForPred,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}

SDerroreval_500_Mn<-SDErrorDistCalc_Mn(STHDdatall2_mn,500,.2,.5)
write.csv(SDerroreval_500_Mn,"SthdErrorDistCV_500_Mn.csv")
SDerroreval_500_Mn<-read.csv("SthdErrorDistCV_500_Mn.csv")

SD.speciesID<-rep("Steelhead",length.out=length(SDerroreval_500_Mn[,"errorDist"]))
SD.methodID<-rep("Neighborhood",length.out=length(SDerroreval_500_Mn[,"errorDist"]))
SthdEDplotDat<-NULL
SthdEDplotDat<-SDerroreval_500_Mn$errorDist
SthdEDplotDat<-as.data.frame(SthdEDplotDat)
SthdEDplotDat$Species<-SD.speciesID
SthdEDplotDat$Method<-SD.methodID
colnames(SthdEDplotDat)<-c("ErrorDist","Species","Method")
write.csv(SthdEDplotDat,"Steelhead error distance data for all species plot_500_Mn.csv")


#-------------------------------------------------------------
##    adjusted error distance function to incorportate
##    upstream means for predicting ULO locations
##------------------------------------------------------------
SDErrorDistCalc_SPcheck<-function(dataset,Iter,per,cut){
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modSD<-glmer(OBSSD_V_2~ArealogZ+SlopeZ+ElevZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],SlopeZ=DatForPred[,"SlopeZ"],ElevZ=DatForPred[,"ElevZ"],Wetland_P=DatForPred[,"Wetland_P"],Geolgoy=DatForPred[,"Geolgoy"],HUC_10=DatForPred[,"HUC_10"])
    #make predictions
    p <- predictInterval(merMod = Best.modSD, newdata = newdat.k, 
                         level = 0.95, n.sims = 1000,
                         stat = "median", which = "all", 
                         include.resid.var = FALSE)
    p.CV<-as.vector(invlogit(p[,2][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect
    p.CVupper<-as.vector(invlogit(p[,3][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
    p.CVlower<-as.vector(invlogit(p[,4][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSSD_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      singlestreamdat<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(singlestreamdat$p.CV) #number of rows (i.e stream reaches) being evaluated 
      if(anyNA(singlestreamdat[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (j in 1:nC) compR[j]<- {
        if(any(singlestreamdat[1:j,"p.CV"]<cut)=="FALSE"){compR[j]<-1}
        else{
          if((singlestreamdat[j,"order"]+5)>max(singlestreamdat[,"order"])){
            if(singlestreamdat[j,"order"]!=max(singlestreamdat[,"order"])){
              if(mean(singlestreamdat[(j+1):max(singlestreamdat[,"order"]),"p.CV"])>cut){compR[j]<-1}else{compR[j]<-0
              break}}
            else{if(singlestreamdat[j,"p.CV"]>cut){compR[j]<-1}else{compR[j]<-0
            break}}}
          else{if(mean(singlestreamdat[(j+1):(j+5),"p.CV"])>cut){compR[j]<-1}else{compR[j]<-0
          break}}
        }}# value t
      segCount<-sum(compR)
      sumdat<-summarizefostream(singlestreamdat,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}


SDerroreval_500_UP<-SDErrorDistCalc_SPcheck(STHDdatall2_mn,500,.2,.5)
write.csv(SDerroreval_500_UP,"SthdErrorDistCV_500_Upstream.csv")
SDerroreval_500_UP<-read.csv("SthdErrorDistCV_500_Upstream.csv")

SD.speciesID<-rep("Steelhead",length.out=length(SDerroreval_500_UP[,"errorDist"]))
SD.methodID<-rep("Upstream",length.out=length(SDerroreval_500_UP[,"errorDist"]))
SthdEDplotDatUP<-NULL
SthdEDplotDatUP<-SDerroreval_500_UP$errorDist
SthdEDplotDatUP<-as.data.frame(SthdEDplotDatUP)
SthdEDplotDatUP$Species<-SD.speciesID
SthdEDplotDatUP$Method<-SD.methodID
colnames(SthdEDplotDatUP)<-c("ErrorDist","Species","Method")
write.csv(SthdEDplotDatUP,"Steelhead error distance data for all species plot_500_Upstream.csv")


#------
##calculate error distance based on last consequtive reach
#------

SDErrorDistCalc_single<-function(dataset,Iter,per,cut){
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modSD<-glmer(OBSSD_V_2~ArealogZ+SlopeZ+ElevZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],SlopeZ=DatForPred[,"SlopeZ"],ElevZ=DatForPred[,"ElevZ"],Wetland_P=DatForPred[,"Wetland_P"],Geolgoy=DatForPred[,"Geolgoy"],HUC_10=DatForPred[,"HUC_10"])
    #make predictions
    p <- predictInterval(merMod = Best.modSD, newdata = newdat.k, 
                         level = 0.95, n.sims = 1000,
                         stat = "median", which = "all", 
                         include.resid.var = FALSE)
    p.CV<-as.vector(invlogit(p[,2][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect
    p.CVupper<-as.vector(invlogit(p[,3][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
    p.CVlower<-as.vector(invlogit(p[,4][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSSD_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop  will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      singlestreamdat<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(singlestreamdat$p.CV) #number of rows (i.e stream reaches) being evaluated 
      if(anyNA(singlestreamdat[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (i in 1:nC) compR[i]<- {
        if(any(singlestreamdat[1:i,"p.CV"]<cut)=="FALSE"){ # value that here is .5 will be generic vairiable "cutpoint" in function
          compR[i]<-1
        }
        else{
          compR[i]<-0
        } 
      }
      segCount<-sum(compR)
      sumdat<-summarizefostream(singlestreamdat,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}

SDerroreval_500_single<-SDErrorDistCalc_single(STHDdatall2_mn,500,.2,.5)


SD.speciesID_single<-rep("Steelhead",length.out=length(SDerroreval_500_single[,"errorDist"]))
SD.methodID_single<-rep("Single",length.out=length(SDerroreval_500_single[,"errorDist"]))
SthdEDplotDatSingle<-NULL
SthdEDplotDatSingle<-SDerroreval_500_single$errorDist
SthdEDplotDatSingle<-as.data.frame(SthdEDplotDatSingle)
SthdEDplotDatSingle$Species<-SD.speciesID_single
SthdEDplotDatSingle$Method<-SD.methodID_single
colnames(SthdEDplotDatSingle)<-c("ErrorDist","Species","Method")
write.csv(SthdEDplotDatSingle,"Steelhead error distance data for all species plot_500_sinlge.csv")

SthdEDallmethods<-rbind(SthdEDplotDat,SthdEDplotDatUP,SthdEDplotDatSingle)


ggplot(data=SthdEDallmethods,aes(x=Species, y=ErrorDist))+
  geom_boxplot(aes(fill=Method))+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-2300,6000))+
  theme_classic()




ggplot(data=SthdEDallmethods,aes(x=Species, y=ErrorDist))+
  geom_boxplot(aes(fill=Method))+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-2300,6000))+
  theme_classic()

COSDgroupedEDmethods<-rbind(CoEDallmethods,SthdEDallmethods)


ggplot(data=COSDgroupedEDmethods,aes(x=Species, y=ErrorDist))+
  geom_boxplot(aes(fill=Method))+
  xlab("Species")+ylab("Error Distance (m)")+
  coord_cartesian(ylim = c(-5000,10000))+
  theme_classic()




##------------
#  Chum
##------------

#LOAD AND FORMAT DATA
CHdatall<-read.csv("ChumAllYears_01s_NHD_HUC_R.csv") #all
CHdatall<-na.omit(CHdatall)#remove stream reaches where there are missing covariate data
row.names(CHdatall)<-1:length(CHdatall[,1]) #change row names for continous sequence
head(CHdatall)
CHdatall$Geolgoy<-as.factor(CHdatall$Geolgoy)
CHdatall$HUC_10<-as.factor(CHdatall$HUC_10)
CHdatall$Wetland_P<-as.factor(CHdatall$Wetland_P)
CHdatall$Arealog<-log10(CHdatall$Area) #log transform drainage area variable
CHdatall$SiteID<-as.factor(CHdatall$SiteID)
CHdatall$Dist<-as.numeric(CHdatall$Dist)
CHdatall$WORK<-as.numeric(CHdatall$WORK)

#Load NHD Data with segments upstream of chum barriers removed
NHD_all_forCHPred<-read.csv("NHD_forCHPreds_noBarMerge_R.csv")
anyNA(NHD_all_forCHPred)

NHD_all_forCHPred<-na.omit(NHD_all_forCHPred)#remove stream reaches where there are missing covariate data
row.names(NHD_all_forCHPred)<-1:length(NHD_all_forCHPred[,1]) #change row names for continous sequence
head(NHD_all_forCHPred)
NHD_all_forCHPred$Geolgoy<-as.factor(NHD_all_forCHPred$Geolgoy)
NHD_all_forCHPred$HUC_10<-as.factor(NHD_all_forCHPred$HUC_10)
NHD_all_forCHPred$LLID<-as.factor(NHD_all_forCHPred$LLID)
NHD_all_forCHPred$Wetland_P<-as.factor(NHD_all_forCHPred$Wetland_P)
NHD_all_forCHPred$Arealog<-log10(NHD_all_forCHPred$Area) #log transform drainage area variable

head(NHD_all_forCHPred)
names(NHD_all_forCHPred)
length(NHD_all_forCHPred$NHD_FID)

NHD_all_forCHPred2<-NHD_all_forCHPred
names(NHD_all_forCHPred2)
NHD_all_forCHPred2[,c(2:6,8,15)]<-scale(NHD_all_forCHPred2[,c(2:6,8,15)]) #standarize continous variables
head(NHD_all_forCHPred2)


colnames(NHD_all_forCHPred2)[c(8,15)]<-c("ElevZ","ArealogZ")

CHdatall2<-CHdatall
names(CHdatall2)
CHdatall2$ArealogZ<-zscore(CHdatall$Arealog,mean(NHD_all_forCHPred$Arealog),sd(NHD_all_forCHPred$Arealog))
CHdatall2$ElevZ<-zscore(CHdatall$Elevation_m,mean(NHD_all_forCHPred$Elevation_m),sd(NHD_all_forCHPred$Elevation_m))
head(CHdatall2)
#-------------------------------------------------------------
##    adjusted error distance function to incorportate
##    neighborhood means for predicting ULO locations
##------------------------------------------------------------

CHpred_withMNs<-read.csv("CH_predictions_neighborhood_values.csv")

#read in file that has predicted neighborhood means
for(i in 1:length(CHdatall2[,1])){
  if(any(any(CHpred_withMNs[,"NHD_FID"]==CHdatall2[i,"NHD_FID"]))){
    CHdatall2[i,"Mn_pOC"]<-CHpred_withMNs[as.numeric(rownames(CHpred_withMNs[which(CHpred_withMNs[,"NHD_FID"]==CHdatall2[i,"NHD_FID"]),])),"Mn_pOC"]}
  else{CHdatall2[i,"Mn_pOC"]<-NA}
}

length(CHdatall2$NHD_FID)-(length(which(is.na(CHdatall2$Mn_pOC)))) #wont be equal because some segments are above barriers that were removed from NHD layer for prediction analysis

CHdatall2_mn<-CHdatall2[!is.na(CHdatall2$Mn_pOC),]

#-------------------------------------------------------------
##    adjusted error distance function to incorportate
##    neighborhood means for predicting ULO locations
##------------------------------------------------------------
CHErrorDistCalc_Mn<-function(dataset,Iter,per,cut){
  cutpointV<-cut
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modCH<-glm(OBSCH_V_2~ArealogZ+ElevZ+Geolgoy,family = binomial(link=logit),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],ElevZ=DatForPred[,"ElevZ"],Geolgoy=DatForPred[,"Geolgoy"])
    #make predictions
    predict.CH <- predict(Best.modCH, newdat.k, se.fit = T)
    predict.CH <- data.frame(predict.CH$fit, predict.CH$se.fit)
    colnames(predict.CH) <- c("fit", "se.fit")
    length(predict.CH$fit)
    p.CV<-exp(predict.CH$fit)/(1+exp(predict.CH$fit))
    p.CVupper<-exp(predict.CH$fit+1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit+1.96*predict.CH$se.fit))
    p.CVlower<-exp(predict.CH$fit-1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit-1.96*predict.CH$se.fit))#reate a vector of predicted values that include random effect lower CI
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSCH_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    Mn_pOC<-as.vector(DatForPred[,"Mn_pOC"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength,Mn_pOC)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      StreamProfileForPred<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(StreamProfileForPred[,"p.CV"])#number of rows (i.e stream reaches) being evaluated 
      if(anyNA(StreamProfileForPred[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (j in 1:nC) compR[j]<- {
        if(any(StreamProfileForPred[1:j,"p.CV"]<cutpointV)=="FALSE"){compR[j]<-1}
        else{
          if((any(StreamProfileForPred[j:nC,"p.CV"]>cutpointV)=="TRUE")&(StreamProfileForPred[j,"Mn_pOC"]>cutpointV)){compR[j]<-1}else{compR[j]<-0
          break}}}# value t
      segCount<-sum(compR)
      sumdat<-summarizefostream(StreamProfileForPred,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}

CHerroreval_500_Mn<-CHErrorDistCalc_Mn(CHdatall2_mn,500,.2,.5)
write.csv(CHerroreval_500_Mn,"ChumErrorDistCV_500_Mn.csv")
CHerroreval_500_Mn<-read.csv("ChumErrorDistCV_500_Mn.csv")

CH.speciesID<-rep("Chum",length.out=length(CHerroreval_500_Mn[,"errorDist"]))
CH.methodID<-rep("NB",length.out=length(CHerroreval_500_Mn[,"errorDist"]))
ChumEDplotDat<-NULL
ChumEDplotDat<-CHerroreval_500_Mn$errorDist
ChumEDplotDat<-as.data.frame(ChumEDplotDat)
ChumEDplotDat$Species<-CH.speciesID
ChumEDplotDat$Method<-CH.methodID
colnames(ChumEDplotDat)<-c("ErrorDist","Species","Method")
write.csv(ChumEDplotDat,"Chum error distance data for all species plot_500_Mn.csv")


#-------------------------------------------------------------
##    adjusted error distance function to incorportate
##    upstream means for predicting ULO locations
##------------------------------------------------------------
CHErrorDistCalc_SPcheck<-function(dataset,Iter,per,cut){
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modCH<-glm(OBSCH_V_2~ArealogZ+ElevZ+Geolgoy,family = binomial(link=logit),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],ElevZ=DatForPred[,"ElevZ"],Geolgoy=DatForPred[,"Geolgoy"])
    #make predictions
    predict.CH <- predict(Best.modCH, newdat.k, se.fit = T)
    predict.CH <- data.frame(predict.CH$fit, predict.CH$se.fit)
    colnames(predict.CH) <- c("fit", "se.fit")
    length(predict.CH$fit)
    p.CV<-exp(predict.CH$fit)/(1+exp(predict.CH$fit))
    p.CVupper<-exp(predict.CH$fit+1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit+1.96*predict.CH$se.fit))
    p.CVlower<-exp(predict.CH$fit-1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit-1.96*predict.CH$se.fit))#reat
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSCH_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      singlestreamdat<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(singlestreamdat$p.CV) #number of rows (i.e stream reaches) being evaluated 
      if(anyNA(singlestreamdat[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (j in 1:nC) compR[j]<- {
        if(any(singlestreamdat[1:j,"p.CV"]<cut)=="FALSE"){compR[j]<-1}
        else{
          if((singlestreamdat[j,"order"]+5)>max(singlestreamdat[,"order"])){
            if(singlestreamdat[j,"order"]!=max(singlestreamdat[,"order"])){
              if(mean(singlestreamdat[(j+1):max(singlestreamdat[,"order"]),"p.CV"])>cut){compR[j]<-1}else{compR[j]<-0
              break}}
            else{if(singlestreamdat[j,"p.CV"]>cut){compR[j]<-1}else{compR[j]<-0
            break}}}
          else{if(mean(singlestreamdat[(j+1):(j+5),"p.CV"])>cut){compR[j]<-1}else{compR[j]<-0
          break}}
        }}# value t
      segCount<-sum(compR)
      sumdat<-summarizefostream(singlestreamdat,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}


CHerroreval_500_UP<-CHErrorDistCalc_SPcheck(CHdatall2_mn,500,.2,.5)
write.csv(CHerroreval_500_UP,"ChumErrorDistCV_500_Upstream.csv")
CHerroreval_500_UP<-read.csv("ChumErrorDistCV_500_Upstream.csv")

CH.speciesID<-rep("Chum",length.out=length(CHerroreval_500_UP[,"errorDist"]))
CH.methodID<-rep("Upstream",length.out=length(CHerroreval_500_UP[,"errorDist"]))
ChumEDplotDatUP<-NULL
ChumEDplotDatUP<-CHerroreval_500_UP$errorDist
ChumEDplotDatUP<-as.data.frame(ChumEDplotDatUP)
ChumEDplotDatUP$Species<-CH.speciesID
ChumEDplotDatUP$Method<-CH.methodID
colnames(ChumEDplotDatUP)<-c("ErrorDist","Species","Method")
write.csv(ChumEDplotDatUP,"Chum error distance data for all species plot_500_Upstream.csv")


#------
##calculate error distance based on last consequtive reach
#------

CHErrorDistCalc_single<-function(dataset,Iter,per,cut){
  streamdat<-dataset
  ListofStreams<-samplemystreams(streamdat,Iter,per)
  nList<-length(ListofStreams)
  print(nList)
  Alllistsum<-data.frame()
  for(i in 1:nList){
    DatForModFit<-subset(dataset,(!(dataset[,"SiteID"]%in%ListofStreams[[i]])))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    DatForPred<-subset(dataset,(dataset[,"SiteID"]%in%ListofStreams[[i]]))#DATA FIT MODEL EXCLUDING STREAMS USED FOR ERROR DISTATNCE cv
    #fit model
    Best.modCH<-glm(OBSCH_V_2~ArealogZ+ElevZ+Geolgoy,family = binomial(link=logit),data = DatForModFit)
    #set up dataframe for making predictions back on NHD segments
    newdat.k<-data.frame(ArealogZ=DatForPred[,"ArealogZ"],ElevZ=DatForPred[,"ElevZ"],Geolgoy=DatForPred[,"Geolgoy"])
    #make predictions
    predict.CH <- predict(Best.modCH, newdat.k, se.fit = T)
    predict.CH <- data.frame(predict.CH$fit, predict.CH$se.fit)
    colnames(predict.CH) <- c("fit", "se.fit")
    length(predict.CH$fit)
    p.CV<-exp(predict.CH$fit)/(1+exp(predict.CH$fit))
    p.CVupper<-exp(predict.CH$fit+1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit+1.96*predict.CH$se.fit))
    p.CVlower<-exp(predict.CH$fit-1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit-1.96*predict.CH$se.fit))#reat
    #extract NHD segments, observed values, and SiteIDs
    NHD_seg<-as.vector(DatForPred[,"NHD_FID"])#NHD segment numbers
    t_val<-as.vector(DatForPred[,"OBSCH_V_2"])#create a vector of observed values
    SiteID<-as.vector(DatForPred[,"SiteID"])# SiteIDs
    R_Meas<-as.vector(DatForPred[,"R_meas"]) #dist for sorting
    Arealog<-as.vector(DatForPred[,"ArealogZ"])
    Elevation_m<-as.vector((DatForPred[,"Elevation_m"]))
    SubBasin<-as.vector(DatForPred[,"HUC_10"])
    reachlength<-as.vector(DatForPred[,"Shape_Leng"])
    pred.all<-cbind.data.frame(NHD_seg,t_val,p.CV,p.CVlower,p.CVupper,R_Meas,Arealog,Elevation_m,reachlength)
    pred.all$SiteID<-SiteID
    pred.all$SubBasin<-SubBasin
    n<-length(pred.all[,1])
    Orderedstreams<-OrdermyStreamData(pred.all)#function to add order to stream
    SumErrorDist<-data.frame()
    WithinlistStreams<-as.vector(ListofStreams[[i]])
    listcount<- (1:nList)[i] #thsi will be which list based on loop  will be which list based on loop through list
    for(i in 1:length(WithinlistStreams)){
      singlestreamdat<-subset(Orderedstreams,(Orderedstreams[,"SiteID"]%in%WithinlistStreams[i])) #think were this dataframe with be coming from
      compR<-NULL
      nC<-length(singlestreamdat$p.CV) #number of rows (i.e stream reaches) being evaluated 
      if(anyNA(singlestreamdat[,"p.CV"])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
      for (i in 1:nC) compR[i]<- {
        if(any(singlestreamdat[1:i,"p.CV"]<cut)=="FALSE"){ # value that here is .5 will be generic vairiable "cutpoint" in function
          compR[i]<-1
        }
        else{
          compR[i]<-0
        } 
      }
      segCount<-sum(compR)
      sumdat<-summarizefostream(singlestreamdat,segCount)
      listnumber<-rep(listcount,length.out=length(sumdat[,1]))
      Cutpoint<-rep(cut,length.out=length(sumdat[,1]))
      sumdat<-cbind(sumdat,listnumber,Cutpoint)
      SumErrorDist<-rbind(SumErrorDist,sumdat)
    }
    Alllistsum<-rbind(Alllistsum,SumErrorDist)
  }
  return(Alllistsum)
}

CHerroreval_500_single<-CHErrorDistCalc_single(CHdatall2_mn,500,.2,.5)
write.csv(CHerroreval_500_single,"ChumErrorDistCV_500_single.csv")


CH.speciesID_single<-rep("Chum",length.out=length(CHerroreval_500_single[,"errorDist"]))
CH.methodID_single<-rep("Single",length.out=length(CHerroreval_500_single[,"errorDist"]))
ChumEDplotDatSingle<-NULL
ChumEDplotDatSingle<-CHerroreval_500_single$errorDist
ChumEDplotDatSingle<-as.data.frame(ChumEDplotDatSingle)
ChumEDplotDatSingle$Species<-CH.speciesID_single
ChumEDplotDatSingle$Method<-CH.methodID_single
colnames(ChumEDplotDatSingle)<-c("ErrorDist","Species","Method")
write.csv(ChumEDplotDatSingle,"Chum error distance data for all species plot_500_sinlge.csv")

