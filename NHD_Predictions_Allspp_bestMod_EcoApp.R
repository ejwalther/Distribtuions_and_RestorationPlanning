#================================================================================
#        Predicting stream reaches within range of occurrence for coho, steelhead, and chum
#   
#     Author: Eric J Walther
#     Last date modified: April 15, 2021
#=================================================================================

#This script will take all the NHD stream reaches in the Chehalis and predict the probablity that each reach is
#within the range of occurrence for the species of interest. First, each dataset for fitting the best glmmm for each species needs to be imported.
#These models were identified during the analysis for CH1. Datasets are found in respect folders for each model and imported below. 

rm(list=ls(all=TRUE)) #clear workspace

library(ggplot2)
library(MuMIn)
library(merTools)
library(lme4)
library(dplyr)
library(ResourceSelection)
library(pROC)
library(lmtest)
setwd(...) #add source to dataset

zscore<-function(x,y,z){
  (x-y)/z
}

#-------------------------------------------
#LOAD AND FORMAT DATA for coho best model
#-------------------------------------------
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


#-------------------------------------------
#LOAD AND FORMAT DATA for steelhead best model
#-------------------------------------------
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


#-------------------------------------------
#LOAD AND FORMAT DATA for chum best model
#-------------------------------------------
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




#--------------------------------
#       Load NHD Data
#--------------------------------
NHD_all<-read.csv("NHD_wAttsStreamID_08122020_R.csv")
anyNA(NHD_all)
NHD_all<-na.omit(NHD_all)#remove stream reaches where there are missing covariate data
row.names(NHD_all)<-1:length(NHD_all[,1]) #change row names for continous sequence
head(NHD_all)
NHD_all$Geolgoy<-as.factor(NHD_all$Geolgoy)
NHD_all$HUC_10<-as.factor(NHD_all$HUC_10)
NHD_all$LLID<-as.factor(NHD_all$LLID)
NHD_all$Wetland_P<-as.factor(NHD_all$Wetland_P)
NHD_all$Arealog<-log10(NHD_all$Area) #log transform drainage area variable

head(NHD_all)
names(NHD_all)
length(NHD_all$NHD_FID)

NHD_all2<-NHD_all
names(NHD_all2)
NHD_all2[,c(2:6,8,15)]<-scale(NHD_all2[,c(2:6,8,15)]) #standarize continous variables
head(NHD_all2)


#--------------------------------
#       Load NHD Data with segments upstream of coho barriers removed
#--------------------------------
NHD_all_forCoPred<-read.csv("NHD_forCoPreds_noBar_R.csv")
anyNA(NHD_all_forCoPred)

NHD_all_forCoPred<-na.omit(NHD_all_forCoPred)#remove stream reaches where there are missing covariate data
row.names(NHD_all_forCoPred)<-1:length(NHD_all_forCoPred[,1]) #change row names for continous sequence
head(NHD_all_forCoPred)
NHD_all_forCoPred$Geolgoy<-as.factor(NHD_all_forCoPred$Geolgoy)
NHD_all_forCoPred$HUC_10<-as.factor(NHD_all_forCoPred$HUC_10)
NHD_all_forCoPred$LLID<-as.factor(NHD_all_forCoPred$LLID)
NHD_all_forCoPred$Wetland_P<-as.factor(NHD_all_forCoPred$Wetland_P)
NHD_all_forCoPred$Arealog<-log10(NHD_all_forCoPred$Area) #log transform drainage area variable

head(NHD_all_forCoPred)
names(NHD_all_forCoPred)
length(NHD_all_forCoPred$NHD_FID)

NHD_all_forCoPred2<-NHD_all_forCoPred
names(NHD_all_forCoPred2)
NHD_all_forCoPred2[,c(2:6,8,15)]<-scale(NHD_all_forCoPred2[,c(2:6,8,15)]) #standarize continous variables
head(NHD_all_forCoPred2)

colnames(NHD_all_forCoPred2)[c(3,8,15)]<-c("MAPZ","ElevZ","ArealogZ")

cohodatNBall2<-cohodatNBall
cohodatNBall2[,c(2:7,15,18)]
cohodatNBall2$ArealogZ<-zscore(cohodatNBall$Arealog,mean(NHD_all_forCoPred$Arealog),sd(NHD_all_forCoPred$Arealog))
cohodatNBall2$MAPZ<-zscore(cohodatNBall$mn_precip,mean(NHD_all_forCoPred$mn_precip),sd(NHD_all_forCoPred$mn_precip))
cohodatNBall2$ElevZ<-zscore(cohodatNBall$Elevation_m,mean(NHD_all_forCoPred$Elevation_m),sd(NHD_all_forCoPred$Elevation_m))
names(cohodatNBall2)


#--------------------------------------------------------------------------
#       Load NHD Data with segments upstream of steelhead barriers removed
#--------------------------------------------------------------------------
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

head(NHD_all_forSDPred)
names(NHD_all_forSDPred)
length(NHD_all_forSDPred$NHD_FID)

NHD_all_forSDPred2<-NHD_all_forSDPred
names(NHD_all_forSDPred2)
NHD_all_forSDPred2[,c(2:6,8,15)]<-scale(NHD_all_forSDPred2[,c(2:6,8,15)]) #standarize continous variables
head(NHD_all_forSDPred2)

colnames(NHD_all_forSDPred2)[c(4,8,15)]<-c("SlopeZ","ElevZ","ArealogZ")

STHDdatall2<-STHDdatall
names(STHDdatall2)
STHDdatall2$ArealogZ<-zscore(STHDdatall$Arealog,mean(NHD_all_forSDPred$Arealog),sd(NHD_all_forSDPred$Arealog))
STHDdatall2$ElevZ<-zscore(STHDdatall$Elevation_m,mean(NHD_all_forSDPred$Elevation_m),sd(NHD_all_forSDPred$Elevation_m))
STHDdatall2$SlopeZ<-zscore(STHDdatall$slope,mean(NHD_all_forSDPred$slope),sd(NHD_all_forSDPred$slope))
head(STHDdatall2)


#----------------------------------------------------------------------
#       Load NHD Data with segments upstream of chum barriers removed
#----------------------------------------------------------------------
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


#-------------------------------------
# fit best models identified in CH 1
#-------------------------------------
#Coho model
Best.modCo<-glmer(OBSCO_V_2~ArealogZ+MAPZ+ElevZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)),data = cohodatNBall2)
#Steelhead model
Best.modSD<-glmer(OBSSD_V_2~ArealogZ+ElevZ+SlopeZ+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),data = STHDdatall2)
#Chum model--note: no random effect since we did not collect data from all subbasins, thus not enouhg levels to use mixed model 
Best.modCH<-glm(OBSCH_V_2~ArealogZ+ElevZ+Geolgoy,family = binomial(link=logit),data = CHdatall2)


#-------------------------------------
# predict reach probabilities for coho model
#-------------------------------------
#covariates for best coho model include: mn_precip, Elevation_m, Arealog, Geolgoy, and Wetland_P
#Make sure these are all included in the predicitons dataframe below for newdat.CoMod object. 

newdat.CoMod<-data.frame(ArealogZ=NHD_all_forCoPred2[,"ArealogZ"],MAPZ=NHD_all_forCoPred2[,"MAPZ"],ElevZ=NHD_all_forCoPred2[,"ElevZ"],Geolgoy=NHD_all_forCoPred2[,"Geolgoy"],Wetland_P=NHD_all_forCoPred2[,"Wetland_P"],HUC_10=NHD_all_forCoPred2[,"HUC_10"])
p<-predictInterval(merMod = Best.modCo, newdata = newdat.CoMod, 
                   level = 0.95, n.sims = 1000,
                   stat = "median", which = "all", 
                   include.resid.var = FALSE)

p.Co<-as.vector(invlogit(p[,2][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect
p.Coupper<-as.vector(invlogit(p[,3][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
p.Colower<-as.vector(invlogit(p[,4][which(p[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
NHD_seg<-as.vector(NHD_all_forCoPred[,"NHD_FID"])#create a vector of NHD segment IDs
length(NHD_seg)
length(p.Co)
length(p.Coupper)
length(p.Colower)
COpred.all<-cbind(NHD_seg,p.Co,p.Colower,p.Coupper)
head(COpred.all)
tail(COpred.all)
COpred.all<-as.data.frame(COpred.all)
length(COpred.all$NHD_seg)

#check fitted values to make sure they are the same among models
fitted.values(Best.modCo_NoStand)[1]
fitted.values(Best.modCo)[1]
rowcheck2<-as.integer(rownames(subset(NHD_all_forCoPred,NHD_all_forCoPred$NHD_FID==cohodatNBall[1,"NHD_FID"])))
COpred.all[rowcheck2,]


NHD_all_forCoPred$pCo<-p.Co
NHD_all_forCoPred$p.CoUp<-p.Coupper
NHD_all_forCoPred$p.CoL<-p.Colower
head(NHD_all_forCoPred)

write.csv(COpred.all,"Co_predictions_NHDall_bestmod.csv") #create output file for coho predictions


#-------------------------------------
# predict reach probabilities for stelhead model
#-------------------------------------
#covariates for best steelhead model include: slope, Elevation_m, Arealog, Geolgoy, and Wetland_P
#Make sure these are all included in the predicitons dataframe below for newdat.CoMod object. 

Best.modSD_nostand<-glmer(OBSSD_V_2~Arealog+Elevation_m+slope+Wetland_P+Geolgoy+(1|HUC_10),family = binomial(link=logit),glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),data = STHDdatall)


newdat.SDMod<-data.frame(ArealogZ=NHD_all_forSDPred2[,"ArealogZ"],SlopeZ=NHD_all_forSDPred2[,"SlopeZ"],ElevZ=NHD_all_forSDPred2[,"ElevZ"],Geolgoy=NHD_all_forSDPred2[,"Geolgoy"],Wetland_P=NHD_all_forSDPred2[,"Wetland_P"],HUC_10=NHD_all_forSDPred2[,"HUC_10"])
pSD<-predictInterval(merMod = Best.modSD, newdata = newdat.SDMod, 
                   level = 0.95, n.sims = 1000,
                   stat = "median", which = "all", 
                   include.resid.var = FALSE)

p.SD<-as.vector(invlogit(pSD[,2][which(pSD[,1]=="combined")]))# create a vector of predicted values that include random effect
p.SDupper<-as.vector(invlogit(pSD[,3][which(pSD[,1]=="combined")]))# create a vector of predicted values that include random effect upper CI
p.SDlower<-as.vector(invlogit(pSD[,4][which(pSD[,1]=="combined")]))# create a vector of predicted values that include random effect lower CI
NHD_seg<-as.vector(NHD_all_forSDPred[,"NHD_FID"])#create a vector of NHD segment IDs
length(NHD_seg)
length(p.SD)
length(p.SDupper)
length(p.SDlower)
SDpred.all<-cbind(NHD_seg,p.SD,p.SDlower,p.SDupper)
head(SDpred.all)
tail(SDpred.all)
SDpred.all<-as.data.frame(SDpred.all)
length(SDpred.all$NHD_seg)

#check fitted values to make sure they are the same among models
fitted.values(Best.modSD_nostand)[1]
fitted.values(Best.modSD)[1]
rowcheck<-as.integer(rownames(subset(NHD_all_forSDPred,NHD_all_forSDPred$NHD_FID==STHDdatall[1,"NHD_FID"])))
SDpred.all[rowcheck,]

NHD_all_forSDPred$pSD<-p.SD
NHD_all_forSDPred$p.SDUp<-p.SDupper
NHD_all_forSDPred$p.SDL<-p.SDlower
write.csv(SDpred.all,"SD_predictions_NHDall_bestmod.csv") #create output file for steelhead predictions

#-------------------------------------
# predict reach probabilities for chum model
#-------------------------------------
#covariates for best steelhead model include: Elevation_m, Arealog,and Geolgoy
#no mixed effects included in this model since not all sub basins were surveyed
#since we are predictign outward we will not standardize the data. Fitted values are the same with and without standardization
# Best.modCH_noStand<-glm(OBSCH_V_2~Arealog+Elevation_m+Geolgoy,family = binomial(link=logit),data = CHdatall)

newdat.CH<-data.frame(ArealogZ=NHD_all_forCHPred2[,"ArealogZ"],ElevZ=NHD_all_forCHPred2[,"ElevZ"],Geolgoy=NHD_all_forCHPred2[,"Geolgoy"])

# Generate predicted values associated with the range of input values (FIXED effects model)
predict.CH <- predict(Best.modCH, newdat.CH, se.fit = T)
predict.CH <- data.frame(predict.CH$fit, predict.CH$se.fit)
colnames(predict.CH) <- c("fit", "se.fit")
length(predict.CH$fit)
predict.CH$fit.invlogit <- exp(predict.CH$fit)/(1+exp(predict.CH$fit))
predict.CH$upr.invlogit <- exp(predict.CH$fit+1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit+1.96*predict.CH$se.fit))
predict.CH$lwr.invlogit <- exp(predict.CH$fit-1.96*predict.CH$se.fit)/(1+exp(predict.CH$fit-1.96*predict.CH$se.fit))
NHD_seg<-as.vector(NHD_all_forCHPred[,"NHD_FID"])
predict.CH$NHD_seg<-NHD_seg
head(predict.CH)
p.CH<-predict.CH$fit.invlogit

write.csv(predict.CH,"CH_predictions_NHDall_bestmod_merged.csv") #create output file for chum predictions

NHD_all_forCHPred$pCH<-p.CH
NHD_all_forCHPred$p.CHUp<-predict.CH$upr.invlogit
NHD_all_forCHPred$p.CHL<-predict.CH$lwr.invlogit
#load functions needed for predicitons

source("E:\\Z167574\\R_data\\Chehalis upper extent\\Predictions\\PredictULOLocation.R")
source("E:\\Z167574\\R_data\\Chehalis upper extent\\Predictions\\OrderStreamProfiles.R")

#---------------------------------------
#    Predict Coho ULO Location
#---------------------------------------
COpred.all<-read.csv("Co_predictions_NHDall_bestmod.csv") #import file exported from text above. If running text continuously, # out this line. This dataset is provided to the user for expedniency. 
COpred_withMNs<-read.csv("CO_predictions_neighborhood_values.csv") #this is a dataset that joins above dataset that joins exported data for ArcMap
length(NHD_all_forCoPred[,1])
length(COpred_withMNs[,1]) #this has RASTERVALU which will need to be renamed "mn_pOC"
length(COpred.all[,1])

for(i in 1:length(COpred.all[,1])){
  NHD_all_forCoPred[which(NHD_all_forCoPred[,"NHD_FID"]==COpred.all[i,"NHD_seg"]),"pCo"]<-COpred.all[i,"p.Co"]
  #NHD_all_forCoPred[which(NHD_all_forCoPred[,"NHD_FID"]==COpred_withMNs[i,"NHD_FID"]),"Mn_pOC"]<-COpred_withMNs[i,"RASTERVALU"]
}

orderedNHD_allCo<-OrdermyStreamData(NHD_all_forCoPred) #he

#for mean neighborhood  pCO analysis----

CohoULOs0.5<-LocationofULOs_mn(orderedNHD_allCo,0.500000,"pCo","PCo_ULO")
write.csv(CohoULOs0.5,"Predicted_CohoULOs_glmm0.5_wNeighborhood.csv")
Co.5_Hab<-(sum(subset(CohoULOs0.5,CohoULOs0.5$Range=="Within")[,"Shape_Leng"]))/1000

CohoULOs0.75<-LocationofULOs_mn(orderedNHD_allCo,0.750000,"pCo","PCo_ULO")
write.csv(CohoULOs0.75,"Predicted_CohoULOs_glmm0.75_wNeighborhood.csv")
Co.75_Hab<-(sum(subset(CohoULOs0.75,CohoULOs0.75$Range=="Within")[,"Shape_Leng"]))/1000

CohoULOs0.25<-LocationofULOs_mn(orderedNHD_allCo,0.250000,"pCo","PCo_ULO")
write.csv(CohoULOs0.25,"Predicted_CohoULOs_glmm0.25_wNeighborhood.csv")
Co.25_Hab<-(sum(subset(CohoULOs0.25,CohoULOs0.25$Range=="Within")[,"Shape_Leng"]))/1000

cutpoint<-c(0.25,0.50,0.75)
estHab<-c(Co.25_Hab,Co.5_Hab,Co.75_Hab)
cohohabest<-as.data.frame(cbind(cutpoint,estHab))

Coho_3cut<-ggplot()+
  geom_point(data=cohohabest,aes(x=cutpoint,y=estHab),color="salmon", size = 2)+
  theme_classic()+
  xlab("cut-pont") + ylab("Rkm within Neighborhood")+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))

#---------------------------
# c<-seq(0.1,0.9,by=0.05)
# habestimate<-NULL
# for(i in c){
#  ULO<-LocationofULOs(orderedNHD_all,i,"p.Co","P.Co.ULO")
#  Hab<-(sum(subset(ULO,ULO$Range=="Within")[,"Shape_Leng"]))/1000
#  habestimate<-c(habestimate,Hab)
# }
# 
# 
# CohoRange.est<-as.data.frame(cbind(c,habestimate))
# 
# Coho_mulitCP<-ggplot()+
#   geom_point(data=CohoRange.est,aes(x=c,y=habestimate),color="salmon", size = 2)+
#   theme_classic()+
#   xlab("cut-pont") + ylab("Rkm within Neighborhood")+
#   theme(axis.title.x = element_text(size = 18))+
#   theme(axis.title.y = element_text(size = 18))+
#   theme(axis.text.x = element_text(size = 16))+
#   theme(axis.text.y = element_text(size = 16))
# 
# ppi <- 600 
# png("Coho_mulitCP.png", width=10, height=8, units = 'in', res = 600) #save as fi
# Coho_mulitCP
# dev.off()

#mean upstream pCO anaylsis-----

#PDT 0.10
CohoULOs0.1UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.100000,"pCo","PCo_ULO")
write.csv(CohoULOs0.1UM,"Predicted_CohoULOs_glmm0.1_wUpstreamMn.csv")
#PDT 0.15
CohoULOs0.15UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.150000,"pCo","PCo_ULO")
write.csv(CohoULOs0.15UM,"Predicted_CohoULOs_glmm0.15_wUpstreamMn.csv")
#PDT 0.20
CohoULOs0.2UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.200000,"pCo","PCo_ULO")
write.csv(CohoULOs0.2UM,"Predicted_CohoULOs_glmm0.2_wUpstreamMn.csv")
#---------------------------------------------------------------------------------------------
#PDT 0.25
CohoULOs0.25UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.250000,"pCo","PCo_ULO")
write.csv(CohoULOs0.25UM,"Predicted_CohoULOs_glmm0.25_wUpstreamMn.csv")
Co.25_HabUM<-(sum(subset(CohoULOs0.25UM,CohoULOs0.25UM$Range=="Within")[,"Shape_Leng"]))/1000
#---------------------------------------------------------------------------------------------
#PDT 0.30
CohoULOs0.3UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.300000,"pCo","PCo_ULO")
write.csv(CohoULOs0.3UM,"Predicted_CohoULOs_glmm0.30_wUpstreamMn.csv")
#PDT 0.35
CohoULOs0.35UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.350000,"pCo","PCo_ULO")
write.csv(CohoULOs0.35UM,"Predicted_CohoULOs_glmm0.35_wUpstreamMn.csv")
#PDT 0.40
CohoULOs0.40UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.400000,"pCo","PCo_ULO")
write.csv(CohoULOs0.40UM,"Predicted_CohoULOs_glmm0.40_wUpstreamMn.csv")
#PDT 0.45
CohoULOs0.45UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.400000,"pCo","PCo_ULO")
write.csv(CohoULOs0.45UM,"Predicted_CohoULOs_glmm0.45_wUpstreamMn.csv")
#---------------------------------------------------------------------------------------------
#PDT 0.50
CohoULOs0.5UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.500000,"pCo","PCo_ULO")
write.csv(CohoULOs0.5UM,"Predicted_CohoULOs_glmm0.5_wUpstreamMn.csv")
Co.5_HabUM<-(sum(subset(CohoULOs0.5UM,CohoULOs0.5UM$Range=="Within")[,"Shape_Leng"]))/1000
#---------------------------------------------------------------------------------------------
#PDT 0.55
CohoULOs0.55UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.550000,"pCo","PCo_ULO")
write.csv(CohoULOs0.55UM,"Predicted_CohoULOs_glmm0.55_wUpstreamMn.csv")
#PDT 0.60
CohoULOs0.60UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.600000,"pCo","PCo_ULO")
write.csv(CohoULOs0.60UM,"Predicted_CohoULOs_glmm0.60_wUpstreamMn.csv")
#PDT 0.65
CohoULOs0.65UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.650000,"pCo","PCo_ULO")
write.csv(CohoULOs0.65UM,"Predicted_CohoULOs_glmm0.65_wUpstreamMn.csv")
#PDT 0.70
CohoULOs0.70UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.700000,"pCo","PCo_ULO")
write.csv(CohoULOs0.70UM,"Predicted_CohoULOs_glmm0.70_wUpstreamMn.csv")
#---------------------------------------------------------------------------------------------
#PDT 0.75
CohoULOs0.75UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.750000,"pCo","PCo_ULO")
write.csv(CohoULOs0.75UM,"Predicted_CohoULOs_glmm0.75_wUpstreamMn.csv")
Co.75_HabUM<-(sum(subset(CohoULOs0.75UM,CohoULOs0.75UM$Range=="Within")[,"Shape_Leng"]))/1000
#---------------------------------------------------------------------------------------------
#PDT 0.80
CohoULOs0.80UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.800000,"pCo","PCo_ULO")
write.csv(CohoULOs0.80UM,"Predicted_CohoULOs_glmm0.80_wUpstreamMn.csv")
#PDT 0.85
CohoULOs0.85UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.850000,"pCo","PCo_ULO")
write.csv(CohoULOs0.85UM,"Predicted_CohoULOs_glmm0.85_wUpstreamMn.csv")
#PDT 0.90
CohoULOs0.90UM<-LocationofULOs_SPcheck(orderedNHD_allCo,0.900000,"pCo","PCo_ULO")
write.csv(CohoULOs0.90UM,"Predicted_CohoULOs_glmm0.90_wUpstreamMn.csv")

#--------------------------------------------------------------------------------------------------------

cutpointUM<-c(0.25,0.50,0.75)
estHabUM<-c(Co.25_HabUM,Co.5_HabUM,Co.75_HabUM)
(cohohabestUM<-as.data.frame(cbind(cutpointUM,estHabUM)))

Coho_3cutUM<-ggplot()+
  geom_point(data=cohohabestUM,aes(x=cutpointUM,y=estHabUM),color="salmon", size = 2)+
  theme_classic()+
  xlab("cut-pont") + ylab("Rkm within Neighborhood")+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))

#examine amount of habitat over range of cutpoint values
c<-seq(0.1,0.9,by=0.05)
habestimate<-NULL
for(i in c){
  ULO<-LocationofULOs_SPcheck(orderedNHD_allCo,i,"pCo","PCo_ULO")
  Hab<-(sum(subset(ULO,ULO$Range=="Within")[,"Shape_Leng"]))/1000
  habestimate<-c(habestimate,Hab)
 }


CohoRange.est<-as.data.frame(cbind(c,habestimate))

Coho_mulitCP<-ggplot()+
   geom_point(data=CohoRange.est,aes(x=c,y=habestimate),color="salmon", size = 2)+
   theme_classic()+
   xlab("cut-pont") + ylab("Rkm within Neighborhood")+
   theme(axis.title.x = element_text(size = 18))+
   theme(axis.title.y = element_text(size = 18))+
   theme(axis.text.x = element_text(size = 16))+
   theme(axis.text.y = element_text(size = 16))




#---------------------------------------------------
#compare results from method 1 and method 2
#---------------------------------------------------
comparedist<-cbind(cohohabestUM,cohohabest$estHab)
colnames(comparedist)<-c("cutpoint","upstreamMn","NB_Mn")

Coho_cutcomp<-ggplot()+
  geom_point(data=comparedist,aes(x=cutpoint,y=upstreamMn),color="salmon", size = 2)+
  geom_point(data=comparedist,aes(x=cutpoint,y=NB_Mn),color="red", size = 2)+
  theme_classic()+
  xlab("cut-pont") + ylab("Rkm within Neighborhood")+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))


#---------------------------------------
#    Predict Steelhead ULO Location
#---------------------------------------
SDpred.all<-read.csv("SD_predictions_NHDall_bestmod.csv")
length(NHD_all_forSDPred[,1])
length(SDpred.all[,1])

#this funcitno allows you to not have to run the prediced data combine the predicted values with the NHD reach values
for(i in 1:length(SDpred.all[,1])){
  NHD_all_forSDPred[which(NHD_all_forSDPred[,"NHD_FID"]==SDpred.all[i,"NHD_seg"]),"pSD"]<-SDpred.all[i,"p.SD"]
}


orderedNHD_allSD<-OrdermyStreamData(NHD_all_forSDPred) #
write.csv(orderedNHD_allSD,"SteelheadPredictions_FINAL.csv")

##mean upstream pSD anaylsis
#PDT 0.10
SteelheadULOs0.10UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.1000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.10UM,"Predicted_SteelheadULOs_glmm0.10_wUpstreamMn.csv")
#PDT 0.15
SteelheadULOs0.15UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.1500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.15UM,"Predicted_SteelheadULOs_glmm0.15_wUpstreamMn.csv")
#PDT 0.20
SteelheadULOs0.20UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.2000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.20UM,"Predicted_SteelheadULOs_glmm0.20_wUpstreamMn.csv")
#--------------------------------------------------------------------------------------------------------
#PDT 0.25
SteelheadULOs0.25UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.2500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.25UM,"Predicted_SteelheadULOs_glmm0.25_wUpstreamMn.csv")
SD.25_HabUM<-(sum(subset(SteelheadULOs0.25UM,SteelheadULOs0.25UM$Range=="Within")[,"Shape_Leng"]))/1000
#--------------------------------------------------------------------------------------------------------
#PDT 0.30
SteelheadULOs0.30UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.3000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.30UM,"Predicted_SteelheadULOs_glmm0.30_wUpstreamMn.csv")
#PDT 0.35
SteelheadULOs0.35UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.3500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.35UM,"Predicted_SteelheadULOs_glmm0.35_wUpstreamMn.csv")
#PDT 0.40
SteelheadULOs0.40UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.4000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.40UM,"Predicted_SteelheadULOs_glmm0.40_wUpstreamMn.csv")
#PDT 0.45
SteelheadULOs0.45UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.4500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.45UM,"Predicted_SteelheadULOs_glmm0.45_wUpstreamMn.csv")
#--------------------------------------------------------------------------------------------------------
#PDT 0.50
SteelheadULOs0.5UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.5UM,"Predicted_SteelheadULOs_glmm0.5_wUpstreamMn.csv")
SD.5_HabUM<-(sum(subset(SteelheadULOs0.5UM,SteelheadULOs0.5UM$Range=="Within")[,"Shape_Leng"]))/1000
#--------------------------------------------------------------------------------------------------------
#PDT 0.55
SteelheadULOs0.55UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.5500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.55UM,"Predicted_SteelheadULOs_glmm0.55_wUpstreamMn.csv")
#PDT 0.60
SteelheadULOs0.60UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.6000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.60UM,"Predicted_SteelheadULOs_glmm0.60_wUpstreamMn.csv")
#PDT 0.65
SteelheadULOs0.65UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.6500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.65UM,"Predicted_SteelheadULOs_glmm0.65_wUpstreamMn.csv")
#PDT 0.70
SteelheadULOs0.70UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.7000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.70UM,"Predicted_SteelheadULOs_glmm0.70_wUpstreamMn.csv")
#--------------------------------------------------------------------------------------------------------
#PDT 0.75
SteelheadULOs0.75UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.7500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.75UM,"Predicted_SteelheadULOs_glmm0.75_wUpstreamMn.csv")
SD.75_HabUM<-(sum(subset(SteelheadULOs0.75UM,SteelheadULOs0.75UM$Range=="Within")[,"Shape_Leng"]))/1000
#--------------------------------------------------------------------------------------------------------
#PDT 0.80
SteelheadULOs0.80UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.8000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.80UM,"Predicted_SteelheadULOs_glmm0.80_wUpstreamMn.csv")
#PDT 0.85
SteelheadULOs0.85UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.8500000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.85UM,"Predicted_SteelheadULOs_glmm0.85_wUpstreamMn.csv")
#PDT 0.90
SteelheadULOs0.90UM<-LocationofULOs_SPcheck(orderedNHD_allSD,0.9000000,"pSD","pSD_ULO")
write.csv(SteelheadULOs0.90UM,"Predicted_SteelheadULOs_glmm0.90_wUpstreamMn.csv")

#--------------------------------------------------------------------------------------------------------
cutpointUM<-c(0.25,0.50,0.75)
SDestHabUM<-c(SD.25_HabUM,SD.5_HabUM,SD.75_HabUM)
(SDhabestUM<-as.data.frame(cbind(cutpointUM,SDestHabUM)))


SD_3cutUM<-ggplot()+
  geom_point(data=SDhabestUM,aes(x=cutpointUM,y=SDestHabUM),color="grey53", size = 2)+
  theme_classic()+
  xlab("cut-pont") + ylab("Rkm within Neighborhood")+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))

Multispp_3cutUM<-ggplot()+
  geom_point(data=cohohabestUM,aes(x=cutpointUM,y=estHabUM),color="salmon", size = 2)+
  geom_point(data=SDhabestUM,aes(x=cutpointUM,y=SDestHabUM),color="grey53", size = 2)+
  theme_classic()+
  xlab("cut-pont") + ylab("Rkm within Neighborhood")+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))



#---------------------------------------
#    Predict Chum ULO Location
#---------------------------------------
CHpred.all<-read.csv("CH_predictions_NHDall_bestmod_merged.csv")
length(NHD_all_forCHPred[,1])
length(CHpred.all[,1])

#this funcitno allows you to not have to run the prediced data combine the predicted values with the NHD reach values
for(i in 1:length(CHpred.all[,1])){
  NHD_all_forCHPred[which(NHD_all_forCHPred[,"NHD_FID"]==CHpred.all[i,"NHD_seg"]),"pCH"]<-CHpred.all[i,"fit.invlogit"]
}


orderedNHD_allCH<-OrdermyStreamData(NHD_all_forCHPred) #

##mean upstream pSD anaylsis-----------------------------------------------------------------------------
#PDT 0.10
ChumULOs0.10UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.1000000,"pCH","pCH_ULO")
write.csv(ChumULOs0.10UM,"Predicted_ChumULOs_glmm0.10_wUpstreamMn.csv")
#PDT 0.15
ChumULOs0.15UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.1500000,"pCH","pCH_ULO")
write.csv(ChumULOs0.15UM,"Predicted_ChumULOs_glmm0.15_wUpstreamMn.csv")
#PDT 0.20
ChumULOs0.20UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.2000000,"pCH","pCH_ULO")
write.csv(ChumULOs0.20UM,"Predicted_ChumULOs_glmm0.20_wUpstreamMn.csv")
#--------------------------------------------------------------------------------------------------------
#PDT 0.25
ChumULOs0.25UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.2500000,"pCH","pCH_ULO")
write.csv(ChumULOs0.25UM,"Predicted_ChumULOs_glmm0.25_wUpstreamMn.csv")
CH.25_HabUM<-(sum(subset(ChumULOs0.25UM,ChumULOs0.25UM$Range=="Within")[,"Shape_Leng"]))/1000
#--------------------------------------------------------------------------------------------------------
#PDT 0.30
ChumULOs0.30UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.3000000,"pCH","pCH_ULO")
write.csv(ChumULOs0.30UM,"Predicted_ChumULOs_glmm0.30_wUpstreamMn.csv")
#PDT 0.35
ChumULOs0.35UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.3500000,"pCH","pCH_ULO")
write.csv(ChumULOs0.35UM,"Predicted_ChumULOs_glmm0.35_wUpstreamMn.csv")
#PDT 0.40
ChumULOs0.40UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.4000000,"pCH","pCH_ULO")
write.csv(ChumULOs0.40UM,"Predicted_ChumULOs_glmm0.40_wUpstreamMn.csv")
#PDT 0.45
ChumULOs0.45UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.4500000,"pCH","pCH_ULO")
write.csv(ChumULOs0.45UM,"Predicted_ChumULOs_glmm0.45_wUpstreamMn.csv")
#--------------------------------------------------------------------------------------------------------
#PDT 0.50
ChumULOs0.5UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.500000,"pCH","pCH_ULO")
write.csv(ChumULOs0.5UM,"Predicted_ChumULOs_glmm0.5_wUpstreamMn.csv")
CH.5_HabUM<-(sum(subset(ChumULOs0.5UM,ChumULOs0.5UM$Range=="Within")[,"Shape_Leng"]))/1000
#--------------------------------------------------------------------------------------------------------
#PDT 0.55
ChumULOs0.55UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.550000,"pCH","pCH_ULO")
write.csv(ChumULOs0.55UM,"Predicted_ChumULOs_glmm0.55_wUpstreamMn.csv")
#PDT 0.60
ChumULOs0.60UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.600000,"pCH","pCH_ULO")
write.csv(ChumULOs0.60UM,"Predicted_ChumULOs_glmm0.60_wUpstreamMn.csv")
#PDT 0.65
ChumULOs0.65UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.650000,"pCH","pCH_ULO")
write.csv(ChumULOs0.65UM,"Predicted_ChumULOs_glmm0.65_wUpstreamMn.csv")
#PDT 0.70
ChumULOs0.70UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.700000,"pCH","pCH_ULO")
write.csv(ChumULOs0.70UM,"Predicted_ChumULOs_glmm0.70_wUpstreamMn.csv")
#--------------------------------------------------------------------------------------------------------
#PDT 0.75
ChumULOs0.75UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.7500000,"pCH","pCH_ULO")
write.csv(ChumULOs0.75UM,"Predicted_ChumULOs_glmm0.75_wUpstreamMn.csv")
CH.75_HabUM<-(sum(subset(ChumULOs0.75UM,ChumULOs0.75UM$Range=="Within")[,"Shape_Leng"]))/1000
#--------------------------------------------------------------------------------------------------------
#PDT 0.80
ChumULOs0.80UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.800000,"pCH","pCH_ULO")
write.csv(ChumULOs0.80UM,"Predicted_ChumULOs_glmm0.80_wUpstreamMn.csv")
#PDT 0.85
ChumULOs0.85UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.850000,"pCH","pCH_ULO")
write.csv(ChumULOs0.85UM,"Predicted_ChumULOs_glmm0.85_wUpstreamMn.csv")
#PDT 0.90
ChumULOs0.90UM<-LocationofULOs_SPcheck(orderedNHD_allCH,0.900000,"pCH","pCH_ULO")
write.csv(ChumULOs0.90UM,"Predicted_ChumULOs_glmm0.90_wUpstreamMn.csv")


