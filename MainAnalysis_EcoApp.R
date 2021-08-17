#================================================================================
#        Predicting stream reaches within range of occurrence for coho, steelhead, and chum
#   
#     Author: Eric J Walther
#     Last date modified: Aug 11, 2021
#=================================================================================

#This script will take all the NHD stream reaches in the Chehalis and predict the probability that each reach is
#within the range of occurrence for the species of interest. First, each data set for fitting the best glmmm for each species needs to be imported.
#These models were identified during the analysis for CH1. Data sets are found in respect folders for each model and imported below. 

rm(list=ls(all=TRUE)) #clear workspace

library(ggplot2)
library(MuMIn)
library(merTools)
library(lme4)
library(dplyr)
library(ResourceSelection)
library(pROC)
library(lmtest)
library(tidyr)
library(Hmisc)
library(reshape2)

setwd(...) #add source to dataset

#function to calculate z-score values 
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

#merge predicted data to main coho NHD dataset
NHD_all_forCoPred$pCo<-p.Co
NHD_all_forCoPred$p.CoUp<-p.Coupper
NHD_all_forCoPred$p.CoL<-p.Colower
head(NHD_all_forCoPred)

#create output file for coho predictions
write.csv(COpred.all,"Co_predictions_NHDall_bestmod.csv") 

#-------------------------------------
# predict reach probabilities for stelhead model
#-------------------------------------
#covariates for best steelhead model include: slope, Elevation_m, Arealog, Geolgoy, and Wetland_P

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

#merge predicted data to main steelhead NHD dataset
NHD_all_forSDPred$pSD<-p.SD
NHD_all_forSDPred$p.SDUp<-p.SDupper
NHD_all_forSDPred$p.SDL<-p.SDlower

#create output file for steelhead predictions
write.csv(SDpred.all,"SD_predictions_NHDall_bestmod.csv") 

#-------------------------------------
# predict reach probabilities for chum model
#-------------------------------------
#covariates for best chum model include: Elevation_m, Arealog,and Geolgoy
#no mixed effects included in this model since not all sub basins were surveyed

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
#Note: user will need to set appriopraite pathway to where these scripts are stored
source("~\\PredictULOLocation.R")
source("~\\OrderStreamProfiles.R")

#------------------------------------------------
#    Prepare data to predict Coho ULO Location
#------------------------------------------------
#import file exported from text above. If running text continuously,
# "#" out line 325 below. This dataset is provided to the user for expediency. 
COpred.all<-read.csv("Co_predictions_NHDall_bestmod.csv") 

#this is a dataset that joins above dataset that joins exported data for ArcMap
COpred_withMNs<-read.csv("CO_predictions_neighborhood_values.csv") 

#function the append predicted values and neighborhood mean values to NHD dataset predicting ULO location using 
for(i in 1:length(COpred.all[,1])){
  NHD_all_forCoPred[which(NHD_all_forCoPred[,"NHD_FID"]==COpred.all[i,"NHD_seg"]),"pCo"]<-COpred.all[i,"p.Co"]
}

#order stream profiles sequentially in the downstream to upstream direction
#save object that stores sorted streams with predicitive values for each NHD segment.
#This object will be used for function below.
orderedNHD_allCo<-OrdermyStreamData(NHD_all_forCoPred)
#export data
write.csv(orderedNHD_allCo,"CohoPredictions_FINAL.csv")

#---------------------------------------
#    Prepare data to pedict Steelhead ULO Location
#---------------------------------------
SDpred.all<-read.csv("SD_predictions_NHDall_bestmod.csv")
length(NHD_all_forSDPred[,1])
length(SDpred.all[,1])

#this funcitno allows you to not have to run the prediced data combine the predicted values with the NHD reach values
for(i in 1:length(SDpred.all[,1])){
  NHD_all_forSDPred[which(NHD_all_forSDPred[,"NHD_FID"]==SDpred.all[i,"NHD_seg"]),"pSD"]<-SDpred.all[i,"p.SD"]
}

#save object that stores sorted streams with predicitive values for each NHD segment.
#This object will be used for function below.
orderedNHD_allSD<-OrdermyStreamData(NHD_all_forSDPred) 
#export data
write.csv(orderedNHD_allSD,"SteelheadPredictions_FINAL.csv")

#-----------------------------------------------
#    Prepare data to predict Chum ULO Location
#-----------------------------------------------
CHpred.all<-read.csv("CH_predictions_NHDall_bestmod_merged.csv")
length(NHD_all_forCHPred[,1])
length(CHpred.all[,1])

#this funcitno allows you to not have to run the prediced data combine the predicted values with the NHD reach values
for(i in 1:length(CHpred.all[,1])){
  NHD_all_forCHPred[which(NHD_all_forCHPred[,"NHD_FID"]==CHpred.all[i,"NHD_seg"]),"pCH"]<-CHpred.all[i,"fit.invlogit"]
}

#save object that stores sorted streams with predicitive values for each NHD segment.
#This object will be used for function below.
orderedNHD_allCH<-OrdermyStreamData(NHD_all_forCHPred) 
#export data
write.csv(orderedNHD_allCH,"ChumPredictions_FINAL.csv")

#----------------------------------
#         Objective One
#----------------------------------

#Function to predict ULO fro all stream profiles using multiple Probability thresholds
#The required fields include: 
#1) dataset: This is the ordered data set for all NHD segments that the model will predict on. 
#            This data frame is generated previously in the script that includes the predicted 
#            values using the best model. 
#2) cutpoint: This is the range of probability decision thresholds that the user would like to use.
#             For this analysis 0.10-0.90 at 0.05 intervals were exploered.
#3) species: This defines the columns in the data frame that contain the predicted values.
#4) prediction: This is the nameof the new column created that will idenfity if the NHD segment is
#               the ULO for a stream profile (binary output: TRUE/FALSE).
#5) filename: This is the name for the file that is exported after each loop of the function.

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

#define the terms used for the MULTITEST_ULO function above for each species 
filenameCo<-"Predicted_CohoULOs_glmm"
filenameSthd<-"Predicted_SthdULOs_glmm"
filenameCH<-"Predicted_ChumULOs_glmm"
cutpointUM<-seq(0.1,0.9,by=0.05)
checks<-c("pCo","p.CoUp","p.CoL")
checksSD<-c("pSD","p.SDUp","p.SDL")
checksCH<-c("pCH","p.CHUp","p.CHL")

#Run the function above for each species
Co_CI<-MULTITEST_ULO(orderedNHD_allCo,cutpointUM,checks,"PCo_ULO",filenameCo)
SD_CI<-MULTITEST_ULO(orderedNHD_allSD,cutpointUM,checksSD,"pSD_ULO",filenameSthd)
CH_CI<-MULTITEST_ULO(orderedNHD_allCH,cutpointUM,checksCH,"pCH_ULO",filenameCH)

#format dataframe of the functino output
CH_all<-cbind(CH_CI,cutpointUM)
colnames(CH_all)<-c("Mean","Upper","Lower","Cut")
SD_all<-cbind(SD_CI,cutpointUM)
colnames(SD_all)<-c("Mean","Upper","Lower","Cut")
Co_all<-cbind(Co_CI,cutpointUM)
colnames(Co_all)<-c("Mean","Upper","Lower","Cut")

#combine all data into one dataframe
allsppMulticut_CI<-rbind(Co_all,SD_all,CH_all)
allsppMulticut_CI$Species<-c(rep("Coho",17),rep("Steelhead",17),rep("Chum",17))
write.csv(allsppMulticut_CI,"PDT_values_Range_allspp.csv") #export data 

#Calculate Percent Change in coho range of occurrence between 0.25 and 0.75 PDTs
Co_all_sub<-subset(allsppMulticut_CI,allsppMulticut_CI$Species=="Coho")
Co_all_sub<-Co_all_sub[,-1]
PerChange(subset(Co_all_sub,Co_all_sub$Cut==0.25)[,1],subset(Co_all_sub,Co_all_sub$Cut==0.75)[,1])
subset(Co_all_sub,Co_all_sub$Cut==0.25)[,1]-subset(Co_all_sub,Co_all_sub$Cut==0.75)[,1]

#Calculate Percent Change in steelhead range of occurrence between 0.25 and 0.75 PDTs
SD_all_sub<-subset(allsppMulticut_CI,allsppMulticut_CI$Species=="Steelhead")
SD_all_sub<-SD_all_sub[,-1]
PerChange(subset(SD_all_sub,SD_all_sub$Cut==0.25)[,1],subset(SD_all_sub,SD_all_sub$Cut==0.75)[,1])
subset(SD_all_sub,SD_all_sub$Cut==0.25)[,1]-subset(SD_all_sub,SD_all_sub$Cut==0.75)[,1]

#Calculate Percent Change in steelhead range of occurrence between 0.25 and 0.75 PDTs
CH_all_sub<-subset(allsppMulticut_CI,allsppMulticut_CI$Species=="Chum")
CH_all_sub<-CH_all_sub[,-1]
PerChange(subset(CH_all_sub,CH_all_sub$Cut==0.25)[,1],subset(CH_all_sub,CH_all_sub$Cut==0.75)[,1])
subset(CH_all_sub,CH_all_sub$Cut==0.25)[,1]-subset(CH_all_sub,CH_all_sub$Cut==0.75)[,1]


###----------------------------###
###       Figure 2             ###
###----------------------------###

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

#coho summary----------------
Coho_SWIFD<-read.csv("Coho_SWIFD_NHDsegs_NoGH.csv")
Coho_SWIFD_NHD<-unique(Coho_SWIFD$NHD_FID)

Co.5within<-subset(CohoULOs0.5UM,CohoULOs0.5UM$Range=="Within")
Co.25within<-subset(CohoULOs0.25UM,CohoULOs0.25UM$Range=="Within")
Co.75within<-subset(CohoULOs0.75UM,CohoULOs0.75UM$Range=="Within")

Co.5within_NHD<-unique(Co.5within$NHD_FID)
Co.25within_NHD<-unique(Co.25within$NHD_FID)
Co.75within_NHD<-unique(Co.75within$NHD_FID)

#calculate the amount of habitat within the current described distribution
Co_SWIFD_length<-sum(Coho_SWIFD[,"Shape_Leng"])/1000

Co.25within_length<-sum(Co.25within[,"Shape_Leng"])/1000
Co.5within_length<-sum(Co.5within[,"Shape_Leng"])/1000
Co.75within_length<-sum(Co.75within[,"Shape_Leng"])/1000

#shared coho
Sharedwithin_COpred0.5<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%in%Co.5within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_COpred0.25<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%in%Co.25within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_COpred0.75<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%in%Co.75within_NHD)[,"Shape_Leng"]))/1000

#only within current described coho distribution
onlySWIFD_co0.5<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.5within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_co0.25<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.25within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_co0.75<-(sum(subset(Coho_SWIFD,Coho_SWIFD$NHD_FID%nin%Co.75within_NHD)[,"Shape_Leng"]))/1000

#only within predicted range of occurrence for coho
onlyPred_co0.5<-(sum(subset(Co.5within,Co.5within$NHD_FID%nin%Coho_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_co0.25<-(sum(subset(Co.25within,Co.25within$NHD_FID%nin%Coho_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_co0.75<-(sum(subset(Co.75within,Co.75within$NHD_FID%nin%Coho_SWIFD_NHD)[,"Shape_Leng"]))/1000

#steelhead summary--------------
SD_SWIFD<-read.csv("Sthd_SWIFD_NHDsegs_NoGH.csv")
SD_SWIFD_NHD<-unique(SD_SWIFD$NHD_FID)
SD.5within<-subset(SteelheadULOs0.5UM,SteelheadULOs0.5UM$Range=="Within")
SD.25within<-subset(SteelheadULOs0.25UM,SteelheadULOs0.25UM$Range=="Within")
SD.75within<-subset(SteelheadULOs0.75UM,SteelheadULOs0.75UM$Range=="Within")

SD.5within_NHD<-unique(SD.5within$NHD_FID)
SD.25within_NHD<-unique(SD.25within$NHD_FID)
SD.75within_NHD<-unique(SD.75within$NHD_FID)

#calculate the amount of habitat within the current described distribution
SD_SWIFD_length<-sum(SD_SWIFD[,"Shape_Leng"])/1000

#shared steelhead
Sharedwithin_SDpred0.5<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%in%SD.5within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_SDpred0.25<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%in%SD.25within_NHD)[,"Shape_Leng"]))/1000
Sharedwithin_SDpred0.75<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%in%SD.75within_NHD)[,"Shape_Leng"]))/1000

#only within current described steelhead distribution
onlySWIFD_SD0.5<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%nin%SD.5within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_SD0.25<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%nin%SD.25within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_SD0.75<-(sum(subset(SD_SWIFD,SD_SWIFD$NHD_FID%nin%SD.75within_NHD)[,"Shape_Leng"]))/1000

#only within predicted range of occurrence for sthd
onlyPred_SD0.5<-(sum(subset(SD.5within,SD.5within$NHD_FID%nin%SD_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_SD0.25<-(sum(subset(SD.25within,SD.25within$NHD_FID%nin%SD_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_SD0.75<-(sum(subset(SD.75within,SD.75within$NHD_FID%nin%SD_SWIFD_NHD)[,"Shape_Leng"]))/1000


#chum summmary--------------
Chum_SWIFD<-read.csv("Chum_SWIFD_NHDsegs_NoGH.csv")
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

#only within current described chum distribution
onlySWIFD_CH<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%nin%CH.5within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_CH.25<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%nin%CH.25within_NHD)[,"Shape_Leng"]))/1000
onlySWIFD_CH.75<-(sum(subset(Chum_SWIFD,Chum_SWIFD$NHD_FID%nin%CH.75within_NHD)[,"Shape_Leng"]))/1000

#only within predicted range of occurrence for chum
onlyPred_CH<-(sum(subset(CH.5within,CH.5within$NHD_FID%nin%CH_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_CH.25<-(sum(subset(CH.25within,CH.25within$NHD_FID%nin%CH_SWIFD_NHD)[,"Shape_Leng"]))/1000
onlyPred_CH.75<-(sum(subset(CH.75within,CH.75within$NHD_FID%nin%CH_SWIFD_NHD)[,"Shape_Leng"]))/1000


###----------------------###
###     Figure 3
###----------------------###
allspp_CI_melt<-melt(allspp_CI,id.vars = c("Species","Cutpoint"),measure.vars = c("Upper", "Lower","Mean"))

#coho
PCohocompV<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Coho")[,"value"]
CoPercDif<-(PCohocompV-Co_SWIFD_length)/((PCohocompV+Co_SWIFD_length)/2)*100
CoPerChng<-((PCohocompV-Co_SWIFD_length)/(Co_SWIFD_length))*100
CohoPreds<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Coho")
CohoPreds$PerDif<-CoPercDif
CohoPreds$PerChng<-CoPerChng
CohoPreds
CohoPreds$KmDiff<-CohoPreds$value-Co_SWIFD_length

#steelhead
PSTHDcompV<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Steelhead")[,"value"]
SthdPercDif<-(PSTHDcompV-SD_SWIFD_length)/((PSTHDcompV+SD_SWIFD_length)/2)*100
SthdoPerChng<-((PSTHDcompV-SD_SWIFD_length)/(SD_SWIFD_length))*100
SthdPreds<-subset(allspp_CI_melt,allspp_CI_melt$Species=="Steelhead")
SthdPreds$PerDif<-SthdPercDif
SthdPreds$PerChng<-SthdoPerChng
SthdPreds
SthdPreds$KmDiff<-SthdPreds$value-SD_SWIFD_length

#chum
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

#chum
ChumPreds.mean<-subset(ChumPreds,ChumPreds$variable=="Mean")
ChumPreds.mean<-ChumPreds.mean[,-c(4,5)]
spread.ChumPreds.mean<-ChumPreds.mean %>% spread(variable,PerChng)
ChumPreds.CI<-subset(ChumPreds,ChumPreds$variable!="Mean")
ChumPreds.CI<-ChumPreds.CI[,-c(4,5)]
spread.ChumPreds.CI<-ChumPreds.CI %>% spread(variable,PerChng)
chumPerDif<-cbind(spread.ChumPreds.CI,spread.ChumPreds.mean[,"Mean"])
colnames(chumPerDif)<-c("Species","Cutpoint","Upper","Lower","Mean")

#steelhead
SthdPreds.mean<-subset(SthdPreds,SthdPreds$variable=="Mean")
SthdPreds.mean<-SthdPreds.mean[,-c(4,5)]
spread.SthdPreds.mean<-SthdPreds.mean %>% spread(variable,PerChng)
SthdPreds.CI<-subset(SthdPreds,SthdPreds$variable!="Mean")
SthdPreds.CI<-SthdPreds.CI[,-c(4,5)]
spread.SthdPreds.CI<-SthdPreds.CI %>% spread(variable,PerChng)
sthdPerDif<-cbind(spread.SthdPreds.CI,spread.SthdPreds.mean[,"Mean"])
colnames(sthdPerDif)<-c("Species","Cutpoint","Upper","Lower","Mean")

#coho
CohoPreds.mean<-subset(CohoPreds,CohoPreds$variable=="Mean")
CohoPreds.mean<-CohoPreds.mean[,-c(4,5)]
spread.CohoPreds.mean<-CohoPreds.mean %>% spread(variable,PerChng)
CohoPreds.CI<-subset(CohoPreds,CohoPreds$variable!="Mean")
CohoPreds.CI<-CohoPreds.CI[,-c(4,5)]
spread.CohoPreds.CI<-CohoPreds.CI %>% spread(variable,PerChng)
cohoPerDif<-cbind(spread.CohoPreds.CI,spread.CohoPreds.mean[,"Mean"])
colnames(cohoPerDif)<-c("Species","Cutpoint","Upper","Lower","Mean")

#combine all species
PerDifAllspp<-rbind(cohoPerDif,sthdPerDif,chumPerDif)

#order and factor species
PerDifAllspp$Species<-factor(PerDifAllspp$Species, levels = c("Coho","Steelhead","Chum"),labels = c("Coho","Steelhead","Chum"))

#create figure
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
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.position = c(.9,.8))

#export figure
ppi <- 600 
tiff("Multispp_PerChng_Fig3.tiff", width=8, height=6, units = 'in', res = 600) #save as fi
Multispp_PerDif
dev.off()


