##-------------------------------------------------------------
##          Predict ULO location v.1
##
##  Note: This version of the function is designed to determine
##        the ULO location with the mean 1 km neighborhood probability
##        values as a control for identifying the ULO. 
##
#--------------------------------------------------------------

LocationofULOs_mn<-function(dataset,cut,species,prediction){
streamdat<-dataset
ListofStreams<-unique(streamdat[,"LLID"])
nList<-length(ListofStreams)
UpdatedNHDdat<-data.frame()
for(i in 1:nList){
StreamProfileForPred<-subset(dataset,(dataset[,"LLID"]%in%ListofStreams[i]))#DATA FIT MOD
compR<-NULL
nC<-length(StreamProfileForPred[,species])#number of rows (i.e stream reaches) being evaluated 
if(anyNA(StreamProfileForPred[,species])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
for (j in 1:nC) compR[j]<- {
  if(any(StreamProfileForPred[1:j,species]<cut)=="FALSE"){compR[j]<-1}
  else{
    if((any(StreamProfileForPred[j:nC,species]>cut)=="TRUE")&(StreamProfileForPred[j,"Mn_pOC"]>cut)){compR[j]<-1}else{compR[j]<-0
    break}}}# value t
segCount<-sum(compR)
if(segCount>0){StreamProfileForPred[segCount,prediction]<-"TRUE";StreamProfileForPred[which(!(StreamProfileForPred[,prediction]%in%"TRUE")),prediction]<-"FALSE"}else{StreamProfileForPred[,prediction]<-"FALSE"}
if(segCount>0){StreamProfileForPred[1:segCount,"Range"]<-"Within";StreamProfileForPred[which(!(StreamProfileForPred[,"Range"]%in%"Within")),"Range"]<-"Outside"}else{StreamProfileForPred[,"Range"]<-"Outside"}
UpdatedNHDdat<-rbind(UpdatedNHDdat,StreamProfileForPred)
}
return(UpdatedNHDdat)
}
##-------------------------------------------------------------
##          Predict ULO location v.2
##
##  Note: This version of the function is designed to determine
##        the ULO location based on condition of upstream reaches
##        from the focal reach based along an indivual stream profile
##
#--------------------------------------------------------------
LocationofULOs_SPcheck<-function(dataset,cut,species,prediction){
  streamdat<-dataset
  ListofStreams<-unique(streamdat[,"LLID"])
  nList<-length(ListofStreams)
  UpdatedNHDdat<-data.frame()
  for(i in 1:nList){
    StreamProfileForPred<-subset(dataset,(dataset[,"LLID"]%in%ListofStreams[i]))#DATA FIT MOD
    compR<-NULL
    nC<-length(StreamProfileForPred[,species])#number of rows (i.e stream reaches) being evaluated 
    if(anyNA(StreamProfileForPred[,species])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
    for (j in 1:nC) compR[j]<- {
      if(any(StreamProfileForPred[1:j,species]<cut)=="FALSE"){compR[j]<-1}
      else{
        if((StreamProfileForPred[j,"order"]+5)>max(StreamProfileForPred[,"order"])){
          if(StreamProfileForPred[j,"order"]!=max(StreamProfileForPred[,"order"])){
            if(mean(StreamProfileForPred[(j+1):max(StreamProfileForPred[,"order"]),species])>cut){compR[j]<-1}else{compR[j]<-0
            break}}
          else{if(StreamProfileForPred[j,species]>cut){compR[j]<-1}else{compR[j]<-0
          break}}}
        else{if(mean(StreamProfileForPred[(j+1):(j+5),species])>cut){compR[j]<-1}else{compR[j]<-0
        break}}
      }}# value t
    segCount<-sum(compR)
    if(segCount>0){StreamProfileForPred[segCount,prediction]<-"TRUE";StreamProfileForPred[which(!(StreamProfileForPred[,prediction]%in%"TRUE")),prediction]<-"FALSE"}else{StreamProfileForPred[,prediction]<-"FALSE"}
    if(segCount>0){StreamProfileForPred[1:segCount,"Range"]<-"Within";StreamProfileForPred[which(!(StreamProfileForPred[,"Range"]%in%"Within")),"Range"]<-"Outside"}else{StreamProfileForPred[,"Range"]<-"Outside"}
    UpdatedNHDdat<-rbind(UpdatedNHDdat,StreamProfileForPred)
  }
  return(UpdatedNHDdat)
}

##-------------------------------------------------------------
##          Predict ULO location v.3
##
##  Note: This version of the function is designed to determine
##        the ULO location based on the last consecutive sttream
##        reach that is above the cutpoint 
##
#--------------------------------------------------------------

LocationofULOs<-function(dataset,cut,species,prediction){
  streamdat<-dataset
  ListofStreams<-unique(streamdat[,"LLID"])
  nList<-length(ListofStreams)
  UpdatedNHDdat<-data.frame()
  for(i in 1:nList){
    StreamProfileForPred<-subset(dataset,(dataset[,"LLID"]%in%ListofStreams[i]))#DATA FIT MOD
    compR<-NULL
    nC<-length(StreamProfileForPred[,species])#number of rows (i.e stream reaches) being evaluated 
    if(anyNA(StreamProfileForPred[,species])=="FALSE"){print("no NAs predicted")}else{print("predicted NAs, error!!! check code");break}
    for (j in 1:nC) compR[j]<- {
      if(any(StreamProfileForPred[1:j,species]<cut)=="FALSE"){compR[j]<-1}
      else{compR[j]<-1}}# value t
    segCount<-sum(compR)
    if(segCount>0){StreamProfileForPred[segCount,prediction]<-"TRUE";StreamProfileForPred[which(!(StreamProfileForPred[,prediction]%in%"TRUE")),prediction]<-"FALSE"}else{StreamProfileForPred[,prediction]<-"FALSE"}
    if(segCount>0){StreamProfileForPred[1:segCount,"Range"]<-"Within";StreamProfileForPred[which(!(StreamProfileForPred[,"Range"]%in%"Within")),"Range"]<-"Outside"}else{StreamProfileForPred[,"Range"]<-"Outside"}
    UpdatedNHDdat<-rbind(UpdatedNHDdat,StreamProfileForPred)
  }
  return(UpdatedNHDdat)
}
