#====================================================
#        Functions for ordering streams
#   
#     Author: Eric J Walther
#     Last date modified: August 19, 2020
#===================================================

# These functions below are written for identifying and sortign stream profiles from
# the "NHD_wAttsStreamID_08122020_R.csv" dataset. This will be used to predict ULO locations along
# each stream profile. **WARNING**: there are column names that are specific this dataset used for this analysis (dataset mentioned above)
# to reuse or modify this code, the column names used in the data frams need to be checked and modified if necessary

#function ot sort streams
SortSream<-function(x){
  if(length(x[,1]>1)){
    d<-max(x[,"Elevation_m"])-min(x[,"Elevation_m"])
    if(d>2){
      if(x[which(x[,"R_Meas"]==(min(x[,"R_Meas"]))),"Arealog"]>x[which(x[,"R_Meas"]==(max(x[,"R_Meas"]))),"Arealog"]&x[which(x[,"R_Meas"]==(min(x[,"R_Meas"]))),"Elevation_m"]<x[which(x[,"R_Meas"]==(max(x[,"R_Meas"]))),"Elevation_m"]){
        x<-x[order(x[,"R_Meas"]),] #column 6 is the Rmeas column of measurement of NHHD along route
      }else{x<-x[order(x[,"R_Meas"],decreasing = TRUE),]}
      n<-nrow(x)#identify number of rows in dataframe
      x$order<-seq(1:n) 
    }
    else{if(x[which(x[,"R_Meas"]==(min(x[,"R_Meas"]))),"Arealog"]>x[which(x[,"R_Meas"]==(max(x[,"R_Meas"]))),"Arealog"]){
      x<-x[order(x[,"R_Meas"]),] #column 6 is the Rmeas column of measurement of NHHD along route
    }
      else{x<-x[order(x[,"R_Meas"],decreasing = TRUE),]}
      n<-nrow(x)#identify number of rows in dataframe
      x$order<-seq(1:n) 
    }
  }
  else{
    n<-nrow(x)#identify number of rows in dataframe
    x$order<-seq(1:n) 
  }
  return(x)
}
#function to order the entire dataframe of sampled streams
OrdermyStreamData<-function(x){
  sortedData<-data.frame()
  EachStreamID<-as.vector(unique(x$LLID))
  nstreams<-length(EachStreamID)
  for(i in 1:nstreams){
    dat<-subset(x,x$LLID==EachStreamID[i])
    Sortedstream<-SortSream(dat)
    sortedData<-rbind(sortedData,Sortedstream)
  }
  return(sortedData)
}

# #test on a small subset of data to make sure it works
# uniqueStreamsfortest<-as.vector(unique(NHD_all2$LLID)[1:5])
# 
# which()
# 
# testsetForsortingLLIDs<-subset(NHD_all2,NHD_all2$LLID%in%uniqueStreamsfortest)
# testsetForsortingLLIDs<-testsetForsortingLLIDs[,c(9,11,12,15)]
# 
# 


