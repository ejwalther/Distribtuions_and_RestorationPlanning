library(ggplot2)
library(cowplot)

setwd("...")#set working directory to appropraite location where downloaded data are stored

GM_Flow<-read.csv("ChehalisGrandMound_TimeSeries.csv") 

GM_Flow$datetime<- as.Date(GM_Flow$datetime, "%m/%d/%Y")
GM_Flow$Year<-format(GM_Flow$datetime,"%Y")
GM_Flow$cfs<-as.numeric(GM_Flow$cfs)
GM_Flow$date<-as.Date(GM_Flow$datetime, "%m/%d/%Y")
GM_Flow$month_yr<-format(as.Date(GM_Flow$date), "%Y-%m")
GM_Flow$Month<-format(as.Date(GM_Flow$date), "%m")

monthlymeans<-as.data.frame(GM_Flow%>%
                              group_by(month_yr)%>%
                              summarise(avg.cfs=mean(cfs),
                                        sd.cfs=sd(cfs),
                                        month=Month[1],
                                        Year=Year[1]))
monthlymeans$month_count<-seq(1:length(monthlymeans$month_yr))

#------------------------
#     Coho
#------------------------
monthlymeansCoSpawn<-subset(monthlymeans,monthlymeans$month%in%c("11","12","01"))
monthlymeansCoSpawn$month_count<-seq(1:length(monthlymeansCoSpawn$month_yr))
monthlymeansCoSpawn<-monthlymeansCoSpawn[-c(277,278),]
index<-rep(1:(length(monthlymeansCoSpawn$avg.cfs)/3),3)
index <- index[order(index)]
Coyearmean<-aggregate(x=monthlymeansCoSpawn$avg.cfs, by = list(index), FUN=mean)

meanCo<-NULL
for(i in 1:length(Coyearmean$x)){
  dat<-rep(Coyearmean[i,2],3)
  meanCo<-c(meanCo,dat)
}

  
monthlymeansCoSpawn$mean<-meanCo
monthlymeansCoSpawnPlot<-subset(monthlymeansCoSpawn,monthlymeansCoSpawn$month%in%c("11"))
monthlymeansCoSpawnPlot$Year<-as.numeric(monthlymeansCoSpawnPlot$Year)


GMFLOW2020<-tail(subset(monthlymeans,monthlymeans$month%in%c("11","12","01")))
GMFLOW2020<-GMFLOW2020[5:6,]
GMFLOW2020$mean<-mean(GMFLOW2020$avg.cfs)
GMFLOW2020$Year<-as.numeric(GMFLOW2020$Year)

monthlymeansCoSpawnPlot<-rbind(monthlymeansCoSpawnPlot,GMFLOW2020[1,])
monthlymeansCoSpawnPlot$Year<-as.numeric(monthlymeansCoSpawnPlot$Year)

median(monthlymeansCoSpawnPlot$mean)

#plot coho 
CoFlowGM<-ggplot()+
  geom_line(data=monthlymeansCoSpawnPlot,aes(x=Year,y=mean),size=1)+
  geom_point(data=subset(monthlymeansCoSpawnPlot,monthlymeansCoSpawnPlot$Year=="2017"),aes(x=Year,y=mean),color="red",size=3)+
  geom_point(data=subset(monthlymeansCoSpawnPlot,monthlymeansCoSpawnPlot$Year=="2018"),aes(x=Year,y=mean),color="gold",size=3)+
  scale_y_continuous(limits = c(0,12500),breaks = seq(0,12500,by=2500))+
  scale_x_continuous(breaks=seq(1930,2020,by=15))+
  xlab(label="Year")+ylab(label = "Discharge (cfs)")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  geom_hline(aes(yintercept=median(monthlymeansCoSpawnPlot$mean)),colour="blue",linetype="dashed",size=1.5)
  
#------------------------
#     Steelhead
#------------------------
monthlymeansSthdSpawn<-subset(monthlymeans,monthlymeans$month%in%c("02","03","04","05"))
monthlymeansSthdSpawn$month_count<-seq(1:length(monthlymeansSthdSpawn$month_yr))

index<-rep(1:(length(monthlymeansSthdSpawn$avg.cfs)/4),4)
index <- index[order(index)]
SDyearmean<-aggregate(x=monthlymeansSthdSpawn$avg.cfs, by = list(index), FUN=mean)

meanSD<-NULL
for(i in 1:length(SDyearmean$x)){
  dat<-rep(SDyearmean[i,2],4)
  meanSD<-c(meanSD,dat)
}

monthlymeansSthdSpawn$mean<-meanSD
monthlymeansSthdSpawnPlot<-subset(monthlymeansSthdSpawn,monthlymeansSthdSpawn$month%in%c("02"))
monthlymeansSthdSpawnPlot$Year<-as.numeric(monthlymeansSthdSpawnPlot$Year)

#plot steelhead
SDFlowGM<-ggplot()+
  geom_line(data=monthlymeansSthdSpawnPlot,aes(x=Year,y=mean),size=1)+
  geom_hline(aes(yintercept=median(monthlymeansSthdSpawnPlot$mean)),colour="blue",linetype="dashed",size=1.5)+
  geom_point(data=subset(monthlymeansSthdSpawnPlot,monthlymeansSthdSpawnPlot$Year=="2018"),aes(x=Year,y=mean),color="red",size=3)+
  geom_point(data=subset(monthlymeansSthdSpawnPlot,monthlymeansSthdSpawnPlot$Year=="2019"),aes(x=Year,y=mean),color="gold",size=3)+
  scale_y_continuous(limits = c(0,12500),breaks = seq(0,12500,by=2500))+
  scale_x_continuous(breaks=seq(1930,2020,by=15))+
  xlab(label="Year")+ylab(label = "Discharge (cfs)")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))

#------------------------
#     Chum
#------------------------
monthlymeansChumSpawn<-subset(monthlymeans,monthlymeans$month%in%c("11"))
monthlymeansChumSpawn$month_count<-seq(1:length(monthlymeansChumSpawn$month_yr))
monthlymeansChumSpawn$Year<-as.numeric(monthlymeansChumSpawn$Year)

CHFlowGM<-ggplot()+
  geom_line(data=monthlymeansChumSpawn,aes(x=Year,y=avg.cfs),size=1)+
  geom_hline(aes(yintercept=median(monthlymeansChumSpawn$avg.cfs)),colour="blue",linetype="dashed",size=1.5)+
  geom_point(data=subset(monthlymeansChumSpawn,monthlymeansChumSpawn$Year=="2017"),aes(x=Year,y=avg.cfs),color="red",size=3)+
  geom_point(data=subset(monthlymeansChumSpawn,monthlymeansChumSpawn$Year=="2018"),aes(x=Year,y=avg.cfs),color="gold",size=3)+
  scale_y_continuous(limits = c(0,12500),breaks = seq(0,12500,by=2500))+
  scale_x_continuous(breaks=seq(1930,2020,by=15))+
  xlab(label="Year")+ylab(label = "Discharge (cfs)")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))


#--------------------------------------
#     Combine plots for Figure S2
#--------------------------------------
all_flow<-plot_grid(CoFlowGM,SDFlowGM,CHFlowGM,ncol=2,align="v",labels = c("a)", "b)", "c)"))
ppi <- 600 
tiff(file = "Allspp_flow_Chehalis_GM.tiff", width=14*ppi, height=10*ppi, res=ppi) 
all_flow
dev.off() 

