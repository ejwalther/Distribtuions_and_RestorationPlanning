library(ggplot2)
library(cowplot)

setwd("...") #set working directory to appropraite location where downloaded data are stored


#Coho Escapement
CohoEsc<-read.csv("Coho_escapement.csv")

CohoEscPlot<-ggplot()+
  geom_line(data=CohoEsc,aes(x=Year,y=Escapement),size=1)+
  geom_point(data=subset(CohoEsc,CohoEsc$Year=="2017"),aes(x=Year,y=Escapement),color="red",size=3)+
  geom_point(data=subset(CohoEsc,CohoEsc$Year=="2018"),aes(x=Year,y=Escapement),color="gold",size=3)+
  scale_y_continuous(limits = c(0,110050),breaks = seq(10000,110000,by=20000))+
  scale_x_continuous(limits = c(1968,2020),breaks=seq(1970,2015,by=15))+
  xlab(label="Year")+ylab(label = "Escapement")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  geom_hline(aes(yintercept=median(na.omit(CohoEsc$Escapement))),colour="blue",linetype="dashed",size=1.5)

#Steelhead Escapement
SthdEsc<-read.csv("Steelhead_escapement.csv")

SthdEscPlot<-ggplot()+
  geom_line(data=SthdEsc,aes(x=Year,y=Escapement),size=1)+
  geom_point(data=subset(SthdEsc,SthdEsc$Year=="2018"),aes(x=Year,y=Escapement),color="red",size=3)+
  geom_point(data=subset(SthdEsc,SthdEsc$Year=="2019"),aes(x=Year,y=Escapement),color="gold",size=3)+
  scale_y_continuous(limits = c(0,110050),breaks = seq(10000,110000,by=20000))+
  scale_x_continuous(limits = c(1968,2020),breaks=seq(1970,2015,by=15))+
  xlab(label="Year")+ylab(label = "Escapement")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  geom_hline(aes(yintercept=median(na.omit(SthdEsc$Escapement))),colour="blue",linetype="dashed",size=1.5)

#Chum Escapement
ChumEsc<-read.csv("Chum_escapement.csv")

ChumEscPlot<-ggplot()+
  geom_line(data=ChumEsc,aes(x=Year,y=Escapement),size=1)+
  geom_point(data=subset(ChumEsc,ChumEsc$Year=="2017"),aes(x=Year,y=Escapement),color="red",size=3)+
  geom_point(data=subset(ChumEsc,ChumEsc$Year=="2018"),aes(x=Year,y=Escapement),color="gold",size=3)+
  scale_y_continuous(limits = c(0,110050),breaks = seq(10000,110000,by=20000))+
  scale_x_continuous(limits = c(1968,2020),breaks=seq(1970,2015,by=15))+
  xlab(label="Year")+ylab(label = "Escapement")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  geom_hline(aes(yintercept=median(na.omit(ChumEsc$Escapement))),colour="blue",linetype="dashed",size=1.5)

#--------------------------------------
#     Combine plots for Figure S3
#--------------------------------------

AllSpp_EscPlot<-plot_grid(CohoEscPlot,SthdEscPlot,ChumEscPlot,ncol=2,align="v",labels = c("a)", "b)", "c)"))
ppi <- 600 
tiff(file = "Allspp_Esc_Chehalis.tiff", width=14*ppi, height=10*ppi, res=ppi) 
AllSpp_EscPlot
dev.off() 
