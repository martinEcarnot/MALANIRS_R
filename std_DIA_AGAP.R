# Stadardization between DIASCOPE SPECTRO and NATURASPEC AGAP

library(spectrolab)
library(rchemo)
library(nirsextra)
source("MALANIRS_utils.R")

# Read DIASCOPE
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ech_DIASCOPE/mld24dia.txt"
spD=read.table(f,skip = 11,header =TRUE,sep=",")
spD$ech=spD$Sample.Number #substr(spD$Sample.Number,6,12)
xDIA=as.matrix(spD[,-c(1:6,ncol(spD),ncol(spD)-1)])
row.names(xDIA)=make.unique(spD$ech)

# Read AGAP
xAGAP=NULL
ld=list.dirs("/media/ecarnot/C8DE-2338/MALANIRS/MALANIRS/MLD24DIASCOPE/")
sp=read_spectra(ld)
sp$ech=substr(sp$names,1,12)
# sp$ech=gsub("-","_",sp$ech)
xAGAP=sp$value
row.names(xAGAP)=sp$ech
colnames(xAGAP)=sp$bands
xAGAP=aggregate(xAGAP,by=list(sp$ech), mean)
rownames(xAGAP)=xAGAP[,1]
xAGAP[,1]=NULL

# 
com=intersect(rownames(xAGAP),rownames(xDIA))
xAGAPr=xAGAP[match(com,rownames(xAGAP)),51:750]
xDIAr=xDIA[match(com,rownames(xDIA)),seq(1,ncol(xDIA),2)]


# PDS from MALANIRS_rad_ring_tets
xDIApred<-as.matrix(xDIAr)%*%as.matrix(mPDS$P)
xDIApred<-sweep(xDIApred, 2, as.numeric(t(mPDS$Intercept)), "+")
plotsp(xDIApred)

