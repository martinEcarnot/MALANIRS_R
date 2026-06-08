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
sp$ech=substr(sp$names,1,13)
# sp$ech=gsub("-","_",sp$ech)
xAGAP1=sp$value
row.names(xAGAP1)=sp$ech
colnames(xAGAP1)=sp$bands
