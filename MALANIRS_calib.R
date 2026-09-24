library(rchemo)
library(nirsextra)
library(readxl)
library(spectrolab)
library(stringr)
source("MALANIRS_list_pre.R")
source("/home/ecarnot/Documents/INRA/Projets/VitaSPEC/vitaspec_R/vitaspec_preCV.R")

## Bioch predicted from DIASCOPE post-harvest chain on Minelandiv2024 maize samples
## and spectra collected on AGAP NIR spectrometer (MalaNIRS project)

# Read AGAP
xAGAP=naturaspec2df("/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ech_DIASCOPE/MLD24DIASCOPE/")
xAGAP$Plot=as.numeric(str_extract(rownames(xAGAP), "(?<=_)\\d+(?=_)"))
sp=aggregate(xAGAP,by=list(xAGAP$Plot), mean)
sp$x=as.matrix(aggregate(xAGAP$x, list(xAGAP$Plot), mean)[,-1])
sp$Plot=sp$Group.1; sp=sp[,-1]

# Read data from Diascope
dat=read_xlsx("/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ech_DIASCOPE/MineLandDiv_T1.2_Trials 2024_INRAE Mauguio_20240729.xlsx", sheet="Data")
dat$Plot=as.numeric(substr(dat$Plot, 6, nchar(dat$Plot)))  

dat=merge(sp,dat,by="Plot")
traits_bioch=colnames(dat)[99] #[c(1,4:6,8:104)]
dat$MLD_DateH=as.Date(dat$MLD_DateH, format="%d/%m/%Y")

n_traits <- length(traits_bioch)
ncols   <- 4
nrows   <- ceiling(n_traits / ncols)
recap <- vector("list", n_traits)
names(recap) <- traits_bioch

pdf("calib/calibration_MLD24_spAGAP.pdf", width = ncols * 4, height = nrows * 4)
par(mfrow = c(nrows, ncols))
for (trait in traits_bioch) {
  cat("R2 CV :", trait, "\n")
  iout=is.na(dat[[trait]])
  if (length(unique(dat[[trait]])) == 1 || length(which(!iout)) < 20) {
    cat("Trait", trait, ": no variance, skipping.\n")
    next
  }
  
  recap[[trait]] <- vitaspec_preCV(dat$x[!iout, ], as.numeric(dat[[trait]][!iout]), list_pre = list_pre, ncomp = 20, titl = trait,
                                   plotLV = TRUE, plotYY = TRUE, verb = FALSE)
}
dev.off()


# Test avec autres spectres
pprot=rbind(list('red',c(500,10,1)),list('snv',''))
fm <- plskern(pre(dat$x[!iout, ],pprot), as.numeric(dat[[trait]][!iout]), nlv = 14)
xnew=naturaspec2df("/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ech_DIASCOPE/MLD24DIASCOPE/test/")
pred=predict(fm, pre(xnew$x,pprot), ncomp = 14)
plot(pred$pred, as.numeric(dat[[trait]][!iout]), xlab="Predicted", ylab="Observed", main=paste("Trait:", trait))

stop()

## From 20 ring test samples measured for biochemistry at CREA, try a calibration

# Calib
ref=read.table("/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/ring_test_bioch.csv", sep=";", header=TRUE, dec=".")
colnames(ref)[1]="ech"
# spA=sp2df(xAGAPrm)
# spA$ech=rownames(xAGAPrm)
# spD=sp2df(xDIAr)
# spD$ech=rownames(xDIAr)
# sp=merge(spA,spD,by="ech")
xAGAPm=aggregate(xAGAP,by=list(rownames(xAGAP)), mean)
rownames(xAGAPm)=xAGAPm[,1]

# DIASCOPE
# sp=sp2df(xDIA)
# sp$ech=rownames(xDIA)
# AGAP
sp=sp2df(as.matrix(xAGAPm[,-1]))
sp$ech=rownames(xAGAPm)

sp=merge(sp,ref,by="ech")

sp$x=sp$x
source("/home/ecarnot/Documents/INRA/Projets/VitaSPEC/vitaspec_R/list_pre.R")
ag=colnames(sp)[4:13]



