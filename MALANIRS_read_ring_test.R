library(spectrolab)
library(rchemo)
library(nirsextra)
source("~/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/PDS.R")

# # Commande bash pour lister les fichiers
# cd /home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/
# find . -type f -name "*.sed" -exec basename {} \; | sed 's/..........$//' | sort -u > list_done_Naturaspec.txt
# grep "MLD" list_done_Naturaspec.txt | awk '{num=substr($0, length($0)-2); print num, $0}' | sort -n | cut -d' ' -f2- > list_done_Naturaspec_MLD.txt
# grep "MRS" list_done_Naturaspec.txt | awk '{num=substr($0, length($0)-2); print num, $0}' | sort -n | cut -d' ' -f2- > list_done_Naturaspec_MRS.txt
# grep -v -e "MLD" -e "MRS" list_done_Naturaspec.txt > list_done_Naturaspec_autres.txt

# Fonction pour tracer courbes 
plotspMALA <-function(x1,x2,unit,n1,n2,tit=""){
  x=rbind(x1,x2)
  col=rep("red",nrow(x))
  col[1:nrow(x1)]="black"
  plotsp(x,col=col, xlab="wavelength (nm)",ylab=unit, lwd=3)
  legend(x="topright",legend=c(n1,n2), col = c("red","black"),lty=1)
  title(tit)
}

xAGAP=NULL
ld=list.dirs("/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/2025_Jun_03-mais_MLD")
sp=read_spectra(ld)
sp$ech=substr(sp$names,1,nchar(sp$names)-10)
sp$ech=gsub("-","_",sp$ech)
xAGAP1=sp$value
row.names(xAGAP1)=sp$ech
colnames(xAGAP1)=sp$bands
xAGAP=rbind(xAGAP,xAGAP1)
  
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/mala_DIA.txt"
spD=read.table(f,skip = 11,header =TRUE,sep=",")
spD$ech=substr(spD$Sample.Number,6,12)
xDIA=as.matrix(spD[,-c(1:6,ncol(spD),ncol(spD)-1)])
row.names(xDIA)=make.unique(spD$ech)
colnames(xDIA)=substr(colnames(spD[,-c(1:6,ncol(spD),ncol(spD)-1)]),2,8)
# save(xAGAP,xDIA, file = "malanirs_ring_test")

com=intersect(rownames(xAGAP),rownames(xDIA))
xAGAPr=xAGAP[which(rownames(xAGAP) %in% com),51:750]
xDIAr=xDIA[which(rownames(xDIA) %in% com),seq(1,ncol(xDIA),2)]
xAGAPrm=aggregate(xAGAPr,by=list(rownames(xAGAPr)), mean)
rownames(xAGAPrm)=xAGAPrm[,1]
xAGAPrm=as.matrix(xAGAPrm[,-1])

idx <- match(rownames(xDIAr), rownames(xAGAPrm))
xAGAPrm=xAGAPrm[idx,]
xAGAPrmp=snv(log(1/xAGAPrm[,130:ncol(xAGAPrm)]))
xDIArp=snv(xDIAr[,130:ncol(xDIAr)])

x=rbind(xAGAPrm,xDIAr)
col=rep("red",nrow(x))
col[1:nrow(xAGAPrm)]="black"
plotsp(x,col=col, xlab="wavelength (nm)",ylab="Absorbance")
legend(x="topright",legend=c("DIASCOPE","AGAP"), col = c("red","black"),lty=1)
title("MALANIRS - Spectra from ring-test (SNV) Common range")

# Save to ChemFlow
xAGAPcf=xAGAPrm
colnames(xAGAPcf)=paste0("V",1:ncol(xAGAPcf))
xDIAcf=xDIAr
colnames(xDIAcf)=paste0("V",1:ncol(xDIAcf))

x=rbind(xAGAPcf, xAGAPcf[rep(1:nrow(xAGAPcf), each = 100), ])
write.table(x, file = "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/xAGAPcf.txt", row.names = 1:nrow(x), sep="\t", quote = F)
x=rbind(xDIAcf, xDIAcf[rep(1:nrow(xDIAcf), each = 100), ])
write.table(x,file="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/xDIAcf.txt", row.names = 1:nrow(x), sep="\t", quote = F)






## Standardisation PDS
ical=1:10
xAGAPcal=xAGAPrm[ical,]
xDIAcal=xDIAr[ical,]
xAGAPval=xAGAPrm[-ical,]
xDIAval=xDIAr[-ical,]
w=10
colw=w+5:(nccol(xAGAPcal)-2*w)

# Presentation des différences
plotspMALA(xAGAPcal,xDIAcal,unit="Absorbance/Reflectance","DIASCOPE","AGAP")
plotspMALA(log(1/xAGAPcal), xDIAcal,unit="Absorbance","DIASCOPE","AGAP", "To abs")
plotspMALA(snv(log(1/xAGAPcal)), snv(xDIAcal),unit="Absorbance","DIASCOPE","AGAP", "To abs + SNV")
plotspMALA(snv(log(1/xAGAPcal[,151:700])), snv(xDIAcal[,151:700]),unit="Absorbance","DIASCOPE","AGAP", "To abs + SNV")

# PDS sur brut
mPDS = PDS(xAGAPcal, xDIAcal,  w, 1)
xDIApred<-as.matrix(xDIAval)%*%as.matrix(mPDS$P)
xDIApred<-sweep(xDIApred, 2, as.numeric(t(mPDS$Intercept)), "+")
plotspMALA(xAGAPval[,colw],xDIApred[,colw],unit="Reflectance","DIASCOPE","AGAP", "DIASCOPE to AGAP")
plotspMALA(snv(xAGAPval)[,colw],snv(xDIApred)[,colw],unit="Reflectance","DIASCOPE","AGAP", "DIASCOPE to AGAP + SNV")
plotsp(snv(xAGAPval)[,colw]-snv(xDIApred)[,colw])

# xDIAcalpred<-as.matrix(xDIAcal)%*%as.matrix(mPDS$P)
# xDIAcalpred<-sweep(xDIAcalpred, 2, as.numeric(t(mPDS$Intercept)), "+")
# plotspMALA(xAGAPcal,xDIAcalpred,unit="Reflectance","DIASCOPE","AGAP", "DIASCOPE to AGAP")


# PDS sur Abs
mPDS = PDS(log(1/xAGAPcal), xDIAcal,  w, 1)
xDIApred<-as.matrix(xDIAval)%*%as.matrix(mPDS$P)
xDIApred<-sweep(xDIApred, 2, as.numeric(t(mPDS$Intercept)), "+")
plotspMALA(log(1/xAGAPval)[,colw],xDIApred[,colw],unit="Absorbance","DIASCOPE","AGAP", "Abs -> DIASCOPE to AGAP")
plotspMALA(snv(log(1/xAGAPval)[,colw]),snv(xDIApred)[,colw],unit="Absorbance","DIASCOPE","AGAP", "Abs -> DIASCOPE to AGAP + SNV")
plotsp(snv(log(1/xAGAPval))[,colw]-snv(xDIApred)[,colw])


# Application de la moyenne des diff
dif=colMeans(log(1/xAGAPcal)-xDIAcal)
difdup=do.call(rbind, replicate(nrow(xDIAval), dif, simplify = FALSE))
xDIApred=xDIAval+difdup
plotspMALA(log(1/xAGAPval), xDIApred,unit="Absorbance","DIASCOPE","AGAP", "STD par difference")




