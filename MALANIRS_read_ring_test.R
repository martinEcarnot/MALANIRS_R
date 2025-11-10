library(spectrolab)
library(rchemo)
library(nirsextra)
source("MALANIRS_utils.R")

# # Bash command to list spectra in directories
# cd /home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/
# find . -type f -name "*.sed" -exec basename {} \; | sed 's/..........$//' | sort -u > list_done_Naturaspec.txt
# grep "MLD" list_done_Naturaspec.txt | awk '{num=substr($0, length($0)-2); print num, $0}' | sort -n | cut -d' ' -f2- > list_done_Naturaspec_MLD.txt
# grep "MRS" list_done_Naturaspec.txt | awk '{num=substr($0, length($0)-2); print num, $0}' | sort -n | cut -d' ' -f2- > list_done_Naturaspec_MRS.txt
# grep -v -e "MLD" -e "MRS" list_done_Naturaspec.txt > list_done_Naturaspec_autres.txt

## Read and organize spectra
xAGAP=NULL
ld=list.dirs("/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/2025_Jun_03-mais_MLD")
sp=read_spectra(ld)
sp$ech=substr(sp$names,1,nchar(sp$names)-10)
sp$ech=gsub("-","_",sp$ech)
xAGAP1=sp$value
row.names(xAGAP1)=sp$ech
colnames(xAGAP1)=sp$bands
xAGAP=xAGAP1 #rbind(xAGAP,xAGAP1)
  
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/mala_DIA.txt"
spD=read.table(f,skip = 11,header =TRUE,sep=",")
spD$ech=substr(spD$Sample.Number,6,12)
xDIA=as.matrix(spD[,-c(1:6,ncol(spD),ncol(spD)-1)])
row.names(xDIA)=make.unique(spD$ech)
colnames(xDIA)=substr(colnames(spD[,-c(1:6,ncol(spD),ncol(spD)-1)]),2,8)

# save(xAGAP,xDIA, file = "malanirs_ring_test")


# Make them fit to the AGAP wavelength range
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

x=rbind(log(1/xAGAPrm),xDIAr)
col=rep("red",nrow(x))
col[1:nrow(xAGAPrm)]="black"
plotsp(x,col=col, xlab="wavelength (nm)",ylab="Absorbance")
legend(x="topright",legend=c("DIASCOPE","AGAP"), col = c("red","black"),lty=1)
title("MALANIRS - Spectra from ring-test (SNV) Common range")

# Save to ChemFlow
# xAGAPcf=xAGAPrm
# colnames(xAGAPcf)=paste0("V",1:ncol(xAGAPcf))
# xDIAcf=xDIAr
# colnames(xDIAcf)=paste0("V",1:ncol(xDIAcf))

## Write data set
write.table(xAGAP, file = "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/sp_AGAP.csv",row.names = F,sep=";", quote = F)
write.table(cbind(rownames(xAGAP),rep("Reflectance_0_1",nrow(xAGAP))), file = "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/sp_AGAP_meta.csv",row.names = F,col.names = c("sample","unit"),sep=";", quote = F)
write.table(xDIA, file = "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/sp_DIA.csv",row.names = F,sep=";", quote = F)
write.table(cbind(rownames(xDIA),rep("Reflectance_0_1",nrow(xDIA))), file = "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/sp_DIA_meta.csv",row.names = F,col.names = c("sample","unit"),sep=";", quote = F)

# x=rbind(xAGAPcf, xAGAPcf[rep(1:nrow(xAGAPcf), each = 100), ])
# write.table(x, file = "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/xAGAPcf.txt", row.names = 1:nrow(x), sep="\t", quote = F)
# x=rbind(xDIAcf, xDIAcf[rep(1:nrow(xDIAcf), each = 100), ])
# write.table(x,file="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ring_test/xDIAcf.txt", row.names = 1:nrow(x), sep="\t", quote = F)






## Standardisation PDS
ical=1:10
xAGAPcal=xAGAPrm[ical,]
xDIAcal=xDIAr[ical,]
xAGAPval=xAGAPrm[-ical,]
xDIAval=xDIAr[-ical,]
w=51
colw=w+5:(ncol(xAGAPcal)-2*w)

# Plot differences between 2 spectro.
# plotspMALA(log(1/xAGAPcal), xDIAcal,unit="Absorbance","DIASCOPE","AGAP", "To abs")
# plotspMALA(snv(log(1/xAGAPcal)), snv(xDIAcal),unit="Absorbance","DIASCOPE","AGAP", "To abs + SNV")
# plotspMALA(snv(log(1/xAGAPcal[,151:700])), snv(xDIAcal[,151:700]),unit="Absorbance","DIASCOPE","AGAP", "To abs + SNV")
plotspgg(rbind(log(1/xAGAPcal),xDIAcal),c(rep("AGAP",nrow(xAGAPcal)),rep("DIASCOPE",nrow(xDIAcal))),"Raw spectra")
plotspgg(rbind(snv(log(1/xAGAPcal)), snv(xDIAcal)),c(rep("AGAP",nrow(xAGAPcal)),rep("DIASCOPE",nrow(xDIAcal))),"Pretreatement SNV")

# PDS on raw data
mPDS = PDS(xAGAPcal, xDIAcal,  w, 1)
xDIApred<-as.matrix(xDIAval)%*%as.matrix(mPDS$P)
xDIApred<-sweep(xDIApred, 2, as.numeric(t(mPDS$Intercept)), "+")

plotspgg(rbind(xAGAPval[,colw],xDIApred[,colw]),c(rep("AGAP",nrow(xAGAPval)),rep("DIASCOPE",nrow(xDIAval))),"DIASCOPE to AGAP Raw")
plotspgg(snv(rbind(xAGAPval[,colw],xDIApred[,colw])),c(rep("AGAP",nrow(xAGAPval)),rep("DIASCOPE",nrow(xDIAval))),"DIASCOPE to AGAP SNV")
plotspgg(rbind(snv(xAGAPval[,colw])-snv(xDIApred[,colw])),title = "DIFF DIASCOPE - AGAP SNV")

# PDS on log(1/R)
mPDS = PDS(log(1/xAGAPcal), xDIAcal,  w, 1)
xDIApred<-as.matrix(xDIAval)%*%as.matrix(mPDS$P)
xDIApred<-sweep(xDIApred, 2, as.numeric(t(mPDS$Intercept)), "+")
plotspgg(rbind(log(1/xAGAPval[,colw]),xDIApred[,colw]),c(rep("AGAP",nrow(xAGAPval)),rep("DIASCOPE",nrow(xDIAval))),"DIASCOPE to AGAP Raw")
plotspgg(snv(rbind(log(1/xAGAPval[,colw]),xDIApred[,colw])),c(rep("AGAP",nrow(xAGAPval)),rep("DIASCOPE",nrow(xDIAval))),"DIASCOPE to AGAP SNV")
plotspgg(rbind(snv(log(1/xAGAPval[,colw]))-snv(xDIApred[,colw])),title = "DIFF DIASCOPE - AGAP SNV")

# Applying difference of average
dif=colMeans(log(1/xAGAPcal)-xDIAcal)
difdup=do.call(rbind, replicate(nrow(xDIAval), dif, simplify = FALSE))
xDIApred=xDIAval+difdup
plotspMALA(log(1/xAGAPval), xDIApred,unit="Absorbance","DIASCOPE","AGAP", "STD par difference")




## All spectra measured in MALANIRS at AGAP -> compute variance and compare to standardisation error

# SD in repetitions from xAGAP 
xAGAPsd=aggregate(snv(xAGAPr[,colw]),by=list(rownames(xAGAPr)), sd)
rownames(xAGAPsd)=xAGAPsd[,1]
xAGAPsd=as.matrix(xAGAPsd[,-1])
plotspgg(xAGAPsd)

# SD between Val from PDS
xPDS=rbind(snv(xAGAPval[,colw]),snv(xDIApred[,colw]))
ad=rep(1:nrow(xAGAPval),2)
xPDSsd=aggregate(xSDPDS,by=list(ad),sd)[,-1]


x=rbind(xPDSsd, xAGAPsd)
class=c(rep("PDS sd",nrow(xAGAPval)),rep("Intra-batch SD",nrow(xAGAPsd)))
plotspgg(x,class,"Comparing SD on SNV spectra")



