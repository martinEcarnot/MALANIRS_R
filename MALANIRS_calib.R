library(rchemo)
library(nirsextra)

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
options(warn=-1)

for (i in 1:length(ag)) {  #
  ag1=ag[i]
  print(ag1)
  pftot=NULL
  ncomp=10
  r2_tt=matrix(nrow=ncomp,ncol=length(list_pre))  # matrix(nrow=ncomp+1,ncol=nseq-6)  #
  fmtt=list()
  fmttn=list()
  
  for (j in 1:length(list_pre)) {  # 1:
    sp$xp=pre(sp$x,list_pre[[j]])
    
    # generate sgm list for Leave-one-out
    segm <- list(rep1 = as.list(1:nrow(sp)))  # segm <- segmkf(n = nrow(datok), K = 5)
    # fm = cvfit(datok$xp, datok[,colnames(datok)==ag1],fun=plsr,segm=segm, ncomp=ncomp)
    fmc = gcvlv(sp$xp, sp[,colnames(sp)==ag1],segm,score = r2, fun = plskern, nlv = 1:ncomp, verb = F)  # !!! pas cor2 avec LOO 
    
    fmtt=append(fmtt, list(fmc)) # fmttn=append(fmttn, list(fm))
    r2_tt[,j] = mser(fmc)$cor2  # r2_tt[,j] = mse(fm, ~ ncomp)$cor2
    pf=t(list_pre[[j]])  # # list(seqlo[j], 2151-seqlo[j+4])
    dim(pf)=c(1,2*ncol(pf))
    pf=paste0(pf, collapse = "_")
    pftot=c(pftot,pf)
  }
  best_pre_lo=which(r2_tt == max(r2_tt[1:ncomp,]), arr.ind = TRUE)[1,]
  # matplot(r2_tt, type = 'l', lty = 1, col = 1:ncol(r2_tt), ylab="R2_Validation_Croisée", xlab="Nombre de Variables Latentes")
  # title(ag1)
  # legend("bottomright", legend = pftot, col = 1:ncol(r2_tt), lty = 1, cex = 0.8, bg = "white")  #cex = 0.6
  fm=fmtt[[best_pre_lo[2]]]  # fmn=fmttn[[best_pre_lo[2]]]
  print(pftot[best_pre_lo[2]], 2)
  print("R2")
  print(round(mser(fm)$cor2,2))  # print(round(mse(fmn, ~ ncomp)$cor2[-1], 2))
  print("SEP")
  print(round(mser(fm)$sep,2))    # print(round(mse(fmn, ~ ncomp)$sep[-1],2))
  cat("\n")
  cat("\n")
  
  fm1=fm$y[fm$y$nlv==best_pre_lo[1],]
  plot(fm1$yp,fm1$yref, xlab="Teneur Prédite par SPIR", ylab="Teneur Mesurée")
  abline(a = 0, b = 1, col = "red", lty = 2)  # a=0, b=1 pour y=x; couleur rouge et ligne en pointillés
  fit <- lm(fm1$yp ~ fm1$yref)
  abline(fit, col = "blue")
  summary_fit <- summary(fit)
  r_squared <- summary_fit$r.squared
  legend("topleft", legend = c("y = x", bquote(Validation_Croisée: ~ R^2 == .(round(r_squared, 2))),bquote(pre :  .(pftot[best_pre_lo[2]])),bquote(ncomp : .(best_pre_lo[[1]])),bquote(n_ech : .(length(fm1$yref)))), col = c("red", "blue", "white", "white", "white"),lty = c(2, 1),bty = "n")
  title(ag1)
}

