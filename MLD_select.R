# MineLanDiv sample from DIASCOPE 2024
library(nirsextra)
library(rchemo)
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/ech_DIASCOPE/mld24dia.txt"
sp0=read.table(f,skip = 11,header =TRUE,sep=",")
sp=sp2df(as.matrix(sp0[,-c(1:6,ncol(sp0),ncol(sp0)-1)]))
colnames(sp$x)=as.numeric(substr(colnames(sp$x),2,nchar(colnames(sp$x))))
sp$ech=sp0$Sample.Number

k=300
ksbrut= sampks(sp$x, k)$train
kssnv= sampks(snv(sp$x), k)$train
# fm <- pcaeigen(sp$x,nlv=6)
# cmp=c(5,6)
# plot(fm$T[,cmp])
# points(fm$T[s,cmp], pch = 19, col = 2, cex = 1.5)

write.csv(sort(sp$ech[kssnv]), "MLD24DIA_list_ks.csv", row.names = FALSE)







