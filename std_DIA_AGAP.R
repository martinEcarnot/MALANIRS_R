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
ld=list.dirs("/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/MLD24DIASCOPE/")
sp=read_spectra(ld)
sp$ech=substr(sp$names,1,12)
# sp$ech=gsub("-","_",sp$ech)
xAGAP=sp$value
row.names(xAGAP)=sp$ech
colnames(xAGAP)=sp$bands
xAGAP=aggregate(xAGAP,by=list(sp$ech), mean)
rownames(xAGAP)=xAGAP[,1]
xAGAP[,1]=NULL
# xAGAP=adj_NaturaSpec(log(1/as.matrix(xAGAP)))
# xAGAP=adj_NaturaSpec(as.matrix(xAGAP))
xAGAP=as.matrix(xAGAP)

# Select common samples and wavelength range
com=intersect(rownames(xAGAP),rownames(xDIA))
xAGAPr=xAGAP[match(com,rownames(xAGAP)),51:750]
xDIAr=xDIA[match(com,rownames(xDIA)),seq(1,ncol(xDIA),2)]

setdiff(rownames(xAGAPr),rownames(xDIAr))


## Tune PDS with 127 samples from MLD24DIASCOPE

# 1. Chargement des packages requis
library(prospectr)

# 2. Chargement des données
# X1 : Spectres sur le spectromètre Master (référence) 
X1 <- xAGAPr 
# X2 : Spectres sur le spectromètre Slave (à corriger) 
X2 <- 1/exp(xDIAr)
# X2 <- xDIAr

# -----------------------------------------------------------------------------
# STEP 1 : Sélection des échantillons d'apprentissage (Kennard-Stone)
# -----------------------------------------------------------------------------
set.seed(3)
n_calibration <- 80 # Nombre d'échantillons retenus pour calculer le transfert
idx_cal <- sample(1:nrow(X1), n_calibration)  # Indices des échantillons de calibration)
idx_val <- setdiff(1:nrow(X1), idx_cal)  # Indices des échantillons de validation
X1_cal <- snv(X1[idx_cal, ])
X2_cal <- snv(X2[idx_cal, ])
X1_val <- snv(X1[idx_val, ])
X2_val <- snv(X2[idx_val, ])

# -----------------------------------------------------------------------------
# STEP 3 : Calcul du transfert et validation de la fenêtre
# -----------------------------------------------------------------------------
# Test de plusieurs tailles de fenêtres glissantes (ex: 3, 5, 9, 13)
window_sizes <- c(11, 31, 51) # divided by 2
rmsd_results <- numeric(length(window_sizes))

# RMSD Avant correction
rmsd_initial <- sqrt(mean((X1_val - X2_val)^2))
cat(sprintf("RMSD initial (sans correction) : %.6f\n", rmsd_initial))

# Recherche de la meilleure fenêtre glissante
for (w in seq_along(window_sizes)) {
  win <- window_sizes[w]

  # PDS from RNIR (cf. MALANIRS_utils.R)
  mPDS <- PDS(masterSpectra = X1_cal, slaveSpectra = X2_cal, MWsize = floor(win/2), Ncomp = 2)
  X2_val_corrected<-X2_val%*%as.matrix(mPDS$P)
  X2_val_corrected<-sweep(X2_val_corrected, 2, as.numeric(t(mPDS$Intercept)), "+")

  # Calcul du RMSD post-correction
  colout=c(1:floor(win/2),(ncol(X1)-floor(win/2)):ncol(X1))
  
  # rmsd_results[w] <- sqrt(mean((X1_val[,-colout] - X2_val_corrected[,-colout])^2))
  # cat(sprintf("RMSD avec fenêtre PDS = %2d : %.6f\n", win, rmsd_results[w]))
  res = compare_spectra((X1_val[,-colout]),(X2_val_corrected[,-colout]))
  cat(sprintf("RMSE = %.5f | NRMSE = %.5f | Bias = %.5f | R² = %.5f | SAM = %.5f rad\n", res$RMSE, res$NRMSE, res$Bias, res$R2, res$SAM_mean_rad))
  plot(plotspgg(rbind(X1_val[,-colout],X2_val_corrected[,-colout]),c(rep("AGAP",nrow(X1_val)),rep("DIASCOPE_corr",nrow(X2_val))),""))
  # plotsp(X1_val[,-colout] - X2_val_corrected[,-colout])
  
  # rmsd_results[w] <- sqrt(mean((snv(X1_val[,-colout]) - snv(X2_val_corrected[,-colout]))^2))
  # cat(sprintf("RMSD avec fenêtre PDS = %2d : %.6f\n", win, rmsd_results[w]))
  # plot(plotspgg(rbind(snv(X1_val[,-colout]),snv(X2_val_corrected[,-colout])),c(rep("AGAP",nrow(X1_val)),rep("DIASCOPE_corr",nrow(X2_val))),""))
  # # plotsp(snv(X1_val[,-colout]) - snv(X2_val_corrected[,-colout]))
}

stop()
# -----------------------------------------------------------------------------
# STEP 4 : Application finale avec la meilleure fenêtre (remplacer X2_val par des nouveaux spectres SNV)
# -----------------------------------------------------------------------------
# Calcul de la matrice F optimale et correction de l'ensemble de la base X2
bestw=15
colout=c(1:bestw,(ncol(X1_cal)-bestw):ncol(X1_cal))
mPDS <- PDS(masterSpectra = X1_cal, slaveSpectra = X2_cal, MWsize = 25, Ncomp = 2)
X2_val_corrected<-X2_val%*%as.matrix(mPDS$P)
X2_val_corrected<-sweep(X2_val_corrected, 2, as.numeric(t(mPDS$Intercept)), "+")[,-colout]



# Conclusions
#
# Utilisation de la PDS
# lot de grains de maïs passé sur NaturaSpec (en reflectance sur table tournante) et sur la chaine post-recolte de Diascope
# 80 lots de calibration et 47 lots de validation
# Maitre: NaturSpec, esclave: Diascope
# Critère d'adéquation entre NaturSpec et Diascope standarisé =  sqrt(mean((X1_val[,-colout] - X2_val_corrected[,-colout])^2))
# Mais ce critère est à véréfier, car est peu sensible aux spectres très bruités
# On utilisera plutot NRMSE, Bias ,R² et  SAM
#
# L'effet des prétraitements (avant ou après la correction PDS) est testés dans la fonction pds_pre de MALANIRS_utils.R
# Globalement, il est préférable de prétraiter avant la correction PDS.
#  (attention aux effets de bords d'un SG qui explosent les diff après pré)
# 
# Pour MALANIRS, dont le but est de faire de la selection phénomique, il vaut mieux obtenir des spectres bruts (voire SNV), 
# donc mieux vaut ne pas prétraiter (ou SNV)
#
# !!!! Le biais est toujours nul qd on compare 2 spectres prétraités SNV (après PDS)
#
# Effet du nombre de composantes de la PDS: sans prétraiter, augmenter ncomp améliore le RMSE, R² et SAM, mais générère dans specrtes très bruités
# Mais si on fait SNV après PDS : résultats différents: c'est meilleur avec 1 seule composante.
#
# C'est mieux si SNV avant PDS
# passer xDIA en reflectance n'améliore pas
# Utiliser adj_Naturaspec n'améliore pas.
#
#




