# Create MALANIRS Database by associating lists of samples and spectra files

library(spectrolab)
library(nirsextra)
library(readxl)
source("MALANIRS_utils.R")

# GQE
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/Liste semences MRS_complète_GQE.xlsx"
smpGQE=read_excel(f)
spGQE=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/GQE/")
spGQE$ech=substr(row.names(spGQE$x),1,nchar(row.names(spGQE$x))-6)
spGQE$code=sub("^[^ ]+ ", "", spGQE$ech)  # remove what is after 1st " "
spGQE$code=gsub(" ","_",spGQE$code)
spGQE$code=gsub("EVA_ZM","EVA_Zm",spGQE$code)
spGQE$code=gsub("_REP[0-9]+$", "",spGQE$code)

check_sp(smpGQE$`MRS / EVA Code`,spGQE$code)

rownames(spGQE$x)=paste0("FRA_",rownames(spGQE$x))
write.csv(spGQE$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_GQE2025.csv", row.names = TRUE, quote = FALSE)

#CREA
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/List_SMTA_28.10.2025_MALANIRS_CREA_INRAE.xlsx"
smpCREA=read_excel(f)
spCREA=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/CREA/")
spCREA$ech=row.names(spCREA$x)
spCREA$code=gsub(" ","",spCREA$ech)
spCREA$code=gsub("_","",spCREA$code)

check_sp(smpCREA$`CREA-Landrace Genebank short CODE`,spCREA$code)

rownames(spCREA$x)=paste0("ITA_",rownames(spCREA$x))
write.csv(spCREA$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CREA2025.csv", row.names = TRUE, quote = FALSE)

# Serbia
spSerb=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/Serbie/")
rownames(spSerb$x)=paste0("SRB_",rownames(spSerb$x))
rownames(spSerb$x)[substr(rownames(spSerb$x),5,7)=="ATA"]=gsub("SRB","TUR",rownames(spSerb$x)[substr(rownames(spSerb$x),5,7)=="ATA"])
write.csv(spSerb$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_Serbia2025.csv", row.names = TRUE, quote = FALSE)


# Portugal
spPor=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/Portugal/")
rownames(spPor$x)=paste0("PRT_",rownames(spPor$x))
write.csv(spPor$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_Portugal2025.csv", row.names = TRUE, quote = FALSE)

# Roumania
spRom=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/Roumanie/")
rownames(spRom$x)=paste0("ROU_",rownames(spRom$x))
write.csv(spRom$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_Romania2025.csv", row.names = TRUE, quote = FALSE)


# CRB_Gamet
spCRB=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/CRB_Gamet/")
rownames(spCRB$x)=paste0("FRA_",rownames(spCRB$x))
write.csv(spCRB$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025.csv", row.names = TRUE, quote = FALSE)
