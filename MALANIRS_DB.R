# Create MALANIRS Database by associating lists of samples and spectra files

library(spectrolab)
library(nirsextra)
library(readxl)
source("MALANIRS_utils.R")

sourcef="/media/ecarnot/C8DE-2338/MALANIRS/MALANIRS/"
# GQE
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/Liste semences MRS_complète_GQE.xlsx"
smpGQE=read_excel(f)
spGQE=naturaspec2df(paste0(sourcef,"GQE"))
spGQE$ech=substr(row.names(spGQE$x),1,nchar(row.names(spGQE$x))-6)
spGQE$code=sub("^[^ ]+ ", "", spGQE$ech)  # remove what is after 1st " "
spGQE$code=gsub(" ","_",spGQE$code)
spGQE$code=gsub("EVA_ZM","EVA_Zm",spGQE$code)
spGQE$code=gsub("_REP[0-9]+$", "",spGQE$code)

check_sp(smpGQE$`MRS / EVA Code`,spGQE$code)

# rownames(spGQE$x)=paste0("FRA_",rownames(spGQE$x))
rownames(spGQE$x)=rownames(spGQE$x)
write.csv(spGQE$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_GQE2025.csv", row.names = TRUE, quote = FALSE)

#CREA
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/List_SMTA_28.10.2025_MALANIRS_CREA_INRAE.xlsx"
smpCREA=read_excel(f)
spCREA=naturaspec2df(paste0(sourcef,"CREA/"))
spCREA$ech=row.names(spCREA$x)
spCREA$code=gsub(" ","",spCREA$ech)
spCREA$code=gsub("_","",spCREA$code)

check_sp(smpCREA$`CREA-Landrace Genebank short CODE`,spCREA$code)

# rownames(spCREA$x)=paste0("ITA_",rownames(spCREA$x))
rownames(spCREA$x)=rownames(spCREA$x)
write.csv(spCREA$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CREA2025.csv", row.names = TRUE, quote = FALSE)

# Serbia
spSerb=naturaspec2df(paste0(sourcef,"Serbie/"))
rownames(spSerb$x)=paste0("SRB_",rownames(spSerb$x))
rownames(spSerb$x)[substr(rownames(spSerb$x),5,7)=="ATA"]=gsub("SRB","TUR",rownames(spSerb$x)[substr(rownames(spSerb$x),5,7)=="ATA"])
write.csv(spSerb$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_Serbia2025.csv", row.names = TRUE, quote = FALSE)


# Portugal
spPor=naturaspec2df(paste0(sourcef,"Portugal/"))
rownames(spPor$x)=paste0("PRT_",rownames(spPor$x))
write.csv(spPor$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_Portugal2025.csv", row.names = TRUE, quote = FALSE)

# Roumania
spRom=naturaspec2df(paste0(sourcef,"Roumanie/"))
rownames(spRom$x)=paste0("ROU_",rownames(spRom$x))
write.csv(spRom$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_Romania2025.csv", row.names = TRUE, quote = FALSE)


# CRB_Gamet
f_crb="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/pop_francaise/feuille prélèvement ZEA-DST2025-001 Compan (Malanirs).xlsx"
smpCRB=read_excel(f_crb)
smpCRB$lot_last5=as.integer(substr(smpCRB$`N° lots`, nchar(smpCRB$`N° lots`)-4, nchar(smpCRB$`N° lots`)))
lookup_crb=setNames(smpCRB$Accessions, smpCRB$lot_last5)

spCRB=naturaspec2df(paste0(sourcef,"CRB_Gamet/"))
rn=rownames(spCRB$x)
prefix=sub("_.*","",rn)
suffix=sub("^[^_]+","",rn)

# Replace numeric prefixes (= last 5 digits of N° lots) with Accessions
is_num=grepl("^[0-9]+$",prefix)
mapped=lookup_crb[prefix[is_num]]
# Keep original rowname if no match found (e.g. prefixes 63 and 161)
new_rn=rn
new_rn[is_num]=ifelse(is.na(mapped), rn[is_num], paste0(mapped, suffix[is_num]))
rownames(spCRB$x)=new_rn

write.csv(spCRB$x, "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025.csv", row.names = TRUE, quote = FALSE)


