# Create MALANIRS Database by associating lists of samples and spectra files

library(spectrolab)
library(nirsextra)
library(readxl)
source("MALANIRS_utils.R")


# GQE
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/Liste semences MRS_complète.xlsx"
smpGQE=read_excel(f)
spGQE=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/GQE/")
spGQE$ech=row.names(spGQE$x)
spGQE$code=sub("^[^ ]+ ", "", spGQE$ech)  # remove what is after 1st " "
spGQE$code=gsub(" ","_",spGQE$code)
spGQE$code=gsub("EVA_ZM","EVA_Zm",spGQE$code)
spGQE$code=gsub("_REP[0-9]+$", "",spGQE$code)

check_sp(smpGQE$`MRS / EVA Code`,spGQE$code)


#CREA
f="/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/List_SMTA_28.10.2025_MALANIRS_CREA_INRAE.xlsx"
smpCREA=read_excel(f)
spCREA=naturaspec2df("/media/ecarnot/Crucial X9/MALANIRS/CREA/")
spCREA$ech=row.names(spCREA$x)
spCREA$code=gsub(" ","",spCREA$ech)
spCREA$code=gsub("_","",spCREA$code)

check_sp(smpCREA$`CREA-Landrace Genebank short CODE`,spCREA$code)

# Serbia
