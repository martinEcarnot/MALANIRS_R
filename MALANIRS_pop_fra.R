# MALANIRS_pop_fra.R
# Association des spectres NIRS (CRB Gamet) avec les données biochimiques
# des populations françaises via les codes PPS
# Explore les traits biochimuqes du fichier table_mean_adjusted_NIRS_biochimie_res_rowcol_GxE.csv

library(readxl)

base <- "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/"

# -------------------------------------------------------------------------
# 1. Lecture des données
# -------------------------------------------------------------------------

# Spectres NIRS CRB Gamet (rownames = "FRAxxxxxxx Ryyyy_0000N")
nirs <- read.csv(paste0(base, "NIRS_CRBGamet2025.csv"),
                 row.names = 1, check.names = FALSE)

# Données biochimiques (séparateur ";" décimale ",", codes PPS)
bio <- read.csv(paste0(base, "pop_francaise/table_mean_adjusted_NIRS_biochimie_res_rowcol_GxE.csv"),
                sep = ";", dec = ",", row.names = 1)

# Table de correspondance Code_FRA <-> Code_PPS
tPop <- read_excel(paste0(base, "pop_francaise/PopFra_262pop_notification_sitecollecte_250220.xlsx"))

# -------------------------------------------------------------------------
# 2. Moyenne des scans NIRS par accession (FRA code)
# -------------------------------------------------------------------------

# Extraire le code FRA depuis les rownames ("FRAxxxxxxx Ryyyy_0000N" -> "FRAxxxxxxx")
nirs_fra <- sub(" .*", "", sub("_[0-9]+$", "", rownames(nirs)))

# Moyenne par accession
nirs_mean <- aggregate(nirs, by = list(FRA = nirs_fra), FUN = mean)
rownames(nirs_mean) <- nirs_mean$FRA
nirs_mean$FRA <- NULL

cat("Accessions NIRS après moyenne:", nrow(nirs_mean), "\n")

# -------------------------------------------------------------------------
# 3. Association Code_FRA -> Code_PPS via tPop
# -------------------------------------------------------------------------

lookup_pps <- setNames(tPop$Code_PPS, tPop$Code_FRA)
nirs_pps   <- lookup_pps[rownames(nirs_mean)]

cat("Accessions avec Code_PPS trouvé:", sum(!is.na(nirs_pps)), "\n")
cat("Accessions sans correspondance PPS:", sum(is.na(nirs_pps)), "\n")

# Garder uniquement les accessions avec un code PPS
nirs_matched <- nirs_mean[!is.na(nirs_pps), ]
rownames(nirs_matched) <- nirs_pps[!is.na(nirs_pps)]

# -------------------------------------------------------------------------
# 4. Jointure avec les données biochimiques
# -------------------------------------------------------------------------

common_pps <- intersect(rownames(nirs_matched), bio$genotype)
cat("Accessions en commun NIRS x biochimie:", length(common_pps), "\n")

bio_matched  <- bio[bio$genotype %in% common_pps, ]
rownames(bio_matched) <- bio_matched$genotype
bio_matched$genotype  <- NULL

nirs_final <- nirs_matched[common_pps, ]
bio_final  <- bio_matched[common_pps, ]

# Vérification alignement
stopifnot(all(rownames(nirs_final) == rownames(bio_final)))

# -------------------------------------------------------------------------
# 5. Export
# -------------------------------------------------------------------------

write.csv(nirs_final, paste0(base, "pop_francaise/NIRS_popFra_matched.csv"),
          row.names = TRUE, quote = FALSE)
write.csv(bio_final, paste0(base, "pop_francaise/bio_popFra_matched.csv"),
          row.names = TRUE, quote = FALSE)

cat("Fichiers exportés:\n")
cat(" -", paste0(base, "pop_francaise/NIRS_popFra_matched.csv"), "\n")
cat(" -", paste0(base, "pop_francaise/bio_popFra_matched.csv"), "\n")


library(dplyr)

# Tous les traits avec _SMH et _MAU
all_traits <- str_remove(smh_cols, "_SMH$")

df_all <- bind_rows(lapply(all_traits, function(trait) {
  bio_final[, c(paste0(trait, "_SMH"), paste0(trait, "_MAU"))] |>
    setNames(c("SMH", "MAU")) |>
    mutate(trait = trait)
}))

r2_all <- df_all |>
  group_by(trait) |>
  summarise(
    r2  = cor(SMH, MAU, use = "complete.obs")^2,
    SMH = max(SMH, na.rm = TRUE),
    MAU = min(MAU, na.rm = TRUE)
  )

ggplot(df_all, aes(x = SMH, y = MAU)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.7) +
  geom_text(data = r2_all, aes(label = sprintf("R²=%.2f", r2)),
            hjust = 1.05, vjust = 0, size = 2.5) +
  facet_wrap(~trait, scales = "free", ncol = 6) +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 7))



# -------------------------------------------------------------------------
# 6. Single Kernel (CRB Gamet - Brimrose)
# -------------------------------------------------------------------------

# Chargement
nirs_sk <- read.csv(paste0(base, "NIRS_CRBGamet2025_SingleKernel.csv"),
                    row.names = 1, check.names = FALSE)

cat("Single kernel chargé:", nrow(nirs_sk), "grains,", ncol(nirs_sk), "longueurs d'onde\n")

# Tirage aléatoire de 10% des grains
set.seed(4721)
idx_sample <- sample(nrow(nirs_sk), size = round(nrow(nirs_sk) * 1))
nirs_sk_s  <- nirs_sk[idx_sample, ]

cat("Après tirage 10%:", nrow(nirs_sk_s), "grains\n")

# Extraction du code FRA depuis les rownames ("FRAxxxxxxx Ryyyy_NNN" -> "FRAxxxxxxx")
sk_fra <- sub(" .*", "", rownames(nirs_sk_s))

# FRA -> PPS via lookup_pps (défini section 3)
sk_pps <- lookup_pps[sk_fra]

cat("Grains avec PPS trouvé:", sum(!is.na(sk_pps)), "\n")
cat("Grains sans correspondance PPS:", sum(is.na(sk_pps)), "\n")

# Garder uniquement les grains avec un code PPS présent dans bio_final
keep <- !is.na(sk_pps) & sk_pps %in% rownames(bio_final)
nirs_sk_matched <- nirs_sk_s[keep, ]
sk_pps_matched  <- sk_pps[keep]

# Association : chaque grain reçoit les valeurs bio de son accession (PPS)
bio_sk <- bio_final[sk_pps_matched, ]
rownames(bio_sk) <- rownames(nirs_sk_matched)   # renommer avec l'id grain

cat("Grains associés à une ligne bio:", nrow(nirs_sk_matched), "\n")
stopifnot(nrow(nirs_sk_matched) == nrow(bio_sk))

write.csv(nirs_sk_matched,
          "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025_SingleKernel_red.csv",
          row.names = TRUE, quote = FALSE)

write.csv(bio_sk,
          "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025_SingleKernel_red_bioch.csv",
          row.names = TRUE, quote = FALSE)

# Moyenne par accession (PPS) des spectres single kernel
nirs_sk_mean <- aggregate(nirs_sk_matched,
                          by = list(PPS = sk_pps_matched),
                          FUN = mean)
rownames(nirs_sk_mean) <- nirs_sk_mean$PPS
nirs_sk_mean$PPS <- NULL

cat("Accessions SK (moyenne):", nrow(nirs_sk_mean), "\n")

# Bio correspondante (une ligne par accession)
bio_sk_mean <- bio_final[rownames(nirs_sk_mean), ]

stopifnot(all(rownames(nirs_sk_mean) == rownames(bio_sk_mean)))

write.csv(nirs_sk_mean,
          "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025_SingleKernel_mean.csv",
          row.names = FALSE, quote = FALSE)

write.csv(bio_sk_mean,
          "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025_SingleKernel_mean_bioch.csv",
          row.names = FALSE, quote = FALSE)

cat("Fichiers moyennes SK exportés:\n")
cat(" -", "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025_SingleKernel_mean.csv", "\n")
cat(" -", "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais/smpl_2025/NIRS_CRBGamet2025_SingleKernel_mean_bioch.csv", "\n")
