library(readxl)

# Compare les mesures biochimiques MBG aux valeurs NIRS, par MLD.
base <- "/home/ecarnot/Documents/INRA/Projets/MalaNIRS_Mais"
f_mbg <- file.path(base, "1_Grain data from MBG_CSIC_MLD_2024_ORIGINAL.xlsx")
f_data <- file.path(
	base,
	"ech_DIASCOPE",
	"MineLandDiv_T1.2_Trials 2024_INRAE Mauguio_20240729.xlsx"
)
f_sortie <- "sorties/comparaison_biochimie_MLD_2024.csv"

mbg <- read_excel(f_mbg)
data <- read_excel(f_data, sheet = "Data")

names(mbg)[names(mbg) == "Code_MLD"] <- "MLD"
names(data)[names(data) == "MineLandDiv Landrace ID"] <- "MLD"

mean_na <- function(x) {
	if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

mbg_summary <- aggregate(
	cbind(Moisture, Protein, Starch) ~ MLD,
	data = mbg,
	FUN = mean_na,
	na.action = na.pass
)
names(mbg_summary)[names(mbg_summary) != "MLD"] <- paste0(
	"mbg_",
	names(mbg_summary)[names(mbg_summary) != "MLD"]
)
mbg_n <- as.data.frame(table(mbg$MLD), stringsAsFactors = FALSE)
names(mbg_n) <- c("MLD", "mbg_n")
mbg_summary <- merge(mbg_summary, mbg_n, by = "MLD", all = TRUE)

data_summary <- aggregate(
	cbind(MLD_OG_NIRS, MLD_PG_NIRS, MLD_SG_NIRS) ~ MLD,
	data = data,
	FUN = mean_na,
	na.action = na.pass
)
names(data_summary)[names(data_summary) != "MLD"] <- paste0(
	"nirs_",
	names(data_summary)[names(data_summary) != "MLD"]
)
data_n <- as.data.frame(table(data$MLD), stringsAsFactors = FALSE)
names(data_n) <- c("MLD", "nirs_n")
data_summary <- merge(data_summary, data_n, by = "MLD", all = TRUE)

comparison <- merge(mbg_summary, data_summary, by = "MLD", all = TRUE)
comparison$delta_moisture <- comparison$nirs_MLD_OG_NIRS - comparison$mbg_Moisture
comparison$delta_protein <- comparison$nirs_MLD_PG_NIRS - comparison$mbg_Protein
comparison$delta_starch <- comparison$nirs_MLD_SG_NIRS - comparison$mbg_Starch
comparison$status <- "matched"
comparison$status[is.na(comparison$mbg_n)] <- "only_in_data"
comparison$status[is.na(comparison$nirs_n)] <- "only_in_mbg"

comparison <- comparison[order(comparison$MLD), ]
write.csv(comparison, f_sortie, row.names = FALSE, na = "")

resume <- function(nom, mbg_col, nirs_col, delta_col) {
	ok <- complete.cases(comparison[, c(mbg_col, nirs_col, delta_col)])
	delta <- comparison[ok, delta_col]
	cat(
		sprintf(
			"%s: n=%d, biais=%.3f, MAE=%.3f, RMSE=%.3f, correlation=%.3f\n",
			nom,
			sum(ok),
			mean(delta),
			mean(abs(delta)),
			sqrt(mean(delta^2)),
			cor(comparison[ok, mbg_col], comparison[ok, nirs_col])
		)
	)
}

cat("CSV écrit : ", f_sortie, "\n", sep = "")
cat("MLD MBG : ", sum(!is.na(unique(mbg$MLD))),
	" | MLD Data : ", sum(!is.na(unique(data$MLD))), "\n", sep = "")
resume("Moisture", "mbg_Moisture", "nirs_MLD_OG_NIRS", "delta_moisture")
resume("Protein", "mbg_Protein", "nirs_MLD_PG_NIRS", "delta_protein")
resume("Starch", "mbg_Starch", "nirs_MLD_SG_NIRS", "delta_starch")

plot(comparison$mbg_Protein,comparison$nirs_MLD_PG_NIRS)