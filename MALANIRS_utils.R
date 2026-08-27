#Piecewise Direct Standardization (PDS) algorithm:

#INPUT:   masterSpectra = Spectra acquired with the master instrument (matrix).
#         slaveSpectra = Spectra acquired with the slave instrument (matrix).
#         MWsize = Half size of the moving window (integer).
#         Ncomp = Number of latent variables used in the PLS model (integer).
#         wavelength = wavelength (numeric vector).

#OUTPUT:  P = the PDS transfer matrix.

# From: https://guifh.github.io/RNIR/PDS.html

PDS<-function(masterSpectra, slaveSpectra, MWsize, Ncomp, wavelength){
  
  require(pls)
  
  #Loop Initialization:
  i<-MWsize
  k<-i-1
  #Creation of an empty P matrix:
  P<-matrix(0,nrow=ncol(masterSpectra),ncol=ncol(masterSpectra)-(2*i)+2)
  InterceptReg<-c()
  
  while(i<=(ncol(masterSpectra)-k)){
    
    #PLS regression:
    fit<- plsr(masterSpectra[,i] ~ as.matrix(slaveSpectra[,(i-k):(i+k)]),
               ncomp=Ncomp, scale=F, method="oscorespls")
    
    #Extraction of the regression coefficients:
    coefReg<-as.numeric(coef(fit, ncomp=Ncomp, intercept = TRUE))
    InterceptReg<-c(InterceptReg,coefReg[1])
    coefReg<-coefReg[2:length(coefReg)]
    
    #Add coefficients to the transfer matrix:
    P[(i-k):(i+k),i-k]<-t(coefReg)
    
    rm(coefReg,fit)
    i<-i+1
    
    #Diplay progression:
    cat("\r",paste(round(i/ncol(masterSpectra)*100)," %",sep=""))}
  
  P<-data.frame(matrix(0,nrow=ncol(masterSpectra),ncol=k), P,
                matrix(0,nrow=ncol(masterSpectra),ncol=k))
  InterceptReg<-c(rep(0,k),InterceptReg,rep(0,k)) 
  
  Output<-list(P = P , Intercept = InterceptReg)
  
  return(Output)}



# Function to draw spectra from 2 groups (with plotsp)
plotspMALA <-function(x1,x2,unit,n1,n2,tit=""){
  x=rbind(x1,x2)
  col=rep("red",nrow(x))
  col[1:nrow(x1)]="black"
  plotsp(x,col=col, xlab="wavelength (nm)",ylab=unit, lwd=3)
  legend(x="topright",legend=c(n1,n2), col = c("red","black"),lty=1)
  title(tit)
}



# Function to draw spectra from 2 groups (with ggplot)
plotspgg <- function(x, class = NULL, title = NULL, ribbon = FALSE) {
  library(tidyverse)
  
  # Passage en long
  df_long <- as.data.frame(x) %>%
    rownames_to_column(var = "spectrum_id") %>%
    pivot_longer(
      cols = -spectrum_id,
      names_to = "wavelength",
      values_to = "intensity"
    ) %>%
    mutate(wavelength = as.numeric(wavelength))
  
  # Statistiques globales (moyenne ± sd)
  df_summary <- df_long %>%
    group_by(wavelength) %>%
    summarise(mean = mean(intensity),
              sd   = sd(intensity), .groups = "drop")
  
  
  if (!is.null(class)) {
    if (length(class) != nrow(x)) stop("class doit avoir la même longueur que nrow(x)")
    
    # Ajouter la classe dupliquée
    df_long$class <- rep(class, each = ncol(x)) |> factor()
    
    p <- ggplot()
    
    # Ajouter ruban si demandé
    if (ribbon) {
      p <- p + geom_ribbon(
        data = df_summary,
        aes(x = wavelength, ymin = mean - sd, ymax = mean + sd),
        fill = "grey70", alpha = 0.3
      )
    }
    
    p <- p +
      geom_line(
        data = df_long,
        aes(x = wavelength, y = intensity, group = spectrum_id, color = class),
        alpha = 0.8, size=1
      ) +
      labs(x = "Wavelength", y = "Intensity", color = "Class", title = title) +
      theme_minimal() +
      theme(
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.title   = element_text(size = 15, face = "bold", hjust = 0.5)
      )
    
  } else {
    p <- ggplot()
    
    p <- p +
      geom_line(
        data = df_long,
        aes(x = wavelength, y = intensity, group = spectrum_id),
        color = "grey50", alpha = 0.8, size=1
      ) +
      labs(x = "Wavelength", y = "Intensity", title = title) +
      theme_minimal()
    
    if (ribbon) {
      p <- p + geom_ribbon(
        data = df_summary,
        aes(x = wavelength, ymin = mean - sd, ymax = mean + sd),
        fill = "grey30", alpha = 0.3
      )
    }
    

  }
  
  
  return(p)
}






# Fonction de  Comparaison entre smp$GQE$`MRS / EVA Code` et spGQE$code ---

check_sp=function(liste_ref,liste_sp) {

  # ---- Comparaison principale ----
  
  # Nombre d'occurrences de chaque code de la liste de référence dans la liste sp
  nb_occurrences <- sapply(liste_ref, function(x) sum(liste_sp == x, na.rm = TRUE))
  
  # Tableau récapitulatif : pour chaque code de référence
  tableau_comparaison <- data.frame(
    Code = liste_ref,
    Occurrences_dans_spGQE = nb_occurrences,
    Présent_dans_spGQE = ifelse(nb_occurrences > 0, "✅ Oui", "❌ Non"),
    stringsAsFactors = FALSE
  )
  
  # ---- Résumé global ----
  nb_total   <- length(liste_ref)
  nb_present <- sum(nb_occurrences > 0)
  nb_absent  <- sum(nb_occurrences == 0)
  
  resume <- data.frame(
    Total = nb_total,
    Présents = nb_present,
    Absents = nb_absent,
    `% Présents` = round(100 * nb_present / nb_total, 1)
  )
  
  # ---- Comparaison inverse ----
  # Codes présents dans spGQE mais absents dans la liste de référence
  codes_sp_non_trouves <- liste_sp[!liste_sp %in% liste_ref]
  codes_sp_non_trouves <- unique(codes_sp_non_trouves)
  
  # ---- Affichage ----
  cat("=== Résumé global ===\n")
  print(resume)
  
  cat("\n=== Codes de référence absents de spGQE ===\n")
  print(unique(liste_ref[nb_occurrences == 0]))
  
  cat("\n=== Codes de spGQE absents de la liste de référence ===\n")
  print(codes_sp_non_trouves)
  
  # ---- Export CSV ----
  write.csv(tableau_comparaison, "comparaison_codes.csv", row.names = FALSE)
  write.csv(data.frame(Code_sp_absent = codes_sp_non_trouves),
            "codes_sp_absents.csv", row.names = FALSE)
  
  cat("\n✅ Fichiers enregistrés :\n")
  cat(" - comparaison_codes.csv (résumé principal)\n")
  cat(" - codes_sp_absents.csv (codes de spGQE absents dans la liste de référence)\n")
}



# PDS: Test si c'est mieux de prétraiter avant la correction 
# Au final,c'est mieux de prétraiter avant
# (attention aux effets de bords
# d'un SG qui explosent les diff après pré)

pds_pre <- function(X1_cal,X2_cal,X1_val,X2_val) {
  
  # Test de plusieurs tailles de fenêtres glissantes (ex: 3, 5, 9, 13)
  window_sizes <-  c(5,11, 15, 51)
  ncomp=1
  rmsd_results <- numeric(length(window_sizes))
  source("MALANIRS_list_pre.R")
  
  # RMSD Avant correction et sans pretraitement
  rmsd_initial <- sqrt(mean((X1_val - X2_val)^2))
  cat(sprintf("RMSD initial (sans PDS) : %.6f\n", rmsd_initial))
  
  for (i in seq_along(list_pre)) {
    X1_calp <- pre(X1_cal, list_pre[[i]])
    X2_calp <- pre(X2_cal, list_pre[[i]])
    X1_valp <- pre(X1_val, list_pre[[i]])
    X2_valp <- pre(X2_val, list_pre[[i]])
  
    # RMSD Avant correction
    rmsd_initial <- sqrt(mean((X1_valp - X2_valp)^2))
    cat(sprintf("RMSD initial (sans PDS) : %.6f\n", rmsd_initial))
    
    # Recherche de la meilleure fenêtre glissante
    for (w in seq_along(window_sizes)) {
      win <- window_sizes[w]
      coloutp=c(1:floor(win),(ncol(X1_calp)-floor(win/2)):ncol(X1_calp))
      colout=c(1:floor(win),(ncol(X1)-floor(win/2)):ncol(X1))
      
      ## Prétraitement avant
      mPDS <- PDS(masterSpectra = X1_calp, slaveSpectra = X2_calp, MWsize = floor(win/2), Ncomp = ncomp)
      X2_valp_corrected<-X2_valp%*%as.matrix(mPDS$P)
      X2_valp_corrected<-sweep(X2_valp_corrected, 2, as.numeric(t(mPDS$Intercept)), "+")
      # Calcul du RMSD post-correction
      rmsd_results[w] <- sqrt(mean((X1_valp[,-coloutp] - X2_valp_corrected[,-coloutp])^2))
      cat(sprintf("RMSD avec fenêtre PDS (Pré avant) = %2d : %.6f\n", win, rmsd_results[w]))
      plot(plotspgg(rbind(X1_valp[,-coloutp] , X2_valp_corrected[,-coloutp]),c(rep("AGAP",nrow(X1_val)),rep("DIASCOPE_corr",nrow(X2_val))),"Pré avant"))
      
      ## Prétraitement après
      mPDS <- PDS(masterSpectra = X1_cal, slaveSpectra = X2_cal, MWsize = floor(win/2), Ncomp = ncomp)
      X2_val_corrected<-X2_val%*%as.matrix(mPDS$P)
      X2_val_corrected<-sweep(X2_val_corrected, 2, as.numeric(t(mPDS$Intercept)), "+")
      
      # Calcul du RMSD post-correction
      rmsd_results[w] <- sqrt(mean((pre(X1_val[,-colout], list_pre[[i]]) - pre(X2_val_corrected[,-colout], list_pre[[i]]))^2))
      cat(sprintf("RMSD avec fenêtre PDS  (Pré après) = %2d : %.6f\n", win, rmsd_results[w]))
      plot(plotspgg(rbind(pre(X1_val[,-colout], list_pre[[i]]) , pre(X2_val_corrected[,-colout], list_pre[[i]])),c(rep("AGAP",nrow(X1_val)),rep("DIASCOPE_corr",nrow(X2_val))),"Pré après"))
      # browser()
    }
  }
}


# Fonction to compare two spectra matrices (master and slave_std) and compute various metrics including RMSE, Bias, NRMSE, R², and SAM.
compare_spectra <- function(master, slave_std) {

  N <- nrow(master)
  P <- ncol(master)
  
  # ============================================================
  # RMSE = sqrt(1/(N*P) * somme((Master - Slave_PDS)^2))
  # ============================================================
  rmse <- sqrt(
    mean((master - slave_std)^2, na.rm = TRUE)
  )
  
  # ============================================================
  # NRMSE : Normalisation par l'écart-type des valeurs du maître
  # ============================================================
  nrmse <- rmse / sd(as.numeric(master), na.rm = TRUE)

  # ============================================================
  # Biais
  # ============================================================
    bias <- mean(
    slave_std - master,
    na.rm = TRUE
  )
  
  # ============================================================
  # R²
  # ============================================================
  valid <- is.finite(master) & is.finite(slave_std)
  r <- cor(
    master[valid],
    slave_std[valid]
  )
  r2 <- r^2
  
  # ============================================================
  # SAM : spectre par spectre
  # ============================================================
  #
  # SAM_i = acos(
  #   sum(Master_i * Slave_i) /
  #   (||Master_i|| * ||Slave_i||)
  # )
  #
  # Puis moyenne des SAM_i
  #
  
  sam_individual <- rep(NA_real_, N)
  
  for (i in seq_len(N)) {
    
    x <- master[i, ]
    y <- slave_std[i, ]

    # Vérification des normes
    norm_x <- sqrt(sum(x^2))
    norm_y <- sqrt(sum(y^2))

    if (length(x) == 0 || norm_x == 0 || norm_y == 0) {
      sam_individual[i] <- NA_real_
      next
    }
    
    # Cosinus de l'angle
    cos_theta <- sum(x * y) / (norm_x * norm_y)
    # Protection contre les erreurs numériques
    cos_theta <- max(-1, min(1, cos_theta))
    
    # SAM en radians
    sam_individual[i] <- acos(cos_theta)
  }
  
  # SAM moyen
  sam_mean_rad <- mean(
    sam_individual,
    na.rm = TRUE
  )

  # ============================================================
  # Résultats
  # ============================================================
  return(list(
    RMSE = rmse,
    NRMSE = nrmse,
    Bias = bias,
    R2 = r2,
    SAM_mean_rad = sam_mean_rad
  ))
}