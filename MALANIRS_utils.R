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


