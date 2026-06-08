# vitaspec_SVR.R
# Calibration SPIR par SVR avec recherche des paramètres par validation croisée
# Analogue de vitaspec_preCV mais utilisant e1071::svm() au lieu de PLS
# Grille : C × epsilon × gamma

r2_fn <- function(pred, obs) {
  obs  <- as.numeric(obs)
  pred <- as.numeric(pred)
  1 - sum((obs - pred)^2, na.rm = TRUE) /
      sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
}


vitaspec_SVR <- function(x, y, list_pre,
                         cost_grid    = 10^seq(-1, 3, by = 0.5),  # 9 valeurs : 0.1 → 1000
                         epsilon_grid = c(0.01, 0.1, 0.5),
                         gamma_grid   = NULL,  # NULL => c(0.1, 1, 10) / ncol(xp)
                         kernel       = "radial",
                         K            = 5,
                         titl         = "",
                         plotCV       = TRUE,
                         plotYY       = TRUE,
                         verb         = FALSE) {

  library(e1071)

  # --- Préparation y (suppression NA) ---
  yv    <- as.numeric(y)
  ok    <- !is.na(yv)
  yv_ok <- yv[ok]
  n_ok  <- sum(ok)

  # Folds fixes pour toutes les itérations (comparaisons équitables)
  segm_folds <- segmkf(n = n_ok, K = K)[[1]]

  npre <- length(list_pre)
  nc   <- length(cost_grid)
  ne   <- length(epsilon_grid)

  r2_by_cost  <- matrix(NA_real_, nrow = nc, ncol = npre)
  r2_best     <- numeric(npre)
  sep_best    <- numeric(npre)
  params_best <- vector("list", npre)
  pftot       <- NULL

  for (j in seq_len(npre)) {
    if (verb) cat(paste0("Prétraitement ", j, "/", npre, "\n"))

    xp    <- as.matrix(pre(x, list_pre[[j]]))
    xp_ok <- xp[ok, , drop = FALSE]

    # Grille gamma : relative à 1/p si non spécifiée
    gv <- if (is.null(gamma_grid)) c(0.1, 1, 10) / ncol(xp_ok) else gamma_grid
    ng <- length(gv)

    # Grille R² 3D : coût × epsilon × gamma
    r2_mat  <- array(NA_real_, dim = c(nc, ne, ng))
    sep_mat <- array(NA_real_, dim = c(nc, ne, ng))

    for (gi in seq_len(ng)) {
      for (ei in seq_len(ne)) {
        for (ci in seq_len(nc)) {
          if (verb) cat(sprintf("  C=%.3g  eps=%.3g  gamma=%.3g\n",
                                cost_grid[ci], epsilon_grid[ei], gv[gi]))

          yp_cv <- rep(NA_real_, n_ok)
          for (k in seq_len(K)) {
            test_idx  <- segm_folds[[k]]
            train_idx <- setdiff(seq_len(n_ok), test_idx)
            fm <- svm(xp_ok[train_idx, , drop = FALSE],
                      yv_ok[train_idx],
                      kernel  = kernel,
                      cost    = cost_grid[ci],
                      epsilon = epsilon_grid[ei],
                      gamma   = gv[gi],
                      scale   = TRUE)
            yp_cv[test_idx] <- predict(fm, xp_ok[test_idx, , drop = FALSE])
          }

          r2_mat[ci, ei, gi]  <- r2_fn(yp_cv, yv_ok)
          sep_mat[ci, ei, gi] <- sqrt(mean((yv_ok - yp_cv)^2))
        }
      }
    }

    # Meilleur combo pour ce prétraitement
    best_idx         <- which(r2_mat == max(r2_mat, na.rm = TRUE),
                              arr.ind = TRUE)[1, ]
    r2_best[j]       <- r2_mat[best_idx[1], best_idx[2], best_idx[3]]
    sep_best[j]      <- sep_mat[best_idx[1], best_idx[2], best_idx[3]]
    params_best[[j]] <- list(
      cost    = cost_grid[best_idx[1]],
      epsilon = epsilon_grid[best_idx[2]],
      gamma   = gv[best_idx[3]]
    )

    # Courbe R² vs log10(C) pour ce prétraitement — max sur (epsilon, gamma)
    r2_by_cost[, j] <- apply(r2_mat, 1, max, na.rm = TRUE)

    pf      <- t(list_pre[[j]])
    dim(pf) <- c(1, 2 * ncol(pf))
    pftot   <- c(pftot, paste0(pf, collapse = "_"))
  }

  # --- Meilleur prétraitement global ---
  best_j   <- which.max(r2_best)
  best_par <- params_best[[best_j]]
  txtpre   <- pftot[best_j]

  cat(
    "Pré    :", txtpre,                        "\n\n",
    "C      :", best_par$cost,                 "\n\n",
    "Epsilon:", best_par$epsilon,              "\n\n",
    "Gamma  :", round(best_par$gamma, 6),      "\n\n",
    "R2     :", round(r2_best[best_j], 2),     "\n\n",
    "SEP    :", round(sep_best[best_j], 2),    "\n\n\n"
  )

  # --- Plot CV : R² vs log10(C), une ligne par prétraitement ---
  if (plotCV) {
    matplot(log10(cost_grid), r2_by_cost,
            type = "l", lty = 1, col = seq_len(npre),
            xlab = "log10(C)", ylab = "R² Validation Croisée (SVR)",
            ylim = range(r2_by_cost, na.rm = TRUE))
    abline(v = log10(best_par$cost), col = "grey60", lty = 2)
    legend("bottomright",
           legend = pftot,
           col    = seq_len(npre),
           lty    = 1, cex = 0.8, bg = "white")
    title(titl)
  }

  # --- Prédictions CV du meilleur modèle (pour plot YY) ---
  if (plotYY) {
    xp_best <- as.matrix(pre(x, list_pre[[best_j]]))
    xp_ok   <- xp_best[ok, , drop = FALSE]

    yp_cv_best <- rep(NA_real_, n_ok)
    for (k in seq_len(K)) {
      test_idx  <- segm_folds[[k]]
      train_idx <- setdiff(seq_len(n_ok), test_idx)
      fm <- svm(xp_ok[train_idx, , drop = FALSE],
                yv_ok[train_idx],
                kernel  = kernel,
                cost    = best_par$cost,
                epsilon = best_par$epsilon,
                gamma   = best_par$gamma,
                scale   = TRUE)
      yp_cv_best[test_idx] <- predict(fm, xp_ok[test_idx, , drop = FALSE])
    }

    plot(yp_cv_best, yv_ok,
         xlab = "Teneur Prédite par SPIR (SVR)",
         ylab = "Teneur Mesurée")
    abline(a = 0, b = 1, col = "red", lty = 2)
    fit       <- lm(yp_cv_best ~ yv_ok)
    abline(fit, col = "blue")
    r_squared <- summary(fit)$r.squared
    legend("topleft",
           legend = c(
             "y = x",
             bquote(Validation_Croisée: ~ R^2 == .(round(r_squared, 2))),
             bquote(C:       .(best_par$cost)),
             bquote(epsilon: .(best_par$epsilon)),
             bquote(gamma:   .(round(best_par$gamma, 5))),
             bquote(n_ech:   .(n_ok))
           ),
           col = c("red", "blue", "white", "white", "white", "white"),
           lty = c(2, 1), bty = "n")
    title(paste("", titl, sep = " - "))
  }

  invisible(list(
    trait   = titl,
    pre     = txtpre,
    cost    = best_par$cost,
    epsilon = best_par$epsilon,
    gamma   = best_par$gamma,
    r2      = round(r2_best[best_j], 4),
    sep     = round(sep_best[best_j], 4)
  ))
}
