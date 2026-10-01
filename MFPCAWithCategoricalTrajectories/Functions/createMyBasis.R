createMyBasis <- function(argvals, type = "spline", n = 5, knots = NULL,
                          order = 3, wavelet = "DaubExPhase", filter.number = 2,
                          bc = "interval",custom_mat=NULL) {
  library(wavethresh)
  library(splines)
  if (!type %in% c("spline", "wavelet","custom"))
    stop("'type' must be either 'spline' or 'wavelet'.")
  if (length(argvals) < 2)
    stop("'argvals' must contain at least two points.")
  if (is.null(n) || length(n) != 1 || n < 1 || n != round(n))
    stop("'n' must be a positive integer.")
  n <- as.integer(n)
  
  # ==========================================================
  # MODE SPLINE
  # ==========================================================
  if (type == "spline") {
    if (is.null(knots)) {
      if(!is.null(n))
      {
        knots1=seq(min(argvals),max(argvals),length.out=n+2)
        knots2=sort(unique(knots))
        knots=knots2[!knots2%in%c(min(argvals),max(argvals))]
        print(knots)
      }
    
      else stop("Please enter n or knots")
      
    } else {
      knots <- sort(unique(knots))
      if (any(knots <= min(argvals) | knots >= max(argvals)))
        stop("All knots must be strictly inside the range of 'argvals'.")
      if (anyDuplicated(knots))
        stop("Knots must be unique.")
    }
    
    if (order == 0) {
      breaks <- c(min(argvals), knots, max(argvals))
      Bmat <- sapply(seq_len(length(breaks) - 1), function(i) {
        if (i < length(breaks) - 1)
          as.numeric(argvals >= breaks[i] & argvals < breaks[i + 1])
        else
          as.numeric(argvals >= breaks[i] & argvals <= breaks[i + 1])
      })
    } else {
      B <- bs(argvals, knots = knots, degree = order, intercept = TRUE)
      Bmat <- matrix(as.numeric(B), nrow = nrow(B), ncol = ncol(B))
    }
    
    if (ncol(Bmat) != n)
      warning(paste0("The constructed spline basis contain exactly" ,ncol(Bmat)," functions."))
    
    return(funData(list(argvals), t(Bmat)))
  }
  
  # ==========================================================
  # MODE WAVELET
  # ==========================================================
  if (type == "wavelet") {
    
    N <- length(argvals)
    
    if (N != 2^round(log2(N)))
      stop("For wavelets, length(argvals) must be a power of 2. Use seq(xmin, xmax, length.out = 2^p).")
    
    if (!isTRUE(all.equal(diff(argvals), rep(diff(argvals)[1], N - 1))))
      stop("For wavelets, 'argvals' must be regularly spaced.")
    
    j <- log2(n + 1)
    
    if (!isTRUE(all.equal(j, round(j)))) {
      warning("'n' must be of the form 2^j - 1 (e.g. 1, 3, 7, 15, 31, 63, ...).")
      return(NULL)
    }
    
    j <- as.integer(round(j))
    J <- log2(N)
    
    if (n > N - 1)
      stop("'n' is too large for the number of evaluation points.")
    
    if (!bc %in% c("interval", "periodic", "symmetric"))
      stop("'bc' must be 'interval', 'periodic' or 'symmetric'.")
    
    if (bc == "interval" && (filter.number < 1 || filter.number > 8))
      stop("For 'bc = interval', 'filter.number' must be between 1 and 8.")
    
    cat("Wavelet basis\n---------------\n",
        "Number of points :", N, "\n",
        "Number of basis functions :", n, "\n",
        "Maximum retained level :", j, "\n",
        "Wavelet :", wavelet, "\n",
        "Filter number :", filter.number, "\n",
        "Boundary :", bc, "\n")
    
    # --------------------------------------------------------
    # Décomposition d'un signal nul pour obtenir la structure
    # des coefficients wavelet
    # --------------------------------------------------------
    d <- wd(rep(0, N), filter.number = filter.number,
            family = wavelet, bc = bc)
    
    # Pour bc = interval, transformed.vector contient :
    # 1 coefficient d'échelle + 1 + 2 + 4 + ... coefficients
    # de détail. Les n premières correspondent donc exactement
    # aux n fonctions voulues.
    if (bc == "interval") {
      
      Bmat <- matrix(0, nrow = N, ncol = n)
      
      for (i in seq_len(n)) {
        di <- d
        di$transformed.vector[] <- 0
        di$transformed.vector[i] <- 1
        Bmat[, i] <- wr(di)
      }
      
    } else {
      
      # Pour les autres conditions aux bords, on conserve ici
      # la construction par base canonique.
      Bfull <- matrix(0, nrow = N, ncol = N)
      
      for (i in seq_len(N)) {
        e <- numeric(N)
        e[i] <- 1
        di <- wd(e, filter.number = filter.number,
                 family = wavelet, bc = bc)
        Bfull[, i] <- wr(di)
      }
      
      Bmat <- Bfull[, seq_len(n), drop = FALSE]
    }
    
    basis <- funData(list(argvals), t(Bmat))
    
    # --------------------------------------------------------
    # Vérification numérique de l'orthonormalité
    # --------------------------------------------------------
    G <- crossprod(Bmat)
    max_off_diag <- max(abs(G - diag(diag(G))))
    diag_error <- max(abs(diag(G) - 1))
    
    cat("Maximum off-diagonal Gram value :", signif(max_off_diag, 4), "\n",
        "Maximum deviation from norm 1 :", signif(diag_error, 4), "\n")
    
    if (max_off_diag > 1e-8 || diag_error > 1e-8)
      warning("The constructed wavelet basis is not numerically orthonormal.")
    return(basis)
  }
  if(type=="custom")
  {
    if (is.null(custom_mat))
      stop("'custom_mat' must be provided when type = 'custom'.")
    
    if (!is.matrix(custom_mat))
      custom_mat <- as.matrix(custom_mat)
    
    N <- length(argvals)
    
    # --------------------------------------------------------
    # Vérification du nombre de points
    # --------------------------------------------------------
    if (nrow(custom_mat) != N) {
      stop(
        paste0(
          "'custom_mat' must have ", N,
          " rows (one row per value of 'argvals')."
        )
      )
    }
    
    # --------------------------------------------------------
    # Nombre de fonctions
    # --------------------------------------------------------
    K <- ncol(custom_mat)
    
    if (K < 1)
      stop("'custom_mat' must contain at least one function.")
    
    cat(
      "Custom basis\n",
      "------------\n",
      "Number of points :", N, "\n",
      "Number of basis functions :", K, "\n"
    )
    
    # --------------------------------------------------------
    # Construction du funData
    # --------------------------------------------------------
    basis <- funData(
      list(argvals),
      t(custom_mat)
    )
    
    return(basis)
  }
 }
# 