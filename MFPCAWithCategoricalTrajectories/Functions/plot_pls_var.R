plotPlsVar=function(respls_opt,resMFPCA,dim_sel)
{
  W <- respls_opt$a[[1]][, dim_sel,drop=F]  # vecteur des coefficients PLS pour la composante choisie
  # Initialisation d'une liste pour stocker les données tidy
  df_list <- list()
  print(dim(W))
  print(W)
  ncomp_in_pls=dim(W)[1]
  # Boucle sur les attributs/fonctions
  for(att in names(resMFPCA$functions)) {
    phi <- resMFPCA$functions[[att]]
    # funData
    X <- phi@X[1:ncomp_in_pls,]
    print(dim(X))
    print(dim(W))# matrice (nb fonctions × nb points)
    # combinaison linéaire avec les poids PLS
    f_recon_values <- t(W) %*% X      # 1 × nb points
    # convertir en data frame tidy pour ggplot
    df <- data.frame(
      t = phi@argvals[[1]],
      value = as.numeric(f_recon_values),
      att = att
    )

    df_list[[att]] <- df
  }

  # Combiner toutes les variables en un seul data frame
  df_all <- bind_rows(df_list)

  # Plot ggplot
  p=ggplot(df_all, aes(x = t, y = value, color = att)) +
    geom_line(size = 1) +
    labs(
      title = paste("Reconstruction PLS dimension", dim_sel),
      x = "t",
      y = "Valeur",
      color = "Attribut"
    ) +
    scale_color_manual(values = colors) +
    theme_minimal()
  return(p)
}