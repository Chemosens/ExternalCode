plotMyBasis <- function(basis,knots) {
  # Récupération des données
  argvals <- basis@argvals[[1]]
  X <- basis@X
  # Passage au format long
  df <- data.frame(
    t = rep(argvals, each = nrow(X)),
    basis = rep(seq_len(nrow(X)), times = length(argvals)),
    value = as.vector(X)
  )
  p=ggplot(df, aes(x = t, y = value, group = basis, color = factor(basis))) +
    geom_line(linewidth = 0.8) +
    labs(
      x = "t",
      y = "",
      color = "Basis"
    ) +
    theme(legend.position="none") + theme_bw()
  return(p)
}

