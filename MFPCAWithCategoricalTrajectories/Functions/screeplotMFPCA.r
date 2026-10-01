screeplotMFPCA=function(resMFPCA,returnValue=TRUE)
{

  df <- data.frame(
    Dimension = seq_along(resMFPCA$values),
    lambda = resMFPCA$values / sum(resMFPCA$values)
  )

  p=ggplot(df, aes(x = Dimension, y = lambda)) +
    geom_point(size = 2) +
    ylim(0, 0.2) +
    labs(
      title = "Eigenvalues (%)",
      x = "Dimension",
      y = expression(lambda)
    ) +
    theme_bw()
  if(!returnValue)
  {
    return(p)
  }
  else{return(df)}
}

