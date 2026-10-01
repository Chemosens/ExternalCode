plotAllEigenFunction <- function(resMFPCA, comp = 1,colors=NULL)
{
  meanFuncMFPCA  <- resMFPCA$meanFunction
  eigenfunMFPCA  <- resMFPCA$functions

  all_df <- do.call(rbind, lapply(names(eigenfunMFPCA), function(perception) {

    Time <- meanFuncMFPCA[[perception]]@argvals[[1]]
    phi  <- t(eigenfunMFPCA[[perception]][comp]@X)

    data.frame(
      Time = Time,
      Eigen = phi,
      Perception = perception
    )
  }))

  p=ggplot(all_df, aes(x = Time, y = Eigen, colour = Perception)) +
    geom_line(linewidth = 1) +
    theme_bw() +
    labs(
      title = paste("Eigenfunctions - component", comp),
      x = "Time",
      y = "Eigenfunction"
    )+
    scale_color_manual(values=colors)
  return(p)
}