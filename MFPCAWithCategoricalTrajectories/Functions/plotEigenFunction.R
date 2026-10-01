plotEigenFunction <- function(resMFPCA,
                              perception,
                              comp = 1,
                              fact = 0.3,
                              ylim = c(0, 0.78),colors=NULL)
{

  meanFuncMFPCA  <- resMFPCA$meanFunction
  eigenfunMFPCA  <- resMFPCA$functions

  df <- data.frame(
    Time = meanFuncMFPCA[[perception]]@argvals[[1]],
    Mean = t(meanFuncMFPCA[[perception]]@X),
    Plus = t(meanFuncMFPCA[[perception]]@X) +
      fact * t(eigenfunMFPCA[[perception]][comp]@X),
    Minus = t(meanFuncMFPCA[[perception]]@X) -
      fact * t(eigenfunMFPCA[[perception]][comp]@X)
  )

  p=ggplot(df, aes(Time)) +
    geom_line(aes(y = Mean, colour = "Mean"),
              linewidth = 1.2) +
    geom_line(aes(y = Plus, colour = "Plus"),
              linewidth = 1,
              linetype = "dashed") +
    geom_line(aes(y = Minus, colour = "Minus"),
              linewidth = 1,
              linetype = "dotted") +
    scale_colour_manual(
      breaks = c("Mean", "Plus", "Minus"),
      values = c(
        Mean = colors[perception],
        Plus = "black",
        Minus = "grey50"
      ),
      labels = c(
        expression(hat(p)),
        bquote(hat(p)+.(fact)*phi[.(comp)]),
        bquote(hat(p)-.(fact)*phi[.(comp)])
      ),
      name = NULL
    )+
    labs(
      title = perception,
      x = "Time",
      y = "Empirical probability"
    ) +
    theme_bw()
  p
  return(p)
}