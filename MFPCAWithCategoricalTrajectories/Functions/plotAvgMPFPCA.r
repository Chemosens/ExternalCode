plotAvgMFPCA=function(resMFPCA)
{
  descriptors <- names(resMFPCA$meanFunction)
  meanFuncMFPCA <- resMFPCA$meanFunction
  df <- bind_rows(lapply(descriptors, function(descr) {
    data.frame(
      Time = meanFuncMFPCA[[descr]]@argvals[[1]],
      Probability = t(meanFuncMFPCA[[descr]]@X),
      Descriptor = descr
    )
  }))
  p=ggplot(df,
           aes(Time, Probability,
               colour = Descriptor,
               linetype = Descriptor)) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = colors) +
    scale_linetype_manual(values = line_type_nb) +
    labs(
      x = "Time",
      y = "Empirical probability"
    ) +
    theme_bw()
  return(p)
}
