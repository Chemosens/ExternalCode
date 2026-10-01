projectFunctionOnMyBasis=function(times_test,values,basis)
{
  f <- funData(list(times_test),matrix(values, nrow = 1)  )
  f_vec <- as.numeric(f@X)
  # matrice de la base
  B <- basis@X
  # pas temporel
  dt <- mean(diff(times_test))
  # matrice de Gram
  G <- B %*% t(B) * dt
  # produits scalaires <phi_j, f>
  b <- B %*% (f_vec * dt)
  # vrais coefficients de projection
  coeff <- solve(G, b)
  # reconstruction
  f_proj <- as.numeric(crossprod(coeff, B))
  df_plot <- data.frame(
    t = times_test,
    f = as.numeric(f@X),
    f_proj = f_proj
  )
  
  p=ggplot(df_plot, aes(x = t)) +
    geom_line(aes(y = f, color = "f(t)"), linewidth = 1) +
    geom_line(aes(y = f_proj, color = "Projection"), linewidth = 1) +
    labs(
      x = "t",
      y = "",
      color = ""
    ) +theme_bw()+
    theme(legend.position="none")
  return(list(proj=f_proj,p=p))
}