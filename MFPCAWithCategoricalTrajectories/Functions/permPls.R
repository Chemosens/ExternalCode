plsPerm=function(pc_pls,product,selec=NULL,nperm=100,sparsity=c(1,1))
{
  if(!is.null(selec))
  {
    pc_pls=pc_pls[product%in%selec,]
    product=product[product%in%selec]
  }
  respls_opt=rgcca(blocks=list(pc=pc_pls,product=product),method="spls",sparsity=sparsity,ncomp=1,response=2)
  pred_opt=rgcca_predict(respls_opt)$score
  pred=rep(NA,nperm)
  for(perm in 1:nperm)
  { print(paste0(perm,"/",nperm))
    respls_perm=rgcca(blocks=list(pc=pc_pls,product=sample(product)),method="spls",response=2,ncomp=1,sparsity=sparsity)
    pred[perm]=rgcca_predict(respls_perm)$score
  }

  df_pred <- data.frame(
    accuracy = pred
  )
  p=ggplot(df_pred, aes(x = accuracy)) +
    geom_histogram(
      bins = 30,
      color = "black"
    ) +
    geom_vline(
      xintercept = pred_opt,
      color = "red",
      linewidth = 1
    ) +
    coord_cartesian(xlim = c(0.3, 1)) +
    labs(
      title = "Histogram of accuracies of permuted models\nand obtained model value",
      x = "Accuracy",
      y = "Count"
    ) +
    theme_minimal()
  return(list(p=p,pred=pred,pval=sum(pred>pred_opt)/length(pred),pls=respls_opt,pred=pred_opt))
}