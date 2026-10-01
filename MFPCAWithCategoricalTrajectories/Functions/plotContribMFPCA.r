plotContribMFPCA=function(resMFPCA,imin=1,imax=3)
{
  eigenfunMFPCA <- resMFPCA$functions
  names(eigenfunMFPCA)
  contribAx1=data.frame()
  for (nom in names(eigenfunMFPCA))
  {
    for(i in imin:imax)
    {
      id=paste0(nom,i)
      contribAx1[id,"state"]=nom
      contribAx1[id,"dim"]=paste0("dim.",i)
      contribAx1[id,'value']=sum(eigenfunMFPCA[[nom]][i,]@X^2)/length(eigenfunMFPCA[[nom]][i,]@X)
    }
  }

  contribAx1_wide=data.frame()
  for (j in 1:length(names(eigenfunMFPCA)))
  {
    for(i in 1:5)
    {
      id=paste0(nom,i)
      nom=names(eigenfunMFPCA)[j]
      contribAx1_wide[j,i]=sum(eigenfunMFPCA[[nom]][i,]@X^2)/length(eigenfunMFPCA[[nom]][i,]@X)
    }
  }
  rownames(contribAx1_wide)=names(eigenfunMFPCA)
  colnames(contribAx1_wide)=paste0("dim ",1:ncol(contribAx1_wide))
  xtable(round(contribAx1_wide,digits=2))

  print(round(contribAx1_wide,digits=2))
  colorsForPlot=colors
  p=ggplot(contribAx1,aes(x=dim,y=value,fill=state))+geom_col()+scale_fill_manual(values=colorsForPlot)+theme_bw()
  return(p)
}

