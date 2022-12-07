#' @param tds tdsObject from chemosensR function tdsRead
#' @param option standardization option (usualStd or sameDurationMax -adding a stop descriptor)
#' @param endTime "max" or a value to be considered as the end of tasting
#' @param tdsdf a data.frame containing the following columns: c("descriptor","start","sequenceDuration","stop","sequenceStop")
convertForCfda=function(tds=NULL,option="usualStd",endTime="max",tdsdf=NULL)
{
  if(!is.null(tds)){tds_df0=tds$df}
  if(!is.null(tdsdf)){tds_df0=tdsdf}
  tds_df0[,"descriptor"]=as.factor(as.character(tds_df0[,"descriptor"]))
  tds_df2=tds_df0[,c("descriptor","start","sequenceDuration","stop","sequenceStop")]
  colnames(tds_df2)=c("state","time","duration","stop","seqStop")
  tds_df2[,"id"]=paste0(tds_df0[,"product"],tds_df0[,"subject"],tds_df0[,"rep"])
  if(endTime=="max")
  {
    globalEnding=max(tds_df2[,"duration"],na.rm=T)
  }
  else(globalEnding=endTime)
  indiv=levels(as.factor(tds_df2[,"id"]))
  tds_df3=tds_df2
  # Standardisation "violente"
  res_tds=NULL
  for(ind in indiv)
  {
    indiv_tds=tds_df2[tds_df2[,"id"]==ind,]
    indiv_tds[,"state"]=as.character(indiv_tds[,"state"])
    indiv_tds=rbind(indiv_tds,c("stop",indiv_tds[1,"seqStop"],indiv_tds[1,"duration"],indiv_tds[1,"seqStop"],indiv_tds[1,"seqStop"],ind))
    indiv_tds[,"time"]=as.numeric(as.character(indiv_tds[,"time"]))
    indiv_tds[,"duration"]=as.numeric(as.character(indiv_tds[,"time"]))
    indiv_tds[,"stop"]=as.numeric(as.character(indiv_tds[,"time"]))
    minim=min(indiv_tds[,"time"],na.rm=T)

    if(option=="usualStd")
    {
      maxim=max(indiv_tds[,"time"],na.rm=T)
      indiv_tds2=indiv_tds
      if(minim!=maxim)
      {
        indiv_tds2[,"time"]=(indiv_tds[,"time"]-minim)/(maxim-minim)
      }
      if(minim==maxim)
      {
        print(ind)
      }
    }
    if(option=="sameDurationMax")
    {
      indiv_tds2=as.data.frame(indiv_tds)
      indiv_tds2[,"time"]=(indiv_tds[,"time"]-minim)
      indiv_tds2=indiv_tds2[indiv_tds2[,"time"]<=globalEnding,]
      #    if(indiv_tds[1,"duration"]!=globalEnding)
      #    {
      line_max=c("stop",globalEnding,0,globalEnding,globalEnding,ind)
      indiv_tds2=rbind(indiv_tds2,line_max)
      #      indiv_tds2[,"time"]=as.numeric(as.character(indiv_tds2[,"time"]))/globalEnding
      indiv_tds2[,"time"]=as.numeric(as.character(indiv_tds2[,"time"]))

      #    }
    }

    res_tds=rbind(res_tds, indiv_tds2[,c("state","time","id")])
  }

  res_tds=res_tds[order(res_tds[,"id"],res_tds[,"time"]),]
  # tds_df3[is.na(tds_df3[,"time"]),]=0
  #  tds_df3=tds_df3[!is.na(tds_df3[,"tim"]),]
  return(res_tds)
}

#'Getting optimal encoding
#'@param tds_tdf4 dataframe
#'@param nbasis number of function basis
#'@param ncores number of cores to be used for parallelisation
#'@param text should the text be displayed in the individual graphs ?
#'@param colors vectors of colors to be used for descriptors
#'@param nBootstrap number of bootstraps
#'@param repel should the text be moved when collision ?
#'@return a list containing (i)fmca results (fmca) (ii) chosen basis (basis),
#'(iii) individual plot (p_comp) (iv)first harmonic plot (p_h1) (v)second harmonic plot (p_h2) (vi) eigenvalue graph (p_eig)
getOptEnc=function(tds_df4,nbasis=6,ncores=6,norder=3,text=TRUE,colors=NULL,nBootstrap=50,repel=FALSE)
{
  basis <- create.bspline.basis(c(0, max(tds_df4[,"time"])), nbasis = nbasis, norder = norder)
  fmca <- compute_optimal_encoding(tds_df4, basis, nCores = ncores,nBootstrap=nBootstrap)
   if(repel==TRUE)
  {
     df=fmca$pc
     colnames(df)=paste0("PC",1:ncol(df))
     df=as.data.frame(df)
     df[,"name"]=rownames(df)
     p_comp=plotComponent(fmca,addNames=F)+geom_vline(xintercept=0,color="darkgrey")+geom_hline(yintercept=0,color="darkgrey")
     p_comp=p_comp+geom_text_repel(x=df[,"PC1"],y=df[,"PC2"],label=df[,"name"])
   }
  if(repel==FALSE)
  {
    p_comp=plotComponent(fmca,addNames=text)+geom_vline(xintercept=0,color="darkgrey")+geom_hline(yintercept=0,color="darkgrey")

  }

  p_eig=plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE)
  p_harm1=plot(fmca,col=colors)
  p_h1=p_harm1+theme_bw()+ggtitle("Harm. 1")
  p_harm2=plot(fmca, harm=2,col=colors)
  p_h2=p_harm2+theme_bw()+ggtitle("Harm. 2")
  p_eig=p_eig+theme_bw()
  p_comp=p_comp+ggtitle("Individual map")+theme_bw()+geom_hline(yintercept=0,color="grey")+geom_vline(xintercept=0,color="grey")
  return(list(fmca=fmca,basis=basis,p_comp=p_comp,p_h1=p_h1,p_h2=p_h2,p_eig=p_eig))
}


# Returns trajectories related to a linear combination of harmonics
#'@param enc_all object from compute_optimal_enconding
#'@param vecteurPropre vector to be used
#'@param nbPts number of points to be considered in the harmonic curve
#'@param colors vector of colors for the descriptors (whose names are descriptors)
#'@param sizeLine size of the harmonic line : "none", "default"(number of citation),"max" or "manual"
#'@param sizeVec if sizeLine=="manual", vector containing the size of the descriptors (whose names are descriptors)
#'@param maxSize number indicating the maximal size of the line
trajectoiresPropres=function(enc_all,vecteurPropre=reslda$scaling[,paste0("LD",harm)],nbPts=100,colors,sizeLine="none",sizeVec=NULL,maxSize=4)
{
  fdlist=list()
  #if(is.null(vecteurPropre)){vecteurPropre=rep(1,ncomp)}
  for(i in 1:length(vecteurPropre))
  {
    alpha <- enc_all$fmca$alpha[[i]] # selection des alpha
    alpha[is.na(alpha)]=0
    fdObj <- fd(alpha, enc_all$fmca$basis)#utilisation de la fonction fd pour obtenir les fonctions propres
    rangex <- fdObj$basis$rangeval
    nBasis <- fdObj$basis$nbasis
    fdlist[[i]] <- eval.fd(seq(0,1,length=nbPts), fdObj)*vecteurPropre[i] #utilisation de la fonction eval.fd pour avoir un nombre discret de points
    if(i==1){matres=fdlist[[1]]}
    if(i>1){matres=matres+fdlist[[i]]}
  }
  # On obtient une matrice avec 100 lignes nbAttributs colonnes (valable pour l'axe 1)
  fdres=as.data.frame(matres);fdres[,"time"]=seq(0,1,length=100)
  ggfdmat=reshape(fdres,direction="long",varying=list(colnames(matres)),times=colnames(matres),timevar="descriptor",v.names="score")
  # p_exp=ggplot(ggfdmat,aes(x=time,y=score,color=descriptor))+geom_line()+scale_color_manual(values=colors)+theme_bw()
  #ggfdmat=reshape(fdres,direction="long",varying=list(colnames(matres)),times=colnames(matres),timevar="descriptor",v.names="score")

  print(sizeLine)
  if(sizeLine!="none")
  {
    df_size=getSizeTable(enc_all,sizeLine=sizeLine,sizeVec=sizeVec,maxSize=maxSize)
    ggfdmat=merge(ggfdmat,df_size,by="descriptor")
    p_exp=ggplot(ggfdmat[ggfdmat[,"descriptor"]!="stop",],aes(x=time,y=score,color=descriptor,size=size))+geom_line()+scale_color_manual(values=colors)+theme_bw()+ scale_size(range = c(0,maxSize))
  }
  else
  {
    p_exp=ggplot(ggfdmat[ggfdmat[,"descriptor"]!="stop",],aes(x=time,y=score,color=descriptor))+geom_line()+scale_color_manual(values=colors)+theme_bw()
  }
  listing=list(p_exp=p_exp,fdres=fdres,ggfdmat=ggfdmat)

  return(listing)
}

#' Plots the harmonics with different sizelines
#'@param enc_all object from compute_optimal_enconding
#'@param harm number of the considered harmonic
#'@param sizeLine size of the harmonic line : "none", "default"(number of citation),"max" or "manual"
#'@param sizeVec if sizeLine=="manual", vector containing the size of the descriptors (whose names are descriptors)
#'@param maxSize number indicating the maximal size of the line
plotHarm=function(enc_all,harm=1,colors,sizeLine="none",sizeVec=NULL,maxSize=4)
{

  fdlist=list()
  # if(class(enc_all)=="pca.fd")
  # {
  #   alpha=enc_all$harmonics$coefs[,harm]
  #   fdObj <- fd(alpha, enc_all$harmonics$basis)#utilisation de la fonction fd pour obtenir les fonctions propres
  # }
  if(class(enc_all)!="pca.fd")
  {
    alpha <- enc_all$fmca$alpha[[harm]] # selection des alpha
    fdObj <- fd(alpha, enc_all$fmca$basis)#utilisation de la fonction fd pour obtenir les fonctions propres
  }
  #if(is.null(vecteurPropre)){vecteurPropre=rep(1,ncomp)}
    alpha[is.na(alpha)]=0
    rangex <- fdObj$basis$rangeval
    nBasis <- fdObj$basis$nbasis
    fdlist[[harm]] <- eval.fd(seq(0,1,length=100), fdObj) #utilisation de la fonction eval.fd pour avoir un nombre discret de points
     matres=fdlist[[harm]]
  # On obtient une matrice avec 100 lignes nbAttributs colonnes (valable pour l'axe 1)
  fdres=as.data.frame(matres);fdres[,"time"]=seq(0,1,length=100)
  ggfdmat=reshape(fdres,direction="long",varying=list(colnames(matres)),times=colnames(matres),timevar="descriptor",v.names="score")
     if(sizeLine!="none")
     {
       df_size=getSizeTable(enc_all,sizeLine=sizeLine,sizeVec=sizeVec,maxSize=maxSize)
       ggfdmat=merge(ggfdmat,df_size,by="descriptor")
       p_exp=ggplot(ggfdmat[ggfdmat[,"descriptor"]!="stop",],aes(x=time,y=score,color=descriptor,size=size))+geom_line()+scale_color_manual(values=colors)+theme_bw()+ scale_size(range = c(0,maxSize))+geom_hline(yintercept=0,color="grey")
     }
     else
     {
       p_exp=ggplot(ggfdmat[ggfdmat[,"descriptor"]!="stop",],aes(x=time,y=score,color=descriptor))+geom_line()+scale_color_manual(values=colors)+theme_bw()+geom_hline(yintercept=0,color="grey")
     }
  listing=list(p_exp=p_exp,fdres=fdres,ggfdmat=ggfdmat)
  return(listing)
}

#' Returns the size of line to be used in plotHarm or trajectoirePropres
#'@param enc_all object from compute_optimal_enconding
#'@param sizeLine size of the harmonic line : "none", "default"(number of citation),"max" or "manual"
#'@param sizeVec if sizeLine=="manual", vector containing the size of the descriptors (whose names are descriptors)
#'@param maxSize number indicating the maximal size of the line
getSizeTable=function(enc_all,sizeLine,sizeVec,maxSize)
{
  if(sizeLine=="manual")
  {
    if(names(sizeVec)!=colnames(enc_all$fmca$alpha[[1]]))
    {stop("size should have the same names as colnames alpha")}else
    {
      df_size=data.frame(descriptor=names(sizeVec),size=sizeVec)
    }
  }

  if(sizeLine=="max")
  {
    size=apply(enc_all$fmca$pt$pt,1,max)
    df_size=data.frame(descriptor=names(size),size=size)
  }
  if(sizeLine=="default")
  {
    size=apply(enc_all$fmca$pt$pt,1,sum)
    df_size=data.frame(descriptor=names(size),size=size)
  }
  # if(sizeLine=="maxSumByProd")
  # {
  #   if(nrow(enc_all$fmca$pt$pt)!=length(products)){ stop("Please enter product parameter with length equalling the number of lines in enc_all$fmca$pt")}
  #   for(prod in unique(products))
  #   {
  #     sizeProd=apply(enc_all$fmca$pt$pt[products==prod,],1,sum)
  #   }
  # }
  if(sizeLine=="none"){
    size=apply(enc_all$fmca$pt$pt,1,sum)
    df_size=data.frame(descriptor=names(size),size=1)
  }
  return(df_size)
}
#' @param ncomp: number of component to use in the reconstruction
#' @param enc_all result of encAll function
#' @param prod_converted
reconstructBarplot=function(ncomp,enc_all,prod_converted)
{
  matres=list()
  for(k in 1:ncomp)
  {
    fdObj <- fd(enc_all$fmca$alpha[[k]], enc_all$fmca$basis)#utilisation de la fonction fd pour obtenir les fonctions propres
    matres[[k]] <- eval.fd(seq(0,1,length=100), fdObj) #utilisation de la fonction eval.fd pour avoir un nombre discret de points
    #fdres=as.data.frame(matres);fdres[[k]][,"time"]=seq(0,1,length=100)
  }
  mat_tot=list()
  p_exp=list()
  real_att_j=att_j=matrix(NA,100,nrow(enc_all$fmca$pc))
  df_res=data.frame()
  #  att_j=real_att_j=matrix(NA,dim(mat_tot_i)[1],nrow(enc_all$fmca$pc))
  for(i in 1:nrow(enc_all$fmca$pc))# Sum on the subjects
  {
    mat_tot_i=matrix(0,dim(matres[[k]])[1],dim(matres[[k]])[2])
    for(k in 1:ncomp)
    {
      mat_tot_i=mat_tot_i+enc_all$fmca$pc[i,k]*matres[[k]]
    }
    mat_tot[[i]]=mat_tot_i
    fdres_tot=as.data.frame(mat_tot_i);fdres_tot[,"time"]=seq(0,1,length=100)
    rshp=reshape(fdres_tot,direction="long",idvar="time",timevar="descr",v.names="score",varying=list(colnames(mat_tot_i)),times=colnames(mat_tot_i))
    p_exp[[i]]=ggplot(rshp,aes(x=time,y=score,color=descr))+geom_line()+theme_bw()+geom_hline(yintercept=0,color="grey")
    prod_i=prod_converted[prod_converted[,"id"]==rownames(enc_all$fmca$pc)[i],]
    prod_i=prod_i[order(prod_i[,"time"]),]
    for(j in 1:(dim(mat_tot_i)[1]))# For each time
    {
      att_j[j,i]=colnames(mat_tot_i)[which.max(mat_tot_i[j,])]
      time2=prod_i[-1,"time"] 
      if(j!=dim(mat_tot_i)[1]-1)
      {
        indice_time=which(prod_i[-nrow(prod_i),"time"]<=fdres_tot[,"time"][j]&time2>=fdres_tot[,"time"][j])
      }
      else
      {
        indice_time=dim(prod_i)[1]
      }
      real_att_j[j,i]=prod_i[indice_time,"state"]
    }
    df_i=data.frame(id=rownames(enc_all$fmca$pc)[i],time=fdres_tot[,"time"],state=att_j[,i])
    df_res=rbind(df_res,df_i)
  }
  names(mat_tot)=rownames(enc_all$fmca$pc)
  return(list(df_res=df_res,real_att_j=real_att_j,att_j=att_j,mat_tot=mat_tot,p_exp=p_exp))
}
#' Requires chemosensR for reading TDS data
#' 
#' @param res result of reconstructBarplot
#'@param indices indexes of res to be recontructeds
getTdsObjectFromReconstruction=function(res,indices)
{
  df_resAll=data.frame()
  for(i in indices)
  {
    df_res=res[[i]]$df_res
    df_res[,"rep"]=1
    if(i<10){i0=paste0(0,i)}else{i0=i}
    df_res[,"product"]=paste0("N",i0)
    df_res[,"subject"]=substr(df_res[,"id"],3,5)
    df_res[,"descriptor"]=df_res[,"state"]
    df_res[,"score"]=1
    df_res[,"period"]=1
    df_res=df_res[,c("subject","product","time","descriptor","score","rep")]
    for(subj in unique(df_res[,"subject"]) )
    {
      for(prod in unique(df_res[,"product"]))
      {
        dfStartStop=data.frame(subject=rep(subj,2),product=rep(prod,2),time=c(0,1),descriptor=c("START","STOP"),score=c(1,1),rep=c(1,1))
        df_res=rbind(df_res,dfStartStop)
      }
    }
    df_resAll=rbind(df_resAll,df_res)
  }
  tdsFictive=tdsRead(df=df_resAll,discretization=0.01)
  return(list(tds=tdsFictive,df=df_resAll))
}