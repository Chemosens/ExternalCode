prepareForMFPCA=function(df_norm1,times.out)
{
  # Computes columns containing acid_0 or acid_1 for example
  df_norm1[,"state"]=paste0(df_norm1[,"descriptor"],"_",df_norm1[,"score"])
  nom_ind <- levels(as.factor( df_norm1$id))
  n_ind <- length(nom_ind) ### 150
  etats <- levels(as.factor(df_norm1$state)) ## modalites des etats
  df_norm1$id <- as.character(df_norm1$id) ## nom des individus
  df_norm1$state <- as.character(df_norm1$state) ## modalites des etats
  n_state <- length(unique(df_norm1$state))
  pos_state <- 1:n_state
  names(pos_state)=levels(factor(df_norm1$state))
  ### Initialisation
  df_mfpca <- list()
  for (i in 1:n_state){df_mfpca[[i]] <- matrix(0,nrow=n_ind, ncol=length(times.out));
  rownames(df_mfpca[[i]])=nom_ind;}
  names(df_mfpca) <- etats
  descriptors=unique(df_norm1[,"descriptor"])
  for (i in 1:n_ind)
  {
    for(descriptor in descriptors)
    {
      df_temp <- df_norm1[df_norm1$id==nom_ind[i]&df_norm1$descriptor==descriptor,]
      mat_temp <- getIndicatrices(df_temp$state,df_temp$time,times.out=times.out)
      state_temp <- unique(df_temp$state)
      for (j in 1:length(state_temp)){
        nom_etat <- state_temp[j]
        pos_liste <- pos_state[etats == nom_etat]
        df_mfpca[[pos_liste]][i,] = mat_temp[j,]
      }
    }
  }
  df_fundata <- list()
  descriptors_1=paste0(descriptors,"_1")
  for (desc in descriptors_1){
    df_fundata[[desc]] <- funData(argvals=times.out,X=df_mfpca[[desc]])
  }
  names(df_fundata)=descriptors_1
  df_multiFun <- multiFunData(df_fundata,values=names(df_fundata))
  names(df_multiFun) <- sub("_1$", "", names(df_multiFun))
  return(df_multiFun)
}
