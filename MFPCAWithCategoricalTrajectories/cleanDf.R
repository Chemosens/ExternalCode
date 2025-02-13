#' @title Clean a dataframe.
#' @description Check and clean a dataframe: column name and type check, removal of special characters, removal of duplicates, completion of missing values by NA.
#' @param df A dataframe.
#' @param cols A list of list, elements of the list being colName, dfName, type, clean, complete.
#' @return The cleaned dataframe.
#' @keywords internal
cleanDf = function(df, cols)
{
  colNames=colnames(df)
  colToComplete = list()
  
  for (i in 1:length(cols)) {
    colName=cols[[i]][["colName"]]
    dfName=cols[[i]][["dfName"]]
    type=cols[[i]][["type"]]
    clean=cols[[i]][["clean"]]
    complete=cols[[i]][["complete"]]
  
    # Colonne présente dans le fichier ?  
    ind=which(colNames==colName)
    if (length(ind)==0) {
      stop(paste("Missing '"),colName, "' column.",sep="")
    }
    
    # Renommage
    colnames(df)[ind]=dfName
  
    # Suppression des caractères spéciaux
    if (clean == TRUE) {
    #  df[,dfName]=cleanString(df[,dfName])
    }
    
    if (complete == TRUE) {
      colToComplete[[dfName]] = as.factor(unique(df[,dfName]))
    }
    
    # Typage
    if (type=="factor") {
      df[,dfName]=as.factor(df[,dfName])
    }
    if (type=="numeric") {
      df[,dfName]=as.numeric(as.character(df[,dfName]))
    }
    if (type=="character") {
      df[,dfName]=as.character(df[,dfName])
    }

  }
  
  # Suppression des doublons
  df=unique(df)
  
  # Complétion des données manquantes par des NA
  if (length(colToComplete)>0) {
    completeDf=as.data.frame(colToComplete[[1]])
    names(completeDf)=names(colToComplete)[[1]]
    
    for (i in 1:length(colToComplete)) {
      tmp=as.data.frame(colToComplete[[i]])
      names(tmp)=names(colToComplete)[[i]]
      completeDf=merge(completeDf,tmp)
    }
    
    completeDf=merge(completeDf,df,by=names(completeDf),all.x=TRUE)
  }

  return (list(completeDf=completeDf))
  
  #TODO : définir max données manquantes
  #TODO : completeMissingValues -> remplacement des NA par valeurs imputées ?
  #TODO : log (nombre, missing, etc)
  
}