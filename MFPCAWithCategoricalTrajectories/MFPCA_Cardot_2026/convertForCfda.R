#' @param tds result of tdsRead
#' @param option "usualStd" or "sameDurationMax"
#' @param endTime "max" or a numerical value indicating the end of the tasting
#' @param tdsdf dataframe to be used instead of tds object (whose colanems are c("descriptor","start","sequenceDuration","stop","sequenceStop"))
convertForCfda=function(tds=NULL,option="usualStd",endTime="max",tdsdf=NULL)
{
  if(!is.null(tds)){tds_df0=tds$df}
  if(!is.null(tdsdf)){tds_df0=tdsdf}
  tds_df0[,"descriptor"]=as.factor(as.character(tds_df0[,"descriptor"]))
  tds_df2=tds_df0[,c("descriptor","start","sequenceDuration","stop","sequenceStop")]
  colnames(tds_df2)=c("state","time","duration","stop","seqStop")
  tds_df2[,"id"]=paste0(tds_df0[,"product"],tds_df0[,"subject"])
  if(endTime=="max")
  {
    globalEnding=max(tds_df2[,"duration"],na.rm=T)
  }
  else(globalEnding=endTime)
  indiv=levels(as.factor(tds_df2[,"id"]))
  tds_df3=tds_df2
  # Standardisation
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
      line_max=c("stop",globalEnding,0,globalEnding,globalEnding,ind)
      indiv_tds2=rbind(indiv_tds2,line_max)
      indiv_tds2[,"time"]=as.numeric(as.character(indiv_tds2[,"time"]))
    }

    res_tds=rbind(res_tds, indiv_tds2[,c("state","time","id")])
  }
  res_tds=res_tds[order(res_tds[,"id"],res_tds[,"time"]),]
  return(res_tds)
}
