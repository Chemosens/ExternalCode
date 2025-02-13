#' @title Split a time variable in periods
#' @description Add a new column corresponding to a uniform time splitting.
#' @param df Dataframe
#' @param periods Numeric. Number of periods
#' @param cols Name of the grouping variables.
#' @return Dataframe with additional column
#' @keywords internal
splitInPeriods=function(df,periods, cols=c("rep","subject","product")) {
  # Découpage en périodes
  df$period=0
  
  f=as.formula(paste("time~",paste(cols, collapse="+")))
  
  start=aggregate(f,df,min)
  stop=aggregate(f,df,max)
  colnames(start)[4]="start"
  colnames(stop)[4]="stop"
  
  df=merge(df,start)
  df=merge(df,stop)
  
  df$incr=(df$stop-df$start)/periods
  for (p in 1:periods) {
    df[df$time>(df$start + (p-1)*df$incr) & df$time <=  (df$start + p*df$incr),"period"]=p
  }
  
  df$start=NULL
  df$stop=NULL
  df$incr=NULL
  
  return (df)
}