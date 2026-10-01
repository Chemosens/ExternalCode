#' @title Read a file and returns a tds object.
#' @description Read a file and returns a tds object.
#' @param file Full path to text delimited file or xls/xlsx file. Imperative columns: product/subject/text.
#' @param df If not NULL, this data.frame with imperative column product, subject, rep is used to produce a tds object
#' @param sep Field separator if the file is a text delimited file.
#' @param cols List of delimited text columns names associated to the expected columns names.
#' @param supCols A vector of additional columns to be included, for example rep or period.
#' @param startWithFirstCitation Logical. If TRUE, time start with the first citation.
#' @param discretization Numeric. Pace of discretisation.
#' @param periods  Numeric. Number of periods to split the sequences.
#' @return list of dataframes. \cr 
#' $canonicalData: product/rep/subject/descriptor/start/score/duration/stop/sequenceStart/sequenceStop/stdStart/stdDuration/stdStop \cr 
#' $dominances: product/rep/subject/descriptor/time/score/period \cr 
#' $standardizedDominances: same as $dominances with time between 0 (start) and 1 (stop) \cr 
#' $durations: product/rep/subject/descriptor/score/period \cr 
#' $citations: product/rep/subject/descriptor/score/period \cr 
#' $behaviour: product/rep/subject/variable/score
#' @export
tdsRead=function(file,df=NULL, sep=";", cols=list(subject="subject", product="product", descriptor="descriptor",time="time",score="score",rep="rep"), supCols="",startWithFirstCitation=FALSE,discretization=1,periods=1) {
  
  expandTimes=function(x,t1,t2,discretization) {
    rep=x[["rep"]]
    product=x[["product"]]
    subject=x[["subject"]]
    descriptor=x[["descriptor"]]
    score=x[["score"]]
    tt1=x[[t1]]
    tt2=x[[t2]]
    time=utils::head(seq(tt1,tt2,discretization),-1)
    m=expand.grid(rep=rep,product=product,subject=subject,descriptor=descriptor,score=score,time=time)
    return (m)
  }
  
  getDominance=function(df, start, stop, interval, maxTime) {
    dominance=do.call(rbind, apply(df,1,expandTimes,start,stop,interval))  
    dominance$score=as.numeric(as.character(dominance$score))
    completeDominance=expand.grid(rep=unique(df$rep), product=unique(df$product),  subject=unique(df$subject), descriptor=unique(df$descriptor),time=seq(0,maxTime,interval))
    dominance=merge(completeDominance,dominance,all=TRUE)
    dominance[is.na(dominance)]=0
    dominance=splitInPeriods(dominance,periods)
    dominance=dominance[dominance$period>0,]
  }
  
  descriptor=NULL
  if(is.null(df))
  {
    df = readXlsxOrCsv(file, sep)
  }

  # Données nettoyées
  res.cleanDf = cleanDf(df, list(list(dfName="subject", colName=cols[["subject"]], type="factor", clean=TRUE, complete=FALSE),list(dfName="product", colName=cols[["product"]], type="factor", clean=TRUE, complete=TRUE),list(dfName="descriptor", colName=cols[["descriptor"]], type="factor", clean=TRUE, complete=TRUE),list(dfName="rep", colName=cols[["rep"]], type="factor", clean=TRUE, complete=TRUE),list(dfName="score", colName=cols[["score"]], type="numeric", clean=FALSE, complete=FALSE),list(dfName="time", colName=cols[["time"]], type="numeric", clean=FALSE, complete=FALSE)))
  df = res.cleanDf$completeDf

  # Tri par subject, rep, product, time
  df=df[order(df$subject,df$rep,df$product,df$time), ]

  # Vérification que les durées existent
  if(dim(df[df$score>0,])[1]==0){
    stop("No data with positive time")
  }
  
  # Vérification que les scores existent
  if(dim(df[df$score>0,])[1]==0){
    stop("No data with positive score")
  }
  
  # Suppression des double Stop
  stops=which(df$descriptor=="STOP")
  
  if(length(stops)==0) {
    stop("No stop in the TDS data")
  }
  
  toRemove=c()
  if(length(stops)!=1)  {
    for(i in 1:(length(stops)-1)) {
      if(stops[i+1]-stops[i]==1) { 
        toRemove=c(toRemove,stops[i+1])
      }
    }
  }
  
  if (length(toRemove)>0)   {
    df=df[-toRemove,]
    warning(paste("Double STOP removed: ",length(toRemove),sep=""))
  }
  
  # Début (start), fin (stop) et durées de dominance (duration) des descripteurs au sein d'une séquence
  indTime=which(colnames(df)=="time")
  colnames(df)[indTime]="start"
  df$duration=c(df$start[-1],0) - df$start
  df$stop=df$start+df$duration
  df=subset(df, (descriptor !="START")&(descriptor!="STOP"))
  
  # Début (start), fin (stop) et durées de dominance (duration) des séquences
  sequenceStart=aggregate(start~subject+product+rep,df[df$score>0,],min)
  colnames(sequenceStart)[4]="sequenceStart"
  df=merge(df,sequenceStart)
  sequenceStop=aggregate(stop~subject+product+rep,df[df$score>0,],max)
  colnames(sequenceStop)[4]="sequenceStop"
  df=merge(df,sequenceStop)
  df$sequenceDuration=df$sequenceStop-df$sequenceStart
  if (startWithFirstCitation==TRUE) 
  {
    #Séquences commencent au premier clic, pas au start
    df$start=df$start-df$sequenceStart
    df$stop=df$stop-df$sequenceStart
    df$sequenceStart = 0
    df$sequenceStop=df$sequenceDuration
  }
  
  df$start=plyr::round_any(df$start,discretization)
  df$stop=plyr::round_any(df$stop,discretization)
  df$stdStart=round(df$start/df$sequenceStop,2)
  df$stdDuration=round(df$duration/df$sequenceStop,2)
  df$stdStop=df$stdStart+df$stdDuration
  
  df=df[df$stdStop>df$stdStart,]
  
 
  # Dominances 0/1, pour tous les sujet/produit/descripteur/rép/temps
  dominances=getDominance(df,"start","stop", discretization, round(max(df$stop,na.rm=TRUE)))
  stdDominances=getDominance(df,"stdStart","stdStop", 0.01, 1)
  
  
  # /!\
  # droplevels(dominances)
  # droplevels(stdDominances)
   print(dominances)
  # Durées de dominance par sujet/produit/descripteur/rép/période
  durations=aggregate(score~subject+descriptor+product+rep+period,dominances,function(x) {return (sum(x)*discretization)} )
  
  # Citations par sujet/produit/descripteur/rép/période
  citations=durations
  citations[citations$score>0,"score"]=1
  
  # Comportement
  behaviour=aggregate(sequenceStart~subject+product+rep,df[df$score>0,],min)
  colnames(behaviour)[4]="sequenceStart"
  sequenceStop=aggregate(sequenceStop~subject+product+rep,df[df$score>0,],max)
  colnames(sequenceStop)[4]="sequenceStop"
  behaviour=merge(behaviour,sequenceStop)
  behaviour$sequenceDuration=behaviour$sequenceStop-behaviour$sequenceStart
  distinctCitations=aggregate(score~subject+product+rep,unique(df[df$score>0,c("subject","product","rep","descriptor","score")]),length)
  colnames(distinctCitations)[4]="descriptors"
  behaviour=merge(behaviour,distinctCitations)
  totalCitations=aggregate(score~subject+product+rep,df,length)
  colnames(totalCitations)[4]="citations"
  behaviour=merge(behaviour,totalCitations)
  meanDuration=aggregate(duration~subject+product+rep,df[df$score>0,],mean)
  colnames(meanDuration)[4]="meanDominanceDuration"
  behaviour=merge(behaviour,meanDuration)
  behaviour$meanDominanceDuration=round(behaviour$meanDominanceDuration,2)
  behaviour=reshape2::melt(data = behaviour, id.vars = c("subject","product","rep"))
  colnames(behaviour)[c(4,5)]=c("variable","score")
  
  # Summary
  summary=list()
  summary$products=unique(as.character(citations$product))
  summary$subjects=unique(as.character(citations$subject))
  summary$descriptors=unique(as.character(citations$descriptor))
  tmp=aggregate(score~subject+product+rep,citations,sum)
  summary$evaluations=aggregate(score~product+subject,tmp[tmp$score>0,],length)
  colnames(summary$evaluations)[3]="n"
  
  
  # Tableau produit / rep + stats sur nombre d'évaluations, de citations et durées
  tmpC=aggregate(score~subject+product+rep,citations,sum)
  productTable=aggregate(score~product,tmpC,mean)
  colnames(productTable)[2]="meanCitations"
  tmp=aggregate(score~product,tmpC,sd)
  productTable=merge(productTable,tmp)
  colnames(productTable)[3]="sdCitations"
  tmpD=aggregate(score~subject+product+rep,durations,sum)
  tmp=aggregate(score~product,tmpD,mean)
  productTable=merge(productTable,tmp)
  colnames(productTable)[4]="meanDuration"
  tmp=aggregate(score~product,tmpD,sd)
  productTable=merge(productTable,tmp)
  colnames(productTable)[5]="sdDuration"
  productTable[,2:5]=round(productTable[,2:5],2)
  tmp=aggregate(score~product,tmpC[tmpC$score>0,],length)
  productTable=merge(productTable,tmp)
  colnames(productTable)[6]="n"
  productTable
  
  # Tableau descripteurs
  tmpA=aggregate(score~subject+product+descriptor+rep,citations,sum)
  tmpA=aggregate(score~product+descriptor+rep,tmpA,mean)
  tmpA=aggregate(score~product+descriptor,tmpA,mean)
  attributeTable=aggregate(score~descriptor,tmpA,mean)
  colnames(attributeTable)[2]="freq"
  tmpA=aggregate(score~subject+product+descriptor+rep,durations,sum)
  tmpA=aggregate(score~product+descriptor+rep,tmpA,mean)
  tmpA=aggregate(score~product+descriptor,tmpA,mean)
  tmpA=aggregate(score~descriptor,tmpA,mean)
  attributeTable=merge(attributeTable,tmpA)
  colnames(attributeTable)[3]="meanDominanceDuration"
  tmpA=aggregate(score~subject+product+descriptor+rep,durations,sum)
  tmpA=aggregate(score~product+descriptor+rep,tmpA[tmpA$score>0,],mean)
  tmpA=aggregate(score~product+descriptor,tmpA,mean)
  tmpA=aggregate(score~descriptor,tmpA,mean)
  attributeTable=merge(attributeTable,tmpA)
  colnames(attributeTable)[4]="meanSojournTime"
  attributeTable[,2:4]=round(attributeTable[,2:4],2)
  attributeTable
  

  result=list(df=df,dominances=dominances,stdDominances=stdDominances, durations=durations, citations=citations, behaviours=behaviour, summary=summary, productInfo=productTable, attributeInfo=attributeTable)
  class(result)="tds"
  
  return (result)
}