#' @title Read a TCATA file (colnames are listed in cols)
#' @param df dataframe whose colnames are subject, product, descriptor, time, score and rep
#' @description Read a file and returns a tcata object.
#' @param file Full path to text delimited file or xls/xlsx file. Imperative columns: product/subject/text.
#' @export
tcataRead=function(df, cols=list(subject="subject", product="product", descriptor="descriptor",time="time",score="score",rep="rep")) 
{
  df[,"id"]=paste0(df[,cols[["subject"]]],"_",df[,cols[["product"]]],"_",df[,cols[["rep"]]])
  descriptors=levels(factor(df[,cols[["descriptor"]]]))
  descriptors_real=descriptors[!descriptors%in%c("START","STOP")]
  ids=unique(df[,"id"])
  df_res=data.frame()
  for(id in ids)
  {
    subject_id=df[df[,"id"]==id,cols[["subject"]]][1]
    product_id=df[df[,"id"]==id,cols[["product"]]][1]
    rep_id=df[df[,"id"]==id,cols[["rep"]]][1]
    df_to_add0=data.frame(subject=rep(subject_id,length(descriptors_real)),
                         product=rep(product_id,length(descriptors_real)),
                         descriptor=descriptors_real,
                         id=rep(id,length(descriptors_real)),
                         time=rep(0,length(descriptors_real)),
                         score=rep(0,length(descriptors_real))
    )
    df_to_add1=data.frame(subject=rep(subject_id,length(descriptors_real)),
                          product=rep(product_id,length(descriptors_real)),
                          descriptor=descriptors_real,
                          id=rep(id,length(descriptors_real)),
                          time=rep(1,length(descriptors_real)),
                          score=rep(0,length(descriptors_real))
    )
    df_id_to_add=df[df[,"id"]==id,c(cols[["subject"]],cols[["product"]],cols[["descriptor"]],"id",cols[["time"]],cols[["score"]])]
    colnames(df_id_to_add)=c("subject","product","descriptor","id","time","score")
    df_id_to_add[,"time"]=df_id_to_add[,"time"]/df_id_to_add[df_id_to_add[,"descriptor"]=="STOP","time"]
    df_id_to_add=df_id_to_add[df_id_to_add[,"descriptor"]%in%descriptors_real,]
    df_id=rbind(df_to_add0,df_id_to_add,df_to_add1)
    df_res=rbind(df_res,df_id)
  }
  return(list(df=df_res))
}
