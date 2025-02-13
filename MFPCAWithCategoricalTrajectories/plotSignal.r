plotSignal <- function(signalSheet,signal,colors=NULL,linetype=NULL)
{
  data=signalSheet[signalSheet[,"Product"]==signal,]
  if(is.null(colors))
  {
    p=ggplot(data,aes(x=Time,y=Score,group=Attribute,color=Attribute,linetype=Attribute))+geom_line(size=2)+theme_bw()+ggtitle(signal)+scale_color_grey()
  }
  else
  {
    if(is.null(linetype))
    p=ggplot(data,aes(x=Time,y=Score,group=Attribute,color=Attribute))+geom_line(size=2)+ggtitle(signal)+scale_color_manual(values=colors)
    if(!is.null(linetype))
    {
      p=ggplot(data,aes(x=Time,y=Score,group=Attribute,color=Attribute,linetype=Attribute))+geom_line(size=2)+ggtitle(signal)+scale_color_manual(values=colors)+scale_linetype_manual(values=linetype)
      
    }
  }
  return(p)
}

plotSignal_nb <- function(signalSheet,signal)
{
  data <-signalSheet[signalSheet[,"Product"]==signal,]
  p <- ggplot(data,aes(x=Time,y=Score,group=Attribute,color=Attribute,linetype=Attribute))+geom_line(size=2)+theme_bw()+ggtitle(signal)+scale_color_grey()
  return(p)
}