# digimouth project
#=================
# This pre-processes Digimouth data for gustometer
rm(list=ls())


# Parameters


# End ====================
library(ggplot2)
library(PTRMSR)
library(plotly)
library(openxlsx)
library(gridExtra)

# Parameters
listRepo=c("20250602","20250603","20250604","20250605",
           "20250610","20250611","20250612","20250616",
           "20250617","20250618",
           "20250623","20250624","20250625",
           "20250630","20250701","20250702","20250703",
           "20250707","20250708","20250710",
           "20250722","20250723","20250724",
           "20250729","20250730","20250731",
           "20250805","20250806","20250807"
           )

# rq; des sol restant dans  15 et le 19
repo=listRepo[7] #"20250425" "20250514"
nameFilePTR=paste0("ptr_",repo,".xlsx")
nameFileSenso=paste0("senso_",repo,".xlsx")
nameFileInt=paste0("int_",repo,".xlsx")
nameFileCom=paste0("com_",repo,".xlsx")
repoDataPTR="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/PreprocessedData"
repoDataSenso="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/Senso/PreprocessedData"
xlsxPtr=read.xlsx(paste0(repoDataPTR,"/",nameFilePTR))
xlsxSenso=read.xlsx(paste0(repoDataSenso,"/",nameFileSenso))
xlsxInt=read.xlsx(paste0(repoDataSenso,"/",nameFileInt))
if(nameFileCom %in% list.files(repoDataSenso))
{
  xlsxCom=read.xlsx(paste0(repoDataSenso,"/",nameFileCom))
}else{xlsxCom=NULL}

if(repo=="20250425")
{
  xlsxSenso[xlsxSenso[,"subject"]=="S00","subject"]="S000"
  xlsxSenso[xlsxSenso[,"subject"]=="S001","subject"]="S001_err"
  xlsxSenso[xlsxSenso[,"subject"]=="S001_2","subject"]="S001"
  xlsxSenso[xlsxSenso[,"subject"]=="S002_2","subject"]="S001"
  xlsxPtr[,"file"]=gsub("-2testsolution", "testbonbon", xlsxPtr[,"file"])
  metaData=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/DonneesBrutes/20250425/metaData.xlsx")
}
if(repo=="20250514")
{
  metaData=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/DonneesBrutes/20250514/metaData.xlsx")
  xlsxSenso[xlsxSenso[,"rep"]=="1","rep"]=0
  xlsxSenso[xlsxSenso[,"rep"]=="2","rep"]=1
  xlsxSenso[xlsxSenso[,"rep"]=="3","rep"]=2
  xlsxSenso[xlsxSenso[,"rep"]=="4","rep"]=3
  xlsxSenso[xlsxSenso[,"rep"]=="5","rep"]=1
  xlsxSenso[xlsxSenso[,"rep"]=="6","rep"]=2
  xlsxSenso[xlsxSenso[,"rep"]=="7","rep"]=3
}
if(repo=="20250603")
{
  xlsxSenso[xlsxSenso[,"subject"]=="s331","subject"]="S331"
  xlsxCom[xlsxCom[,"subject"]=="s331","subject"]="S331"
  xlsxInt[xlsxInt[,"subject"]=="s331","subject"]="S331"
  
}
if(repo=="20250604")
{
  xlsxSenso[xlsxSenso[,"subject"]=="NO20","subject"]="N020"
  xlsxSenso[xlsxSenso[,"subject"]=="y075","subject"]="Y075"
  xlsxCom[xlsxCom[,"subject"]=="NO20","subject"]="N020"
  xlsxCom[xlsxCom[,"subject"]=="y075","subject"]="Y075"
  xlsxInt[xlsxInt[,"subject"]=="NO20","subject"]="N020"
  xlsxInt[xlsxInt[,"subject"]=="y075","subject"]="Y075"
  
}
if(repo=="20250624")
{
  xlsxSenso[xlsxSenso[,"subject"]=="p628","subject"]="P628"
  xlsxInt[xlsxInt[,"subject"]=="p628","subject"]="P628"
  xlsxCom[xlsxCom[,"subject"]=="p628","subject"]="P628"
}
if(repo=="20250708")
{
  xlsxSenso[xlsxSenso[,"subject"]=="963","subject"]="N963"
  xlsxInt[xlsxInt[,"subject"]=="963","subject"]="N963"
  xlsxCom[xlsxCom[,"subject"]=="963","subject"]="N963"
}
if(repo=="20250710")
{
  xlsxSenso[xlsxSenso[,"subject"]=="N2020","subject"]="N020"
  xlsxInt[xlsxInt[,"subject"]=="N2020","subject"]="N020"
  xlsxCom[xlsxCom[,"subject"]=="N2020","subject"]="N020"
}
if(repo=="20250723")
{
  xlsxSenso[xlsxSenso[,"subject"]=="c015","subject"]="C015"
  xlsxInt[xlsxInt[,"subject"]=="c015","subject"]="C015"
  xlsxCom[xlsxCom[,"subject"]=="c015","subject"]="C015"
  xlsxSenso[xlsxSenso[,"subject"]=="p628","subject"]="P628"
  xlsxInt[xlsxInt[,"subject"]=="p628","subject"]="P628"
  xlsxCom[xlsxCom[,"subject"]=="p628","subject"]="P628"
  
}
metaData=read.xlsx(paste0("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/DonneesBrutes/",repo,"/metaData.xlsx"))
metaData=metaData[!is.na(metaData[,"rep"]),]
metaData[metaData[,"fop"]=="nothing","fop"]="noth"
for (subject in unique(metaData[,"subject"]))
{
  for(product in unique(metaData[,"type"]))
  {
    fop_warmup=metaData[metaData[,"subject"]==subject&metaData[,"type"]==product&metaData[,"rep"]==4,"fop"]
    metaData[metaData[,"subject"]==subject &metaData[,"type"]==product& metaData[,"fop"]==fop_warmup ,"rep"]=metaData[metaData[,"subject"]==subject &metaData[,"type"]==product& metaData[,"fop"]==fop_warmup,"rep"]-1
  }
}

plotFile=function(xlsxPtr,xlsxSenso,subject,metaData,product,rep,fop,timeMEB=NULL,xlim=c(0,120),ylim=NULL)
{
  filePTR=metaData[!is.na(metaData[,"file"])&metaData[,"subject"]==subject&metaData[,"rep"]==rep&metaData["type"]==product&metaData[,"fop"]==fop,"file"]
  if(is.null(timeMEB))
  {
    timeMEB=metaData[!is.na(metaData[,"file"])&metaData[,"subject"]==subject&metaData[,"rep"]==rep&metaData["type"]==product&metaData[,"fop"]==fop,"start.(s)"]
    timeMEB=timeMEB[!is.na(timeMEB)]
  }
  filePTR=filePTR[!is.na(filePTR)]
  timeMEB=timeMEB[!is.na(timeMEB)]
  ptrDataIon=xlsxPtr[xlsxPtr[,"file"]==paste0(filePTR,".h5")&xlsxPtr[,"ion"]==ion,]
  ptrDataIon[,"time"]=ptrDataIon[,"time"]-timeMEB
  p=ggplot(ptrDataIon,aes(x=time,y=intensity))+geom_line()+theme_bw()
  dataSenso=xlsxSenso[xlsxSenso[,"product"]==product&xlsxSenso[,"rep"]==rep&xlsxSenso[,"fop"]==fop&xlsxSenso[,"subject"]==subject,]
  dataInt=xlsxInt[xlsxInt[,"product"]==product&xlsxInt[,"rep"]==rep&xlsxInt[,"fop"]==fop&xlsxInt[,"subject"]==subject,]
  com=NULL
  if(!is.null(xlsxCom))
  {
    dataCom=xlsxCom[xlsxCom[,"product"]==product&xlsxCom[,"rep"]==rep&xlsxCom[,"fop"]==paste0(fop,"_")&xlsxCom[,"subject"]==subject,]
    
    com=dataCom[,"com"]
  }
  int=round(dataInt[,"int"],digits=2)
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="chew","time"], linetype = "dashed", color = "red")
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="swallow","time"], linetype = "dotted", color = "blue")
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="swallowImposed","time"], color = "blue")
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="TotalSwallow","time"], color = "blue",linewidth=1.5)
  p=p+geom_vline(xintercept=30, color = "red",linewidth=1.5)
  p=p+ggtitle(paste0(ion,":",subject,", ",product,", rep ",rep, ",\n fop = ",fop, ",int=",int))+xlim(xlim[1],xlim[2])
  if(!is.null(com)){print(com);if(length(com)>0){p=p+xlab(paste0("time (",com,")"))}}
  if(!is.null(ylim)){p=p+ylim(ylim[1],ylim[2])}
  return(p)
}


# merging files
subjects=unique(metaData[,"subject"])
subjects
unique(xlsxSenso[,"subject"])
subject=subjects[2]

product="Sol"
if(product=="Sol"){fops=c("fast","long")}
if(product=="Gel"){fops=c("chew","succ")}
if(product=="Gus"){fops=c("cont","noth")}
ion="isoamylAcetate" #isoamylAcetate" "isoprene"  m71

ylim=c(0,max(xlsxPtr[xlsxPtr[,"ion"]==ion,"intensity"]))
ylim=c(0,max(xlsxPtr[xlsxPtr[,"ion"]==ion,"intensity"]))

p_list=list()
i=1
p_list[[1]]=plotFile(xlsxPtr=xlsxPtr,xlsxSenso=xlsxSenso,subject=subject,
                     metaData=metaData,product=product,rep=0,
                     fop=metaData[metaData[,"subject"]==subject&metaData[,"rep"]==0,"fop"][1],
                     timeMEB=NULL,xlim=c(0,130),ylim=ylim)
i=i+1
fop=fops[1]

for(rep in 1:3)
{
  p_list[[i]]=plotFile(xlsxPtr=xlsxPtr,xlsxSenso=xlsxSenso,subject=subject,metaData=metaData,product=product,rep=rep,fop=fop,timeMEB=NULL,xlim=c(0,130),ylim=ylim)
  i=i+1
}

grid.arrange(grobs=p_list,nrow=2)

p_list2=NULL
fop=fops[2]
i=1
for(rep in 1:3)
{
  p_list2[[i]]=plotFile(xlsxPtr=xlsxPtr,xlsxSenso=xlsxSenso,subject=subject,metaData=metaData,product=product,rep=rep,fop=fop,timeMEB=NULL,xlim=c(0,130),ylim=ylim)
  i=i+1
}

grid.arrange(grobs=p_list2,nrow=2)

