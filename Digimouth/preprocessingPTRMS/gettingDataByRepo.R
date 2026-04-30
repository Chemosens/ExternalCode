# digimouth project
#=================
# This pre-processes Digimouth data for gustometer
rm(list=ls())


# Parameters

repo="20250923"
duration=0.1000083 #0.5000145 

# End ====================
library(ggplot2)
library(PTRMSR)
library(plotly)
library(openxlsx)
nameFile=paste0("ptr_",repo)
repoData=paste0("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/DonneesBrutes/",repo)
saveRepo="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/PreprocessedData"
paramList=c("date","repo","repoData","saveRepo","duration")
date=format(Sys.time(), "%Y-%m-%d %H:%M:%S")
paramValues=c(date,repo,repoData,saveRepo,duration)
df_param=data.frame(param=paramList,value=paramValues)


setwd(repoData)
md=F;it=F
if("metaData.xlsx" %in% list.files()){md=T;metaData=read.xlsx("metaData.xlsx",sep=";")}
if("integrationTable.csv"%in% list.files()){it=T;integrationTable=read.csv("integrationTable.csv",sep=";")}

#metaData=read.csv("metaData.csv",sep=";")
filesIA=list.files(pattern="*.h5")
t0=Sys.time()
ptrRes=ptrReadListIBT(filesIA,sumSpectraOnly=FALSE,mz=NULL,rt=NULL,integrationTable=integrationTable)
t1=Sys.time()
colnames(ptrRes)[4]="ion"
ptrRes[,"duration"]=duration


# Saving excel
setwd(saveRepo)
wb <- createWorkbook()
addWorksheet(wb, "PTR") # Feuille de données

#addWorksheet(wb, "metaData") 
addWorksheet(wb, "param") 
writeData(wb, "PTR", ptrRes)
if(it){addWorksheet(wb, "integrationTable") ;writeData(wb, "integrationTable", integrationTable)}
if(md){addWorksheet(wb, "metaData");writeData(wb, "metaData", metaData)}
writeData(wb, "param", df_param)
saveWorkbook(wb, file=paste0(nameFile,".xlsx"), overwrite = TRUE)


# saving RData
listPtrResults=list(ptr=ptrRes,date_analysis=date,date_data=repo,integrationTable=integrationTable)
save(listPtrResults,file=paste0(nameFile,".RData"))

