library(openxlsx)
library(tidyr)
library(dplyr)
listRepoSol=c("20250602","20250603","20250604","20250605",
              "20250610","20250611","20250612","20250616",
              "20250617","20250618")


listRepoGus=c( "20250623","20250624","20250625",
              "20250630","20250701","20250702","20250703",
              "20250707","20250708","20250710"
              )

listRepoGel=c( "20250722","20250723","20250724",
               "20250729","20250730","20250731",
               "20250805","20250806","20250807"
)
listRepoInVitro=c(listRepo=c( "20250407",      "20250408",
                              "20250409","20250410","20250411-vitro","20250416","20250424-vitro","20250428","20250429","20250430","20250512")
)

gettingInstrumentalDataFromListRepo=function(listRepo)
{
  bigIns=data.frame()
  for(i in 1:length(listRepo))
  {
    repo=listRepo[i] #"20250425" "20250514"
    print(repo)
    repoDataIns="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/PreprocessedData"

    nameFileIns=paste0("ins_",repo,".csv")
    metaData=read.xlsx(paste0("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/PTRMS/DonneesBrutes/",repo,"/metaData.xlsx"))
    metaData=metaData[!is.na(metaData[,"rep"]),]
    metaData[metaData[,"fop"]=="nothing","fop"]="noth"
    metaData=metaData[!is.na(metaData[,"subject"]),]
    for (subject in unique(metaData[,"subject"]))
    {
      fop_warmup=metaData[metaData[,"subject"]==subject&metaData[,"rep"]==4,"fop"]
      metaData[metaData[,"subject"]==subject & metaData[,"fop"]==fop_warmup,"rep"]=metaData[metaData[,"subject"]==subject & metaData[,"fop"]==fop_warmup,"rep"]-1
    }
    for(file in unique(metaData[,"file"]))
    {
      timeMEB=metaData[!is.na(metaData[,"file"])&metaData[,"file"]==file,"start.(s)"]
      timeMEB=timeMEB[!is.na(timeMEB)]
      xlsxIns=read.csv(paste0(repoDataIns,"/",nameFileIns),sep=";")
      #print(head(xlsxIns))
      insDataIon=xlsxIns[xlsxIns[,"file"]==paste0(file,".h5"),]
      #print(timeMEB)
      if(!is.null(timeMEB))
      {
        insDataIon[,"time"]=insDataIon[,"time"]-timeMEB
      }
      bigIns=rbind(bigIns,insDataIon)
    }
  }
  return(bigIns)
}

#==========================================
# Getting instrumental data for all data
#============================================
saveRepo="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf"
insGusto=gettingInstrumentalDataFromListRepo(listRepoGus)
write.csv(insGusto,file=paste0(saveRepo,"/insGusto.csv"),row.names=FALSE)

insSol=gettingInstrumentalDataFromListRepo(listRepoSol)
write.csv(insSol,file=paste0(saveRepo,"/insSol.csv"),row.names=FALSE)
insSol=read.csv(paste0(saveRepo,"/insSol.csv"),row.names=FALSE)

insGel=gettingInstrumentalDataFromListRepo(listRepoGel)
write.csv(insGel,file=paste0(saveRepo,"/insGel.csv"),row.names=FALSE)

#===========================================
# Merging instrumental data with PTR data
#===========================================

repoSave=saveRepo
repoSave="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf"
matrice="Gusto"
if(matrice=="Gel")
{
  ptr_gel=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/resultGel.xlsx",sheet="ptr")
  ins_gel=read.csv(paste0(repoSave,"/insGel.csv"))
  ptr=ptr_gel
  ins=ins_gel
  filename="conc_aci_gel.csv"
}
if(matrice=="Sol")
{
  ptr_sol=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/resultSol.xlsx",sheet="ptr")
  ins_sol=read.csv(paste0(repoSave,"/insSol.csv"))
  ptr=ptr_sol
  ins=rbind(ins_sol,ins_gus)
  ins=ins[substr(ins[,"file"],13,15)!="gus",]
  filename="conc_aci_sol.csv"
}
if(matrice=="Gusto")
{
  ptr_gus=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/resultGusto.xlsx",sheet="ptr")
  ins_gus=read.csv(paste0(repoSave,"/insGusto.csv"))
  ptr=ptr_gus
  ins=ins_gus
  filename="conc_aci_gus.csv"
}

# TODO integrate in vitro
sol_wide = ptr%>%
  group_by(time, file,ion,duration) %>%
  summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = ion,
    values_from = intensity
  )
length(unique(substr(sol_wide$file,8,12)))
sol_ins <- sol_wide %>%
  left_join(ins, by = c("file", "time"))
length(unique(substr(sol_ins$file,8,12)))

sol_ins=sol_ins[!is.na(sol_ins[,"ratio"]),]
unique(substr(sol_ins$file,8,12))
sol_ins[,"conc_ACI"]=sol_ins[,"isoamylAcetate"]*sol_ins[,"ratio"]/sol_ins[,"H3O+"]
sol_ins=as.data.frame(sol_ins)
write.csv(sol_ins,file=paste0(repoSave,"/",filename),row.names=FALSE)


