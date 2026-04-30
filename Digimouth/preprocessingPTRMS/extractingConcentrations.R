library(tidyr)

ptr_gel=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/resultGel.xlsx",sheet="ptr")
ptr_gus=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/resultGusto.xlsx",sheet="ptr")
ptr_sol=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/resultSol.xlsx",sheet="ptr")
repoSave="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf"
ins_sol=read.csv(paste0(repoSave,"/insSol.csv"))
ins_gus=read.csv(paste0(repoSave,"/insGusto.csv"))
ins_gel=read.csv(paste0(repoSave,"/insGel.csv"))

ptr=ptr_gel
ins=ins_gel
filename="conc_aci_gel.csv"
sol_wide = ptr%>%
  group_by(time, file,ion) %>%
  summarise(intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = ion,
    values_from = intensity
)

sol_ins <- sol_wide %>%
  left_join(ins, by = c("file", "time"))

sol_ins=sol_ins[!is.na(sol_ins[,"ratio"]),]
summary(sol_ins[,"file"])
sol_ins[,"conc_ACI"]=sol_ins[,"isoamylAcetate"]*sol_ins[,"ratio"]/sol_ins[,"H3O+"]
sol_ins=as.data.frame(sol_ins)
write.csv(sol_ins,file=paste0(repoSave,"/",filename),row.names=FALSE)

