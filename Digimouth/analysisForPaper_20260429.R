# Functions to load
#=====================
getGroupCurves = function(df_summary_fop, groups, name_id="subject", 
                          ci = TRUE, ci_level = 0.95)
{
  df_summary_fop = as.data.frame(df_summary_fop)
  df_summary_fop[,"group"] = groups[paste0(df_summary_fop[,name_id])]
  
  df_summary_group <- df_summary_fop %>%
    group_by(time_bin, group) %>%
    summarise(
      mean_intensity2 = mean(mean_intensity, na.rm = TRUE),
      sd_intensity2   = sd(mean_intensity, na.rm = TRUE),
      n               = sum(!is.na(mean_intensity)),
      .groups = "drop"
    )
  
  # --- Calcul du ribbon ---
  if (ci) {
    # intervalle de confiance
    alpha = 1 - ci_level
    tval  = qt(1 - alpha/2, df_summary_group$n - 1)
    df_summary_group$lower = df_summary_group$mean_intensity2 - tval * df_summary_group$sd_intensity2 / sqrt(df_summary_group$n)
    df_summary_group$upper = df_summary_group$mean_intensity2 + tval * df_summary_group$sd_intensity2 / sqrt(df_summary_group$n)
  } else {
    # écart-type (ancienne version)
    df_summary_group$lower = df_summary_group$mean_intensity2 - df_summary_group$sd_intensity2
    df_summary_group$upper = df_summary_group$mean_intensity2 + df_summary_group$sd_intensity2
  }
  
  # --- plot ---
  p_fop = ggplot(df_summary_group, aes(x = time_bin, y = mean_intensity2, col = group)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.1) +
    labs(
      x = "Temps (s)",
      y = "average intensity",
      title = "Mean intensity with confidence intervals"
    ) +
    theme_minimal() +
   # geom_vline(xintercept = 15) +
    geom_vline(xintercept = 30, color = "red") +
    geom_vline(xintercept = 60, color = "blue")
  return(p_fop)
}
getIntensityCurves=function(arome_df,format="wide",fop,timePeriod=c(0,120),ion="conc_ACI",files=NULL,main_diff=NULL)
{
  if(!is.null(files)){arome_df=arome_df[arome_df[,"file"]%in%files,]}
  if(format=="long")
  {
    arome_df=arome_df[arome_df[,"ion"]==ion,]
  }
  if(format=="wide")
  {
    new_colnames=colnames(arome_df)
    new_colnames[new_colnames==ion]="intensity"
    colnames(arome_df)=new_colnames
  }
  
  arome=arome_df[arome_df[,"fop"]==fop[1]&arome_df[,"rep"]!=0,]
  print(arome[1,])
  df_summary_chew <- arome[arome[,"time"]>timePeriod[1]&arome[,"time"]<timePeriod[2],] %>%
    mutate(time_bin = round(time)) %>%       # arrondi au sec le plus proche
    group_by(time_bin) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      sd_intensity   = sd(intensity, na.rm = TRUE),
      .groups = "drop"
    )
  n_files_chew=length(unique(arome[,"file"]))
  p_fop1=ggplot(df_summary_chew, aes(x = time_bin, y = mean_intensity)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = mean_intensity - sd_intensity,
                    ymax = mean_intensity + sd_intensity),
                alpha = 0.2, fill = "blue") +
    labs(
      x = "Temps (s)",
      y = paste0("average intensity"),
      title = paste0("Mean intensity with standard deviations ",fop[2],"\n (n=",n_files_chew,")") )+
    theme_minimal()+
    geom_vline(xintercept=30,color="blue")+geom_vline(xintercept=60,color="red")
  
  arome=arome_df[arome_df[,"fop"]==fop[2]&arome_df[,"rep"]!=0,]
  df_summary_succ <- arome[arome[,"time"]>timePeriod[1]&arome[,"time"]<timePeriod[2],] %>%
    mutate(time_bin = round(time)) %>%       # arrondi au sec le plus proche
    group_by(time_bin) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      sd_intensity   = sd(intensity, na.rm = TRUE),
      .groups = "drop"
    )
  n_files_succ=length(unique(arome[,"file"]))
  p_fop2=ggplot(df_summary_succ, aes(x = time_bin, y = mean_intensity)) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = mean_intensity - sd_intensity,
                    ymax = mean_intensity + sd_intensity),
                alpha = 0.2, fill = "blue") +
    labs(
      x = "time (s)",
      y = paste0("mean intensity"),
      title = paste0("Mean intensity with standard deviations",fop[2],"\n (n=",n_files_succ,")")
    ) +
    theme_minimal()
   # geom_vline(xintercept=15,color="grey")+geom_vline(xintercept=30,color="magenta")+geom_vline(xintercept=60,color="red")
  
  # Comparison of fops
  df_summary_chew=as.data.frame(df_summary_chew)
  df_summary_succ=as.data.frame(df_summary_succ)
  df_summary_chew[,"fop"]=fop[1]
  df_summary_chew[,"std_err"]=df_summary_chew[,"sd_intensity"]/sqrt(n_files_chew)
  df_summary_succ[,"std_err"]=df_summary_succ[,"sd_intensity"]/sqrt(n_files_succ)
  df_summary_succ[,"fop"]=fop[2]
  df_summary_gel=rbind(df_summary_chew,df_summary_succ)
  p_diff=ggplot(df_summary_gel, aes(x = time_bin, y = mean_intensity,color=fop,fill=fop)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean_intensity - std_err,
                    ymax = mean_intensity + std_err),
                alpha = 0.2) +
    labs(
      x = "time (s)",
      y = paste0("mean intensity"),
      title = paste0("Mean intensity with standard errors \n ", "n_",fop[1],"=",n_files_chew,", n_",fop[2],"=",n_files_succ)
    ) +
    theme_minimal()+
    geom_vline(xintercept=30,color="red")+geom_vline(xintercept=60,color="blue")
  p_diff
  return(list(p_fop1=p_fop1,p_fop2=p_fop2,p_diff=p_diff))
}

plotOneFileFromPtr=function(ptr,format="wide",ion="conc_ACI",xlsxSenso,xlsxInt,subject,metaData,product,rep,fop,timeMEB=NULL,xlim=c(0,120),ylim=NULL)
{
  if(format=="long")
  {
    ptr=ptr[ptr[,"ion"]==ion,]
  }
  if(format=="wide")
  {
    new_colnames=colnames(ptr)
    new_colnames[new_colnames==ion]="intensity"
    colnames(ptr)=new_colnames
  }
  filePTR=metaData[!is.na(metaData[,"file"])&metaData[,"subject"]==subject&metaData[,"rep"]==rep&metaData["type"]==product&metaData[,"fop"]==fop,"file"]
  filePTR=filePTR[!is.na(filePTR)]
  timeMEB=timeMEB[!is.na(timeMEB)]
  
  ptrDataIon=ptr[ptr[,"file"]==paste0(filePTR,".h5"),]
  #ptrDataIon[,"time"]=ptrDataIon[,"time"]-timeMEB
  p=ggplot(ptrDataIon,aes(x=time,y=intensity))+geom_line()+theme_bw()
  dataSenso=xlsxSenso[xlsxSenso[,"product"]==product&xlsxSenso[,"rep"]==rep&xlsxSenso[,"fop"]==fop&xlsxSenso[,"subject"]==subject,]
  dataInt=xlsxInt[xlsxInt[,"product"]==product&xlsxInt[,"rep"]==rep&xlsxInt[,"fop"]==fop&xlsxInt[,"subject"]==subject,]
  int=round(dataInt[,"int"],digits=2)
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="chew","time"], linetype = "dashed", color = "red")
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="swallow","time"], linetype = "dotted", color = "blue")
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="swallowImposed","time"], color = "blue")
  p=p+geom_vline(xintercept = dataSenso[dataSenso[,"sw"]=="TotalSwallow","time"], color = "blue",linewidth=1.5)
  if(fop %in% c("chew","succ","fast","long")){p=p+geom_vline(xintercept=30, color = "red",linewidth=1.5)}
  p=p+ggtitle(paste0(ion,":",subject,", ",product,", rep ",rep, ",\n fop = ",fop, ",int=",int))+xlim(xlim[1],xlim[2])
  #if(!is.null(com)){print(com);if(length(com)>0){p=p+xlab(paste0("time (",com,")"))}}
  if(!is.null(ylim)){p=p+ylim(ylim[1],ylim[2])}
  return(p)
}
comparisonBoxplot=function(res_period,ion="isoamylAcetate")
{
  new_col=colnames(res_period)
  new_col[new_col==ion]="intensity"
  colnames(res_period)=new_col
  print(res_period[1,])
  moyennes_par_sujet_fop <- res_period %>%
    group_by(subject, fop) %>%
    summarise(moyenne_isoamyl = mean(intensity, na.rm = TRUE)) %>%
    ungroup()
  
  donnees_wide <- moyennes_par_sujet_fop %>%
    pivot_wider(names_from = fop, values_from = moyenne_isoamyl)
  
  # Faire le test de Wilcoxon pairé
  pval=wilcox.test(donnees_wide[[2]], donnees_wide[[3]], paired = TRUE)$p.value
  
  p_box=ggplot(moyennes_par_sujet_fop, aes(x = fop, y = moyenne_isoamyl)) +
    geom_boxplot(outlier.shape = NA, fill = "lightgray") +
    geom_line(aes(group = subject), color = "gray40", alpha = 0.6) +
    geom_point(aes(color = subject), size = 2) +
    theme_minimal() +
    labs(
      title = paste0("Paired comparison (p=",round(pval,digits=3),")"),
      x = "FOP",
      y = "Averaged AUC isoamylAcetate before swallow"
    ) +
    theme(legend.position = "none")
  return(p_box)
}

gettingIndicators=function(matrice="Solution",stat="a",reposave="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/csvIndicators/"
)
{
  if(matrice=="Solution")
  {
    res_period_not_a=read.csv(file=paste0(reposave,"area_sol_not.csv"),sep=";",header=T)
    res_period_mov_a=read.csv(file=paste0(reposave,"area_sol_mov.csv"),sep=";",header=T)
    res_period_deg_a=read.csv(file=paste0(reposave,"area_sol_deg.csv"),sep=";",header=T)
    res_period_tot_a=read.csv(file=paste0(reposave,"area_sol_tot.csv"),sep=";",header=T)
    res_period_not_t=read.csv(file=paste0(reposave,"tmax_sol_not.csv"),sep=";",header=T)
    res_period_mov_t=read.csv(file=paste0(reposave,"tmax_sol_mov.csv"),sep=";",header=T)
    res_period_deg_t=read.csv(file=paste0(reposave,"tmax_sol_deg.csv"),sep=";",header=T)
    res_period_tot_t=read.csv(file=paste0(reposave,"tmax_sol_tot.csv"),sep=";",header=T)
    res_period_not_m=read.csv(file=paste0(reposave,"imax_sol_not.csv"),sep=";",header=T)
    res_period_mov_m=read.csv(file=paste0(reposave,"imax_sol_mov.csv"),sep=";",header=T)
    res_period_deg_m=read.csv(file=paste0(reposave,"imax_sol_deg.csv"),sep=";",header=T)
    res_period_tot_m=read.csv(file=paste0(reposave,"imax_sol_tot.csv"),sep=";",header=T)
    fops=c("fast","long")
  }
  if(matrice=="Gusto")
  {
    res_period_not_a=read.csv( file=paste0(reposave,"area_gus_not.csv"),sep=";",header=T)
    res_period_deg_a=read.csv( file=paste0(reposave,"area_gus_deg.csv"),sep=";",header=T)
    res_period_tot_a=read.csv( file=paste0(reposave,"area_gus_tot.csv"),sep=";",header=T)
    res_period_not_t=read.csv( file=paste0(reposave,"tmax_gus_not.csv"),sep=";",header=T)
    res_period_deg_t=read.csv( file=paste0(reposave,"tmax_gus_deg.csv"),sep=";",header=T)
    res_period_tot_t=read.csv( file=paste0(reposave,"tmax_gus_tot.csv"),sep=";",header=T)
    res_period_not_m=read.csv( file=paste0(reposave,"imax_gus_not.csv"),sep=";",header=T)
    res_period_deg_m=read.csv( file=paste0(reposave,"imax_gus_deg.csv"),sep=";",header=T)
    res_period_tot_m=read.csv( file=paste0(reposave,"imax_gus_tot.csv"),sep=";",header=T)
    fops=c("cont","noth")
  }
  if(matrice=="Gel")
  {
    res_period_not_a=read.csv( file=paste0(reposave,"area_gel_not.csv"),sep=";",header=T)
    res_period_mov_a=read.csv( file=paste0(reposave,"area_gel_mov.csv"),sep=";",header=T)
    res_period_deg_a=read.csv( file=paste0(reposave,"area_gel_deg.csv"),sep=";",header=T)
    res_period_tot_a=read.csv( file=paste0(reposave,"area_gel_tot.csv"),sep=";",header=T)
    res_period_not_t=read.csv( file=paste0(reposave,"tmax_gel_not.csv"),sep=";",header=T)
    res_period_mov_t=read.csv( file=paste0(reposave,"tmax_gel_mov.csv"),sep=";",header=T)
    res_period_deg_t=read.csv( file=paste0(reposave,"tmax_gel_deg.csv"),sep=";",header=T)
    res_period_tot_t=read.csv( file=paste0(reposave,"tmax_gel_tot.csv"),sep=";",header=T)
    res_period_not_m=read.csv( file=paste0(reposave,"imax_gel_not.csv"),sep=";",header=T)
    res_period_mov_m=read.csv( file=paste0(reposave,"imax_gel_mov.csv"),sep=";",header=T)
    res_period_deg_m=read.csv( file=paste0(reposave,"imax_gel_deg.csv"),sep=";",header=T)
    res_period_tot_m=read.csv( file=paste0(reposave,"imax_gel_tot.csv"),sep=";",header=T)
    fops=c("chew","succ")
  }
  
  if(stat=="a")
  {
    df_to_use=list()
    df_to_use[["nothing"]]=res_period_not_a[res_period_not_a[,"rep"]!=0,]
    if(matrice!="Gusto")
    {
      df_to_use[["moving"]]=res_period_mov_a[res_period_mov_a[,"rep"]!=0,]
    }
    df_to_use[["swallow"]]=res_period_deg_a[res_period_deg_a[,"rep"]!=0,]
    df_to_use[["total"]]=res_period_tot_a[res_period_tot_a[,"rep"]!=0,]
  }
  if(stat=="m")
  {
    df_to_use=list()
    df_to_use[["nothing"]]=res_period_not_m[res_period_not_m[,"rep"]!=0,]
    if(matrice!="Gusto")
    {
      df_to_use[["moving"]]=res_period_mov_m[res_period_mov_m[,"rep"]!=0,]
    }
    df_to_use[["swallow"]]=res_period_deg_m[res_period_deg_m[,"rep"]!=0,]
    df_to_use[["total"]]=res_period_tot_m[res_period_tot_m[,"rep"]!=0,]
  }
  if(stat=="t")
  {
    df_to_use=list()
    df_to_use[["nothing"]]=res_period_not_t[res_period_not_t[,"rep"]!=0,]
    if(matrice!="Gusto")
    {
      df_to_use[["moving"]]=res_period_mov_t[res_period_mov_t[,"rep"]!=0,]
    }
    df_to_use[["swallow"]]=res_period_deg_t[res_period_deg_t[,"rep"]!=0,]
    df_to_use[["total"]]=res_period_tot_t[res_period_tot_t[,"rep"]!=0,]
  }
  
  summary(df_to_use[["swallow"]])
  df_to_use[["swallow"]][,"subject"]=as.factor(df_to_use[["swallow"]][,"subject"])
  df_to_use[["swallow"]][,"fop"]=as.factor(df_to_use[["swallow"]][,"fop"])
  df_to_use[["swallow"]]=unique(df_to_use[["swallow"]])
  df_to_use[["swallow"]][df_to_use[["swallow"]][,"subject"]=="P628",]
  summary(df_to_use[["swallow"]][,"subject"])
  return(df_to_use)
}

# Libraries and datasets to load
#=================================
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
repoFigure="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Valorisation/expe_paper/figures"

repoSave="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/"
ptr_sol=read.csv(file=paste0(repoSave,"conc_aci_sol.csv"))
metadata_sol=read.csv(paste0(repoSave,"/metadata_sol.csv"),sep=";",header=T)
ptr_gus=read.csv(file=paste0(repoSave,"conc_aci_gus.csv"))
metadata_gus=read.csv(paste0(repoSave,"/metadata_gus.csv"),sep=";",header=T)
senso_sol=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/senso_sol.xlsx")
int_sol=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/int_sol.xlsx")
senso_gus=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/senso_gus.xlsx")
int_gus=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/int_gus.xlsx")
int_gus=int_gus[!duplicated(int_gus),]

ptr_gel=read.csv(file=paste0(repoSave,"conc_aci_gel.csv"))
metadata_gel=read.csv(paste0(repoSave,"/metadata_gel.csv"),sep=";",header=T)
senso_gel=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/senso_gel.xlsx")
int_gel=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/int_gel.xlsx")

reposave="P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/csvIndicators/"
# Merging data with metadata
arome_sol=merge(ptr_sol[,c("file","duration","time","conc_ACI")],metadata_sol,by.x="file",by.y="file_h5")
arome_gus=merge(ptr_gus[,c("file","duration","time","conc_ACI")],metadata_gus,by.x="file",by.y="file_h5")
arome_gel=merge(ptr_gel[,c("file","duration","time","conc_ACI")],metadata_gel,by.x="file",by.y="file_h5")

arome_sol[arome_sol[,"fop"]=="long","fop"]="rare"
arome_sol[arome_sol[,"fop"]=="fast","fop"]="freq"
unique(arome_sol[,"subject"])[!unique(arome_sol[,"subject"])%in%unique(arome_gus[,"subject"])]
unique(arome_sol[,"subject"])[!unique(arome_sol[,"subject"])%in%unique(arome_gel[,"subject"])]

# Analysis of one subject
#=========================
# Gel analysis
subject="Q407" #Q407 T453 Z089
q407_chew1=plotOneFileFromPtr(ptr=ptr_gel,xlsxSenso=senso_gel,subject=subject,xlsxInt=int_gel,metaData=metadata_gel,product="Gel",rep=1,fop="chew",xlim=c(0,120),ylim=NULL)
q407_chew2=plotOneFileFromPtr(ptr=ptr_gel,xlsxSenso=senso_gel,subject=subject,xlsxInt=int_gel,metaData=metadata_gel,product="Gel",rep=2,fop="chew",xlim=c(0,120),ylim=NULL)
q407_chew3=plotOneFileFromPtr(ptr=ptr_gel,xlsxSenso=senso_gel,subject=subject,xlsxInt=int_gel,metaData=metadata_gel,product="Gel",rep=3,fop="chew",xlim=c(0,120),ylim=NULL)
q407_succ1=plotOneFileFromPtr(ptr=ptr_gel,xlsxSenso=senso_gel,subject=subject,xlsxInt=int_gel,metaData=metadata_gel,product="Gel",rep=1,fop="succ",xlim=c(0,120),ylim=NULL)
q407_succ2=plotOneFileFromPtr(ptr=ptr_gel,xlsxSenso=senso_gel,subject=subject,xlsxInt=int_gel,metaData=metadata_gel,product="Gel",rep=2,fop="succ",xlim=c(0,120),ylim=NULL)
q407_succ3=plotOneFileFromPtr(ptr=ptr_gel,xlsxSenso=senso_gel,subject=subject,xlsxInt=int_gel,metaData=metadata_gel,product="Gel",rep=3,fop="succ",xlim=c(0,120),ylim=NULL)
#grid.arrange(q407_chew1,q407_succ1,q407_chew2,q407_succ2,q407_chew3,
#q407_succ3)

q407_long1=plotOneFileFromPtr(ptr=ptr_sol,xlsxSenso=senso_sol,subject=subject,xlsxInt=int_sol,metaData=metadata_sol,product="Sol",rep=1,fop="long",xlim=c(0,120),ylim=NULL)
q407_long2=plotOneFileFromPtr(ptr=ptr_sol,xlsxSenso=senso_sol,subject=subject,xlsxInt=int_sol,metaData=metadata_sol,product="Sol",rep=2,fop="long",xlim=c(0,120),ylim=NULL)
q407_long3=plotOneFileFromPtr(ptr=ptr_sol,xlsxSenso=senso_sol,subject=subject,xlsxInt=int_sol,metaData=metadata_sol,product="Sol",rep=3,fop="long",xlim=c(0,120),ylim=NULL)
q407_fast1=plotOneFileFromPtr(ptr=ptr_sol,xlsxSenso=senso_sol,subject=subject,xlsxInt=int_sol,metaData=metadata_sol,product="Sol",rep=1,fop="fast",xlim=c(0,120),ylim=NULL)
q407_fast2=plotOneFileFromPtr(ptr=ptr_sol,xlsxSenso=senso_sol,subject=subject,xlsxInt=int_sol,metaData=metadata_sol,product="Sol",rep=2,fop="fast",xlim=c(0,120),ylim=NULL)
q407_fast3=plotOneFileFromPtr(ptr=ptr_sol,xlsxSenso=senso_sol,subject=subject,xlsxInt=int_sol,metaData=metadata_sol,product="Sol",rep=3,fop="fast",xlim=c(0,120),ylim=NULL)
#grid.arrange(q407_long1,q407_fast1,q407_long2,q407_fast2,q407_long3,
#             q407_fast3)

q407_noth1=plotOneFileFromPtr(ptr=ptr_gus,xlsxSenso=senso_gus,subject=subject,xlsxInt=int_gus,metaData=metadata_gus,product="Gus",rep=1,fop="noth",xlim=c(0,120),ylim=NULL)
q407_noth2=plotOneFileFromPtr(ptr=ptr_gus,xlsxSenso=senso_gus,subject=subject,xlsxInt=int_gus,metaData=metadata_gus,product="Gus",rep=2,fop="noth",xlim=c(0,120),ylim=NULL)
q407_noth3=plotOneFileFromPtr(ptr=ptr_gus,xlsxSenso=senso_gus,subject=subject,xlsxInt=int_gus,metaData=metadata_gus,product="Gus",rep=3,fop="noth",xlim=c(0,120),ylim=NULL)
q407_cont1=plotOneFileFromPtr(ptr=ptr_gus,xlsxSenso=senso_gus,subject=subject,xlsxInt=int_gus,metaData=metadata_gus,product="Gus",rep=1,fop="cont",xlim=c(0,120),ylim=NULL)
q407_cont2=plotOneFileFromPtr(ptr=ptr_gus,xlsxSenso=senso_gus,subject=subject,xlsxInt=int_gus,metaData=metadata_gus,product="Gus",rep=2,fop="cont",xlim=c(0,120),ylim=NULL)
q407_cont3=plotOneFileFromPtr(ptr=ptr_gus,xlsxSenso=senso_gus,subject=subject,xlsxInt=int_gus,metaData=metadata_gus,product="Gus",rep=3,fop="cont",xlim=c(0,120),ylim=NULL)
#grid.arrange(q407_noth1,q407_cont1,q407_noth2,q407_cont2,q407_noth3,
#             q407_cont3)
size_text=10
size_axis=9
size_axis_2=9
m1=5;m2=5
m3=m4=5
#==================
# Figure 2
#=================
pdf(paste0(repoFigure,"/Figure_2.pdf"),width=10,height=10)
grid.arrange(q407_fast1+ggtitle("a. Solution in a cup (FOP='freq')")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("Aroma concentration (VPM)")+scale_x_continuous(breaks = seq(0,120, by = 15))+ylim(0,65000)+xlim(0,120),
             q407_long1+ggtitle("b. Solution in a cup (FOP='rare')")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("Aroma concentration (VPM)")+scale_x_continuous(breaks = seq(0,120, by = 15))+ylim(0,65000)+xlim(0,120),
             q407_noth1+ggtitle("c. Solution delivered by gustometer (FOP='noth')")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("Aroma concentration (VPM)")+scale_x_continuous(breaks = seq(0,120, by = 15))+ylim(0,81000)+xlim(0,120),
             q407_cont1+ggtitle("d. Solution delivered by gustometer (FOP='cont')")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("Aroma concentration (VPM)")+scale_x_continuous(breaks = seq(0,120, by = 15))+ylim(0,81000)+xlim(0,120),
             q407_succ1+ggtitle("e. Gelatin gummy disc (FOP='succ')")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("Aroma concentration (VPM)")+scale_x_continuous(breaks = seq(0,120, by = 15))+ylim(0,810)+xlim(0,120),
             q407_chew1+ggtitle("f. Gelatin gummy disc (FOP='chew')")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("Aroma concentration (VPM)")+scale_x_continuous(breaks = seq(0,120, by = 15))+ylim(0,810)+xlim(0,120)
)
dev.off()

arome_sol[arome_sol[,"fop"]=="fast","fop"]="freq"
arome_sol[arome_sol[,"fop"]=="long","fop"]="rare"
res_sol=getIntensityCurves(arome_df=arome_sol,fop=c("freq","rare"))
res_gus=getIntensityCurves(arome_df=arome_gus,fop=c("noth","cont"))
res_gel=getIntensityCurves(arome_df=arome_gel,fop=c("succ","chew"))

# Breathing
ptr_sol[,"isoprene"]
ptr_gel[,"isoprene"]
ptr_gus[,"isoprene"]

# Figure 3 Averaged concentration curves
#============================
pdf(paste0(repoFigure,"/Figure_3.pdf"),width=10,height=10)
grid.arrange(res_sol$p_diff+ ggtitle("Solution in a cup (n=90)")+theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("mean aroma concentration")+scale_x_continuous(breaks = seq(0,120, by = 15)),
             res_gus$p_diff+ ggtitle("Solution delivered by gustometer (n=87)")+theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("mean aroma concentration"),
             res_gel$p_diff+ geom_vline(xintercept=15,col="grey")+ ggtitle("Gelatin gummy discs (n=87)")+theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+ylab("mean aroma concentration"),
                                   nrow=3)
dev.off()

# Variability
indic_sol_a=gettingIndicators(matrice="Solution",stat="a")
indic_sol_m=gettingIndicators(matrice="Solution",stat="m")
indic_gel_a=gettingIndicators(matrice="Gel",stat="a")
indic_gel_m=gettingIndicators(matrice="Gel",stat="m")
indic_gus_a=gettingIndicators(matrice="Gusto",stat="a")
indic_gus_m=gettingIndicators(matrice="Gusto",stat="m")

indic=indic_gus_a[["total"]] # decline over fop and phases
fop="cont"
df_indicator=indic[indic[,"fop"]==fop,]
df_indicator_avg <- df_indicator %>% 
  group_by(subject) %>%
  summarise(
    mean_g = mean(conc_ACI, na.rm = TRUE),
    sd_intra   = sd(conc_ACI, na.rm = TRUE),
    .groups = "drop"
  )

m=mean(df_indicator[,"conc_ACI"])
n=nrow(df_indicator_avg)
sd_total=sd(df_indicator[,"conc_ACI"])
sd_inter=sd(df_indicator_avg$mean_g)*sqrt(n-1)/sqrt(n)
sd_intra=mean(df_indicator_avg$sd_intra)*sqrt(2)/sqrt(3)
print("m=")
print(m)
print("sd=")
print(sd_total);
print("sd inter =")
print(sd_inter);
print("sd intra=")
print(sd_intra)


# Area under curves
#====================
df_to_use_sol=gettingIndicators(matrice="Solution",stat="a")
df_to_use_gel=gettingIndicators(matrice="Gel",stat="a")
df_to_use_gus=gettingIndicators(matrice="Gusto",stat="a")

long_df=df_to_use_sol$total[df_to_use_sol$total[,"fop"]=="long" ,]
fast_df=df_to_use_sol$total[df_to_use_sol$total[,"fop"]=="fast" ,]
long_df_o=long_df[order(long_df[,"subject"],long_df[,"rep"]),]
fast_df_o=fast_df[order(fast_df[,"subject"],fast_df[,"rep"]),]
sum(!long_df_o[,"subject"]==fast_df_o["subject"])
sum(!long_df_o[,"rep"]==fast_df_o["rep"])
wilcox.test(long_df_o[,"conc_ACI"],fast_df_o[,"conc_ACI"],paired=T)
         
ion="conc_ACI"
#ion="isoprene"
p=list()
periods=c("nothing","moving","swallow","total")
#periods=c("nothing","swallow","total")
k=1
fops=c("chew","long")
for(fop in fops)
{
  for(i in periods)
  {
    dfi=df_to_use[[i]][df_to_use[[i]][,"fop"]==fop, c("file", "conc_ACI", "subject")]
    p[[k]]=ggplot(dfi,
                  aes(x = subject, y = conc_ACI, color = subject)) +
      geom_boxplot() +
      geom_jitter(height = 0) +theme_bw()+
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )+ggtitle(paste0(i,", ",fop))
    k=k+1
  }
}

grid.arrange(grobs=p,nrow=2)
grid.arrange(grobs=p[c(1,2,5,6)],nrow=2)
grid.arrange(grobs=p[c(3,4,7,8)],nrow=2)

# comparison of conditions
#==============
# Figure 4
#==============
df_to_use_sol[["swallow"]][,"fop"]=as.character(df_to_use_sol[["swallow"]][,"fop"])
df_to_use_sol[["swallow"]][df_to_use_sol[["swallow"]][,"fop"]=="fast","fop"]="freq"
df_to_use_sol[["swallow"]][!is.na(df_to_use_sol[["swallow"]][,"fop"])&df_to_use_sol[["swallow"]][,"fop"]=="long","fop"]="rare"
p_sol=comparisonBoxplot(res_period=df_to_use_sol[["swallow"]],ion="conc_ACI")
p_gus=comparisonBoxplot(res_period=df_to_use_gus[["nothing"]],ion="conc_ACI")
p_gel=comparisonBoxplot(res_period=df_to_use_gel[["moving"]],ion="conc_ACI")
pdf(paste0(repoFigure,"/Figure_4.pdf"),width=10,height=6)
grid.arrange(p_sol+ggtitle("a. AUC of aroma curves \n(solution,60-120s)")+ylab("AUC"),
             p_gus+ggtitle("b. AUC of aroma curves \n(gustometer, 0-60s),")+ylab("AUC"),
             p_gel+ggtitle("c. AUC of aroma curves\n (gelatin gummy discs,30-60s)")+ylab("AUC")
             ,nrow=1)
dev.off()
#===========
# Figure 6
#===========

p_sol_int=comparisonBoxplot(res_period=int_sol,ion="int")+ylab("Declared intensity")
p_gus_int=comparisonBoxplot(res_period=int_gus,ion="int")+ylab("Declared intensity")
p_gel_int=comparisonBoxplot(res_period=int_gel,ion="int")+ylab("Declared intensity")

pdf(paste0(repoFigure,"/Figure_6.pdf"),width=10,height=6)
grid.arrange(p_sol_int+ggtitle("Paired comparison \n(Solution,pval=0.915)"),
             p_gus_int+ggtitle("Paired comparison \n(Gusto, pval=0.949)"),
             p_gel_int+ggtitle("Paired comparison \n(Gelatin gummy, pval=0.36)"),nrow=1)
dev.off()


auc_sol=df_to_use_sol[["total"]]
auc_gel=df_to_use_gel[["total"]]
auc_gus=df_to_use_gus[["total"]]
int_sol_o=int_sol[order(int_sol[,"fop"],int_sol[,"subject"],int_sol[,"rep"]),]
auc_sol_o=auc_sol[order(auc_sol[,"fop"],auc_sol[,"subject"],auc_sol[,"rep"]),]
int_gel_o=int_gel[order(int_gel[,"fop"],int_gel[,"subject"],int_gel[,"rep"]),]
auc_gel_o=auc_gel[order(auc_gel[,"fop"],auc_gel[,"subject"],auc_gel[,"rep"]),]
int_gus_o=int_gus[order(int_gus[,"fop"],int_gus[,"subject"],int_gus[,"rep"]),]
auc_gus_o=auc_gus[order(auc_gus[,"fop"],auc_gus[,"subject"],auc_gus[,"rep"]),]




auc_tot=rbind(auc_sol_o,auc_gus_o,auc_gel_o)
int_tot=rbind(int_sol,int_gus,int_gel)
auc_tot[,"id"]=paste0(auc_tot[,"subject"], auc_tot[,"fop"],auc_tot[,"rep"])
int_tot[,"id"]=paste0(int_tot[,"subject"], int_tot[,"fop"],int_tot[,"rep"])
auc_tot=auc_tot[auc_tot[,"id"]!=c("V154noth2","V154noth4"),]
rownames(auc_tot)=auc_tot[,"id"]
rownames(int_tot)=int_tot[,"id"]
subjects=unique(int_tot[,"subject"])

cor_subj=rep(NA,length(subjects));names(cor_subj)=subjects
cor_subj_gel=cor_subj_gus=cor_subj_sol=p_subj=p_subj_gel=p_subj_sol=p_subj_gus=cor_subj
subjects
for( subject in subjects)
{ print(subject)
  subj_auc_data=auc_tot[auc_tot[,"subject"]==subject,]
  cor_subj[subject]=cor(subj_auc_data[,"conc_ACI"],int_tot[rownames(subj_auc_data),"int"],method="spearman",use="pairwise.complete.obs")
  #p_subj=cor.test(subj_auc_data[,"conc_ACI"],int_tot[rownames(subj_auc_data),"int"],method="spearman")$p.value
  # Gel

  subj_auc_data_gel=auc_tot[auc_tot[,"subject"]==subject&auc_tot[,"fop"]%in%c("chew","succ"),]
    if(dim(subj_auc_data_gel)[1]!=0)
    {
    cor_subj_gel[subject]=cor(subj_auc_data_gel[,"conc_ACI"],int_tot[rownames(subj_auc_data_gel),"int"],method="spearman",use="pairwise.complete.obs")
  }    # Gus
  subj_auc_data_gus=auc_tot[auc_tot[,"subject"]==subject&auc_tot[,"fop"]%in%c("cont","noth"),]
  if(dim(subj_auc_data_gus)[1]!=0)
    cor_subj_gus[subject]=cor(subj_auc_data_gus[,"conc_ACI"],int_tot[rownames(subj_auc_data_gus),"int"],method="spearman",use="pairwise.complete.obs")
  #sol
  subj_auc_data_sol=auc_tot[auc_tot[,"subject"]==subject&auc_tot[,"fop"]%in%c("fast","long"),]
  cor_subj_sol[subject]=cor(subj_auc_data_sol[,"conc_ACI"],int_tot[rownames(subj_auc_data_sol),"int"],method="spearman",use="pairwise.complete.obs")
}
mean((cor_subj),na.rm=T)
mean((cor_subj_sol),na.rm=T)
mean((cor_subj_gus),na.rm=T)
mean((cor_subj_gel),na.rm=T)

cor(auc_tot[,"conc_ACI"],int_tot[rownames(auc_tot),"int"], method="spearman",use="pairwise.complete.obs")
auc_gel_o[,"id"]=paste0(auc_gel_o[,"subject"], auc_gel_o[,"fop"],auc_gel_o[,"rep"])
int_gel_o[,"id"]=paste0(int_gel_o[,"subject"], int_gel_o[,"fop"],int_gel_o[,"rep"])
rownames(auc_gel_o)=auc_gel_o[,"id"]
rownames(int_gel_o)=int_gel_o[,"id"]
all(auc_gel_o[,"subject"]==int_gel_o[,"subject"])
cor(auc_gel_o[,"conc_ACI"],int_gel_o[rownames(auc_gel_o),"int"],method="spearman")
all(auc_gus_o[,"subject"]==int_gus_o[,"subject"])
all(auc_sol_o[,"subject"]==int_sol_o[,"subject"])
cor(auc_gus_o[,"conc_ACI"],int_gus_o[,"int"],method="spearman")
cor(auc_sol_o[,"conc_ACI"],int_sol_o[,"int"],method="spearman")


# Clustering
#====================

# Same method on AUC
#=====================
# Classification on AUC
# calculate auc = sum of intensities for constant steptime
arome_df=arome_sol
fop=c("rare","freq")
timePeriod=c(30,60)
#timePeriod=c(0,30)
arome=arome_df[arome_df[,"fop"]%in%fop&arome_df[,"rep"]!=0,]
df_summary_fop <- arome[arome[,"time"]>timePeriod[1]&arome[,"time"]<timePeriod[2],] %>%
  mutate(time_bin = round(time)) %>%       # arrondi au sec le plus proche
  group_by(time_bin,fop,subject,rep) %>%
  summarise(
    mean_intensity = mean(log(conc_ACI+1,base=10), na.rm = TRUE),
    sd_intensity   = sd(log(conc_ACI+1,base=10), na.rm = TRUE),
    .groups = "drop"
  )

df_summary_fop_auc=df_summary_fop%>%
  group_by(rep,fop,subject) %>%
  summarise(
    auc = mean(mean_intensity, na.rm = TRUE),
    .groups = "drop"
  )

df_avg_auc=df_summary_fop_auc%>%
  group_by(subject) %>%
  summarise(
    avg_auc = mean(auc, na.rm = TRUE),
    .groups = "drop"
  )

hist(df_summary_fop_auc$auc,breaks=20)
hist(df_avg_auc$avg_auc,breaks=10)
plot(df_avg_auc$avg_auc,df_avg_auc$avg_auc)

m_aci=df_avg_auc$avg_auc
names(m_aci)=as.data.frame(df_avg_auc)[,"subject"]
m_dist=dist(m_aci,method="euclidean");
resh=hclust(m_dist,method="ward.D2")
p_dend=plot(resh)
groups_auc=paste0("G",cutree(resh,k=2)); 
names(groups_auc)=names(m_aci)

# group curves
groups_auc
summary(factor(groups_auc))

df_summary_subject <- arome[arome[,"time"]>timePeriod[1]&arome[,"time"]<timePeriod[2],] %>%
  mutate(time_bin = round(time)) %>%       # arrondi au sec le plus proche
  group_by(subject,time_bin) %>%
  summarise(
    mean_intensity = mean(log(conc_ACI+1,base=10), na.rm = TRUE),
    sd_intensity   = sd(log(conc_ACI+1,base=10), na.rm = TRUE),
    .groups = "drop"
  )

df_summary_subject[,"group"]=groups_auc[as.data.frame(df_summary_subject)[,"subject"]]
p_subj=ggplot(df_summary_subject, aes(x = time_bin, y = mean_intensity,col=group,group=subject)) +
  geom_line(aes(color = group))+
  geom_ribbon(aes(ymin = mean_intensity - sd_intensity,
                  ymax = mean_intensity + sd_intensity,fill=group),
              alpha = 0.05,color = NA) +
  labs(
    x = "Time (s)",
    y = paste0("log-averaged intensity"),
    title = paste0("Mean log concentrations with standard deviations ",fop[2]) )+
  theme_minimal()+
  geom_vline(xintercept=15,color="grey")+geom_vline(xintercept=30,color="blue")+geom_vline(xintercept=60,color="red")
p_subj+xlim(timePeriod[1],timePeriod[2])+ggtitle("a. Heterogeneity of ACI release and clustering")

p_auc=getGroupCurves(df_summary_subject, groups=groups_auc)
p_auc+xlim(timePeriod[1],timePeriod[2])

# Raw data
df_summary_eval <- arome[arome[,"time"]>timePeriod[1]&arome[,"time"]<timePeriod[2],] 
df_summary_eval[,"group"]=groups_auc[as.data.frame(df_summary_eval)[,"subject"]]
df_summary_eval[,"conc_ACI_log"]=log(df_summary_eval[,"conc_ACI"]+1,base=10)
p_eval=ggplot(df_summary_eval, aes(x = time, y = conc_ACI,col=group,group=file)) +
  geom_line(aes(color = group),linewidth=0.1,alpha=0.1)+
  labs(
    x = "Time (s)",
    y = paste0("log- intensity"),
    title = paste0("Mean intensity with standard deviations ",fop[2]) )+
  theme_minimal()+
  geom_vline(xintercept=15,color="grey")+geom_vline(xintercept=30,color="red")+geom_vline(xintercept=60,color="blue")
p_eval+xlim(timePeriod[1],timePeriod[2])


df_summary_eval[,"group"]=groups_auc[as.data.frame(df_summary_eval)[,"subject"]]
df_summary_eval[,"conc_ACI_log"]=log(df_summary_eval[,"conc_ACI"]+1,base=10)

files_G1=unique(df_summary_eval[df_summary_eval[,"subject"]=="Q407","file"])
files_G2=unique(df_summary_eval[df_summary_eval[,"subject"]=="T386","file"])

files_examples_G1=(files_G1)
files_examples_G2=(files_G2)
p_eval_G1=ggplot(df_summary_eval[df_summary_eval[,"group"]=="G1"&df_summary_eval[,"file"]%in%files_examples_G1,], aes(x = time, y = conc_ACI,col=file,group=file)) +
  geom_line(aes(color = file),linewidth=1)+
  labs(
    x = "Time (s)",
    y = paste0("intensity"),
    title = paste0("Raw data for Q407" ))+
  theme_minimal()+
  geom_vline(xintercept=15,color="grey")+geom_vline(xintercept=30,color="red")+geom_vline(xintercept=60,color="blue")

p_eval_G2=ggplot(df_summary_eval[df_summary_eval[,"group"]=="G2"&df_summary_eval[,"file"]%in%files_examples_G2,], aes(x = time, y = conc_ACI,col=file,group=file)) +
  geom_line(aes(color = file),linewidth=1)+
  labs(
    x = "Time (s)",
    y = paste0("intensity"),
    title = paste0("Raw data for T387"))+
  theme_minimal()+
  geom_vline(xintercept=15,color="grey")+geom_vline(xintercept=30,color="red")+geom_vline(xintercept=60,color="blue")
p_G1=p_eval_G1+xlim(timePeriod[1],timePeriod[2])+ylim(0,200000)+theme(legend.position="none")
p_G2=p_eval_G2+xlim(timePeriod[1],timePeriod[2])+ylim(0,200000)+theme(legend.position="none")

# group eval sur juste après la déglutition
df_summary_eval <- arome[arome[,"time"]>timePeriod[2]&arome[,"time"]<75,] 
df_summary_eval[,"group"]=groups_auc[as.data.frame(df_summary_eval)[,"subject"]]
df_summary_eval[,"conc_ACI_log"]=log(df_summary_eval[,"conc_ACI"]+1,base=10)

p_eval2=getGroupCurves(df_summary_eval, groups=groups_auc)
p_eval2+xlim(timePeriod[1],75)

# A cheval sur les deux periodes
# group curves
period=c(0,120)
df_summary_subject <- arome[arome[,"time"]>period[1]&arome[,"time"]<period[2],] %>%
  mutate(time_bin = round(time)) %>%       # arrondi au sec le plus proche
  group_by(subject,time_bin) %>%
  summarise(
    mean_intensity = mean(log(conc_ACI+1,base=10), na.rm = TRUE),
    sd_intensity   = sd(log(conc_ACI+1,base=10), na.rm = TRUE),
    .groups = "drop"
  )
df_summary_subject[,"group"]=groups_auc[as.data.frame(df_summary_subject)[,"subject"]]
p_subj_2=ggplot(df_summary_subject, aes(x = time_bin, y = mean_intensity,col=group,group=subject)) +
  geom_line(aes(color = group))+
  geom_ribbon(aes(ymin = mean_intensity - sd_intensity,
                  ymax = mean_intensity + sd_intensity,fill=group),
              alpha = 0.05,color = NA) +
  labs(
    x = "Time (s)",
    y = paste0("log-averaged intensity"),
    title = paste0("Mean intensity with standard deviations ",fop[2]) )+
  theme_minimal()+
  geom_vline(xintercept=30,color="red")+geom_vline(xintercept=60,color="blue")
p_subj_2+xlim(period[1],period[2])

df_summary_subject=as.data.frame(df_summary_subject)
p_eval2=getGroupCurves(df_summary_subject, groups=groups_auc)
p_eval2

size_text=10
size_axis=9
size_axis_2=9
m1=4;m2=5
m3=m4=4


# Figure 7
#=================

int_sol=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/int_sol.xlsx")
g1=names(groups_auc[groups_auc=="G1"])
g2=names(groups_auc[groups_auc=="G2"])

senso_sol=read.xlsx("P:/Chemosens/Projets/Projets_En_Cours/1-Projets-Aromes/P00354_Digimouth/Donnees/AcquisitionAvrilSeptembre2025/finalDf/senso_sol.xlsx")
result_sol <- senso_sol %>%
  filter(rep != 0) %>% 
  filter(time >30 & time<60) %>% 
  #  filter(fop == "chew") %>% 
  filter(sw == "chew") %>% 
  group_by(subject, product, rep, fop, sw) %>%
  summarise(count = n(), .groups = "drop")
result_sol=as.data.frame(result_sol)
result_sol[,"group"]=groups_auc[result_sol[,"subject"]]
wilcox.test(result_sol[result_sol[,"group"]=="G1","count"],result_sol[result_sol[,"group"]=="G2","count"])
t.test(result_sol[result_sol[,"group"]=="G1","count"],result_sol[result_sol[,"group"]=="G2","count"])

p_chew=ggplot(result_sol,aes(x=group,y=count,color=group))+geom_boxplot()+geom_jitter(height=0)+theme_bw()+  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+ggtitle("Number of chews (gel)")
p_chew
pdf(paste0(repoFigure,"/Figure 7.pdf"),width=8,height=8)
grid.arrange(
  p_G1 +ggtitle("a. Example of raw data for a subject in G1")+ylab("Aroma concentration")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4)),
  p_G2+ggtitle("b. Example of raw data for a subject in G2")+ylab("Aroma concentration")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4)),
  p_eval2+ggtitle("c. Aroma release kinetics according to clustering results")+ylab("Aroma log-concentration.")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))+xlab("Time (s)"),
  p_chew+ggtitle("d. Number of chews (solution in a cup)")+ theme(plot.title = element_text(size = size_text), axis.title = element_text(size = size_axis),axis.text  = element_text(size = size_axis_2),plot.margin = margin(m1, m2, m3, m4))
)
dev.off()