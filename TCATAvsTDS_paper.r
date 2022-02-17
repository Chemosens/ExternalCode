# Analysis Cantin data
#=========================
library(ggplot2)
library(CSUtils)
library(MSnbase)
library(reshape2)
library(PTRMSR)
library(gridExtra)
#====================
# Pretreatment (Figure 1)
#==================
#wd="C:/INRA/Data/Donnees Cantin/Cantin-DTS-PTRviewer"
#listFilesTest=list.files()
#setwd(wd)
datasetTest=read.table(file="CSCA098_S002_2.txt",sep='\t',header=T,dec=",")
ion_to_useTest=unique(colnames(datasetTest)[201:397])
referenceBreathLong="m69.06906..69.06906...Conc."
ion_to_useTest=ion_to_useTest[!ion_to_useTest%in%referenceBreathLong]

datasetTest2=datasetTest[,c("AbsTime",   "RelTime",   "Cycle",ion_to_useTest,referenceBreathLong  )]
colnames(datasetTest2)=substr(colnames(datasetTest2),1,7)
report=ptrvReport(dataset=datasetTest2,selecIons="namely",
                  listIons=substr(c("m71.08535..71.08535...Conc.","m115.11308..115.11...Conc." ),1,7),
                  referenceBreath=substr("m69.06906..69.06906...Conc.",1,7),
                  methodDetectStart="startPeakProportion",
                  noisePeriodIBT=c(0,25),noisePeriodSig=c(0,25),noisePeriodDS=c(0,25),
                  proportionOfMax=0.3,halfWindowSize=5,maxPeaks=NULL,
                  startPeriod=c(20,60),detectingStart=30)

p1=report$gg$p_breath$p_cyclelimits
p2=report$gg$p_curves$p_raw
p3=report$gg$p_curves$p_cycle


datasetTestB=read.table(file="CSCA100_S001_3.txt",sep='\t',header=T,dec=".")
metaData_tds[metaData_tds[,"file"]%in%listFilesTest[c(1,7)],]
datasetTest2B=datasetTestB[,c("AbsTime",   "RelTime",   "Cycle",ion_to_useTest,referenceBreathLong  )]
colnames(datasetTest2B)=substr(colnames(datasetTest2B),1,7)
reportB=ptrvReport(dataset=datasetTest2B,selecIons="namely",
                  listIons=substr(c("m71.08535..71.08535...Conc.","m115.11308..115.11...Conc." ),1,7),
                  referenceBreath=substr("m69.06906..69.06906...Conc.",1,7),
                  methodDetectStart="startPeakProportion",
                  noisePeriodIBT=c(0,25),noisePeriodSig=c(0,25),noisePeriodDS=c(0,25),
                  proportionOfMax=0.3,halfWindowSize=5,maxPeaks=NULL,
                  startPeriod=c(20,60),detectingStart=30)

p1B=reportB$gg$p_breath$p_cyclelimits
p2B=reportB$gg$p_curves$p_raw
p3B=reportB$gg$p_curves$p_cycle


grid.arrange(p1+ggtitle("a."),
             p1B+ggtitle("d."),
             p2+theme(legend.position="none") + ggtitle("b.")+theme(legend.position="bottom"),
             p2B+theme(legend.position="none") + ggtitle("e.")+theme(legend.position="bottom"),
             p3+ ggtitle("c. ")+theme(legend.position="bottom"),
             p3B+ ggtitle("f.")+theme(legend.position="bottom")
             )


# Relevant ion selection
#======================
#setwd("C:/INRA/Data/Donnees Cantin/Donnees&SeanceTimeSens")
dec_vec_tds=rep(".",96);dec_vec_tds[1]=",";
noisePeriod=c(0,25)
halfWindowSize=5

# Selecting only relevant ions for each evaluations
setwd("C:/INRA/Data/Donnees Cantin/Cantin-DTS-PTRviewer")
listFilesTDS=list.files(pattern="*.txt")
metaData_tds=read.table("metaData_tds.csv",sep=";",header=T)
ion_tds=ptrvListSignificantSNRIons(listFilesTDS,metaData=metaData_tds,dec_vec=dec_vec_tds,
                                   noisePeriod=noisePeriod,halfWindowSize=halfWindowSize)

setwd("C:/INRA/Data/Donnees Cantin/Cantin-TCATA-PTRviewer")
listFilesTCATA=list.files(pattern="*.txt")
metaData_tcata=read.table("metaData_tcata.csv",sep=";",header=T)
ion_tcata=ptrvListSignificantSNRIons(listFilesTCATA,metaData=metaData_tcata,
                                     noisePeriod=noisePeriod,halfWindowSize=halfWindowSize)

ions_total=names(ion_tcata$resIons[[1]])
ions_total2=ions_total[sapply(ions_total,function(x)return(substr(x,nchar(x)-5,nchar(x))))==".Conc."]

matSigTds=Reduce(rbind,ion_tds$resIons)
matSigTcata=Reduce(rbind,ion_tcata$resIons)
matSig=rbind(matSigTds,matSigTcata)
vec_prod=c(ion_tds$product,ion_tcata$product)

resA=apply(matSig[vec_prod=="A",],2,function(x){return(sum(x>3,na.rm=T))})
resB=apply(matSig[vec_prod=="B",],2,function(x){return(sum(x>3,na.rm=T))})
resC=apply(matSig[vec_prod=="C",],2,function(x){return(sum(x>3,na.rm=T))})
numberOfEvalByProduct=dim(matSig[vec_prod=="A",])[1]
minimalNumberOfSelection=1*numberOfEvalByProduct
ionA=names(resA[resA>=minimalNumberOfSelection])
ionB=names(resB[resB>=minimalNumberOfSelection])
ionC=names(resC[resC>=minimalNumberOfSelection])
ionSigUnique=unique(c(ionA,ionB,ionC))
listeIonsIndesirables=c("m39.02314..39.02314...Conc.",
                        "m69.06906..69.06906...Conc.")
ionSigUnique2=ionSigUnique[!ionSigUnique%in%listeIonsIndesirables]
# Removing doubles
ion_to_use=ionSigUnique2[sapply(ionSigUnique2,function(x)return(substr(x,nchar(x)-5,nchar(x))))==".Conc."]
ion_to_use

#==========================
# Pretreatment of data
#==========================
# wd="C:/INRA/Data/Donnees Cantin/Cantin-DTS-PTRviewer"
# setwd(wd)
# Pretreatemnt of TDS and TCATA data
tds=tdsRead(file="TDS3.csv", cols=list(subject="SubjectCode", product="ProductCode", descriptor="Descriptor",time="Time",score="Score",rep="Replicate"), supCols="", sep=";",startWithFirstCitation=FALSE,discretization=1,periods=1)
tcata=tcataRead(file="TCATA3.csv", cols=list(subject="SubjectCode", product="ProductCode", descriptor="Descriptor",time="Time",score="Score",rep="Replicate"), supCols="", sep=";",startWithFirstCitation=FALSE,discretization=1,periods=1)

# Pretreatment of PTR-TDS data
listFilesTDS=list.files(pattern="*.txt")
metaData_tds=read.table("metaData_tds2.csv",sep=";",header=T)
res_intensity_tds=ptrvListIntensityByTime(listFilesTDS,removeBlank=TRUE,dec_vec=dec_vec_tds,metaData=metaData_tds,ions=ion_to_use,correction = "cycle",halfWindowSize=halfWindowSize,maxPeaks=NULL,method="SuperSmoother",total=FALSE,breathRatio=FALSE,timeBlank = noisePeriod)
res_intensity_tds$totalIntensity[,"ion"]=substr(res_intensity_tds$totalIntensity$ion,1,7)
res_intensity_tds$time[,"ion"]=substr(res_intensity_tds$time$ion,1,7)
colnames(res_intensity_tds$listRes)[-c(1:4)]=substr(colnames(res_intensity_tds$listRes)[-c(1:4)],1,7)
bloc_auc_tds=res_intensity_tds$listRes[,-c(1:4)];rownames(bloc_auc_tds)=paste0(res_intensity_tds$listRes[,"subject"],"_",res_intensity_tds$listRes[,"product"],"_",res_intensity_tds$listRes[,"rep"])
bloc_auc_tds=bloc_auc_tds[,-which(colnames(bloc_auc_tds)=="m69.069")]
rownamesToUse=rownames(bloc_auc_tds)
bloc_auc_tds_wo=bloc_auc_tds[rownamesToUse,]
bloc_auc_tdsN=100*sweep(bloc_auc_tds,1,apply(bloc_auc_tds,1,sum),"/")

# Pretreatment of PTR-TCATA data
# wd="C:/INRA/Data/Donnees Cantin/Cantin-TCATA-PTRviewer"
# setwd(wd)
listFilesTCATA=list.files(pattern="*.txt")
metaData_tcata=read.table("metaData_tcata2.csv",sep=";",header=T)
res_intensity_tcata=ptrvListIntensityByTime(listFilesTCATA,removeBlank=TRUE,dec_vec=rep(".",length(listFilesTCATA)),metaData=metaData_tcata,ions=ion_to_use,correction = "cycle",halfWindowSize=halfWindowSize,maxPeaks=NULL,method="SuperSmoother",total=FALSE,breathRatio=FALSE,timeBlank = noisePeriod)
res_intensity_tcata$totalIntensity[,"ion"]=substr(res_intensity_tcata$totalIntensity$ion,1,7)
res_intensity_tcata$time[,"ion"]=substr(res_intensity_tcata$time$ion,1,7)
colnames(res_intensity_tcata$listRes)[-c(1:4)]=substr(colnames(res_intensity_tcata$listRes)[-c(1:4)],1,7)

bloc_auc_tcata=res_intensity_tcata$listRes[,-c(1:4)];rownames(bloc_auc_tcata)=paste0(res_intensity_tcata$listRes[,"subject"],"_",res_intensity_tcata$listRes[,"product"],"_",res_intensity_tcata$listRes[,"rep"])
bloc_auc_tcata=bloc_auc_tcata[,-which(colnames(bloc_auc_tcata)=="m69.069")]
bloc_auc_tcata=bloc_auc_tcata[rownamesToUse,]
bloc_auc_tcataN=100*sweep(bloc_auc_tcata,1,apply(bloc_auc_tcata,1,sum),"/")

df_durations_tds=dcast(tds$durations,subject+product+rep~descriptor,fun.aggregate=sum,value.var="score")
df_durations_tds=df_durations_tds[df_durations_tds[,"product"]!="WarmUp",]

bloc_durations_tds=df_durations_tds[,-c(1:3)];rownames(bloc_durations_tds)=paste0(df_durations_tds[,"subject"],"_",df_durations_tds[,"product"],"_",df_durations_tds[,"rep"])
bloc_durations_tds_wo=bloc_durations_tds[rownamesToUse,]
bloc_durations_tdsN=100*sweep(bloc_durations_tds,1,apply(bloc_durations_tds,1,sum),"/")

df_durations_tcata=dcast(tcata$durations,subject+product+rep~descriptor,fun.aggregate=sum,value.var="score")
df_durations_tcata=df_durations_tcata[df_durations_tcata[,"product"]!="WarmUp",]

bloc_durations_tcata=df_durations_tcata[,-c(1:3)];rownames(bloc_durations_tcata)=paste0(df_durations_tcata[,"subject"],"_",df_durations_tcata[,"product"],"_",df_durations_tcata[,"rep"])
bloc_durations_tcata_wo=bloc_durations_tcata[rownamesToUse,]
bloc_durations_tcataN=100*sweep(bloc_durations_tcata,1,apply(bloc_durations_tcata,1,sum),"/")


# Correction by total durations
behav_tcata=dcast(tcata$behaviours, subject+product+rep~variable,value.var="score",fun.aggregate=mean)
behav_tcata2=behav_tcata[behav_tcata[,"product"]!="WarmUp",c("subject","product","rep","sequenceDuration")]
df_durations_tcata_corr=merge(behav_tcata2,df_durations_tcata,by=c("product","subject","rep"))
bloc_durations_tcata_corr=df_durations_tcata_corr[,-c(1:4)];
bloc_durations_tcata_corr=sweep(bloc_durations_tcata_corr,1,df_durations_tcata_corr[,"sequenceDuration"],"/")
rownames(bloc_durations_tcata_corr)=paste0(df_durations_tcata_corr[,"subject"],"_",df_durations_tcata_corr[,"product"],"_",df_durations_tcata_corr[,"rep"])
bloc_durations_tcata_corr=bloc_durations_tcata_corr[rownamesToUse,]

behav_tds=dcast(tds$behaviours, subject+product+rep~variable,value.var="score",fun.aggregate=mean)
behav_tds2=behav_tds[behav_tds[,"product"]!="WarmUp",c("subject","product","rep","sequenceDuration")]
df_durations_tds_corr=merge(behav_tds2,df_durations_tds,by=c("product","subject","rep"))
bloc_durations_tds_corr=df_durations_tds_corr[,-c(1:4)];
bloc_durations_tds_corr=sweep(bloc_durations_tds_corr,1,df_durations_tds_corr[,"sequenceDuration"],"/")
rownames(bloc_durations_tds_corr)=paste0(df_durations_tds_corr[,"subject"],"_",df_durations_tds_corr[,"product"],"_",df_durations_tds_corr[,"rep"])
bloc_durations_tds_corr=bloc_durations_tds_corr[rownamesToUse,]


rownames(bloc_auc_tds)==rownames(bloc_auc_tcata)
rownames(bloc_auc_tds)==rownames(bloc_durations_tcata_corr)
rownames(bloc_durations_tds_corr)==rownames(bloc_durations_tcata_corr)


# Getting data in long format (ggplot type)
auc_tcata_final=reshape(bloc_auc_tcataN,varying=list(1:dim(bloc_auc_tcataN)[2]),times=colnames(bloc_auc_tcataN),direction="long",v.names="intensity",timevar="ion",ids=rownames(bloc_auc_tcataN))
auc_tcata_final[,"product"]=substr(auc_tcata_final[,"id"],6,6)
auc_tcata_final[,"subject"]=substr(auc_tcata_final[,"id"],1,4)
auc_tcata_final[,"rep"]=substr(auc_tcata_final[,"id"],8,8)

auc_tds_final=reshape(bloc_auc_tdsN,varying=list(1:dim(bloc_auc_tdsN)[2]),times=colnames(bloc_auc_tdsN),direction="long",v.names="intensity",timevar="ion",ids=rownames(bloc_auc_tdsN))
auc_tds_final[,"product"]=substr(auc_tds_final[,"id"],6,6)
auc_tds_final[,"subject"]=substr(auc_tds_final[,"id"],1,4)
auc_tds_final[,"rep"]=substr(auc_tds_final[,"id"],8,8)

tds_final=reshape(bloc_durations_tds_corr,varying=list(1:dim(bloc_durations_tds_corr)[2]),times=colnames(bloc_durations_tds_corr),direction="long",v.names="intensity",timevar="descriptor",ids=rownames(bloc_durations_tds_corr))
tds_final[,"product"]=substr(tds_final[,"id"],6,6)
tds_final[,"subject"]=substr(tds_final[,"id"],1,4)
tds_final[,"rep"]=substr(tds_final[,"id"],8,8)

tcata_final=reshape(bloc_durations_tcata_corr,varying=list(1:dim(bloc_durations_tds_corr)[2]),times=colnames(bloc_durations_tcata_corr),direction="long",v.names="intensity",timevar="descriptor",ids=rownames(bloc_durations_tcata_corr))
tcata_final[,"product"]=substr(tcata_final[,"id"],6,6)
tcata_final[,"subject"]=substr(tcata_final[,"id"],1,4)
tcata_final[,"rep"]=substr(tcata_final[,"id"],8,8)

# Using RGCCA 
#==========================================
library(RGCCA)
connection=matrix(0,4,4)
connection[1,2]=connection[2,1]=connection[3,4]=connection[4,3]=1
resPLS=rgcca(list(ptrtds=bloc_auc_tdsN,tds=bloc_durations_tds_corr,ptrtcata=bloc_auc_tcataN,tcata=bloc_durations_tcata_corr),tau=1,connection=connection,ncomp=2,scheme="horst")
p_mpls1=plot(resPLS,block=1,type="sample",resp=(substr(rownames(bloc_auc_tds),6,6)),colors=c(  "#71ad65" ,  "#cd5b45",  "darkgoldenrod3"),text_ind=F,title="PTR-TDS")
p_mpls2=plot(resPLS,block=2,type="sample",resp=(substr(rownames(bloc_auc_tds),6,6)),colors=c(  "#71ad65" ,  "#cd5b45",  "darkgoldenrod3"),text_ind=F,title="TDS")
p_mpls3=plot(resPLS,block=3,type="sample",resp=(substr(rownames(bloc_auc_tds),6,6)),colors=c(  "#71ad65" ,  "#cd5b45",  "darkgoldenrod3"),text_ind=F,title="PTR-TCATA")
p_mpls4=plot(resPLS,block=4,type="sample",resp=(substr(rownames(bloc_auc_tds),6,6)),colors=c(  "#71ad65" ,  "#cd5b45",  "darkgoldenrod3"),text_ind=F,title="TCATA")

ave_ptrtds1=100*round(resPLS$AVE$AVE_X$ptrtds[1],digits=4)
ave_ptrtds2=100*round(resPLS$AVE$AVE_X$ptrtds[1],digits=4)
ave_ptrtcata1=100*round(resPLS$AVE$AVE_X$ptrtcata[1],digits=4)
ave_ptrtcata2=100*round(resPLS$AVE$AVE_X$ptrtcata[2],digits=4)
ave_tds1=100*round(resPLS$AVE$AVE_X$tds[1],digits=4)
ave_tds2=100*round(resPLS$AVE$AVE_X$tds[1],digits=4)
ave_tcata1=100*round(resPLS$AVE$AVE_X$tcata[1],digits=4)
ave_tcata2=100*round(resPLS$AVE$AVE_X$tcata[2],digits=4)

F01_tds=summary(manova(resPLS$Y[[1]][,1:2]~substr(rownames(resPLS$Y[[1]]),6,6)))$stats[1,"approx F"]
F02_tds=summary(manova(resPLS$Y[[2]][,1:2]~substr(rownames(resPLS$Y[[2]]),6,6)))$stats[1,"approx F"]
F01_tcata=summary(manova(resPLS$Y[[3]][,1:2]~substr(rownames(resPLS$Y[[3]]),6,6)))$stats[1,"approx F"]
F02_tcata=summary(manova(resPLS$Y[[4]][,1:2]~substr(rownames(resPLS$Y[[4]]),6,6)))$stats[1,"approx F"]

p01=ggplot(data.frame(comp1=resPLS$Y[[1]][,1],comp2=resPLS$Y[[1]][,2],product=(substr(rownames(resPLS$Y[[1]]),6,6)),name=rownames(resPLS$Y[[1]])),aes(x=comp1,y=comp2,col=product,name=name))+geom_point()+theme_bw()+ggtitle(paste0("PTR-TDS (F=",round(F01_tds,digits=2),")"))+theme(legend.position="none")+xlab(paste0("Axis 1 (",ave_ptrtds1,"%)"))+ylab(paste0("Axis 2 (",ave_ptrtds2,"%)"))
p02=ggplot(data.frame(comp1=resPLS$Y[[2]][,1],comp2=resPLS$Y[[2]][,2],product=(substr(rownames(resPLS$Y[[2]]),6,6)),name=rownames(resPLS$Y[[2]])),aes(x=comp1,y=comp2,col=product,name=name))+geom_point()+theme_bw()+ggtitle(paste0("TDS (F=",round(F02_tds,digits=2),")"))+theme(legend.position="none")+xlab(paste0("Axis 1 (",ave_tds1,"%)"))+ylab(paste0("Axis 2 (",ave_tds2,"%)"))
p03=ggplot(data.frame(comp1=resPLS$Y[[3]][,1],comp2=resPLS$Y[[3]][,2],product=(substr(rownames(resPLS$Y[[3]]),6,6)),name=rownames(resPLS$Y[[3]])),aes(x=comp1,y=comp2,col=product,name=name))+geom_point()+theme_bw()+ggtitle(paste0("PTR-TCATA (F=",round(F01_tcata,digits=2),")"))+theme(legend.position="none")+xlab(paste0("Axis 1 (",ave_ptrtcata1,"%)"))+ylab(paste0("Axis 2 (",ave_ptrtcata2,"%)"))
p04=ggplot(data.frame(comp1=resPLS$Y[[4]][,1],comp2=resPLS$Y[[4]][,2],product=(substr(rownames(resPLS$Y[[4]]),6,6)),name=rownames(resPLS$Y[[4]])),aes(x=comp1,y=comp2,col=product,name=name))+geom_point()+theme_bw()+ggtitle(paste0("TCATA (F=",round(F02_tcata,digits=2),")"))+theme(legend.position="none")+xlab(paste0("Axis 1 (",ave_tcata1,"%)"))+ylab(paste0("Axis 2 (",ave_tcata2,"%)"))


grid.arrange(p01,p02,p03,p04)
grid.arrange(p01+ggtitle("a."),p02+ggtitle("c."),p03+ggtitle("b."),p04+ggtitle("d."))

# correlations
round(cor(cbind(ptr_tds=resPLS$Y[[1]][,1],tds=resPLS$Y[[2]][,1],ptr_tcata=resPLS$Y[[3]][,1],tcata=resPLS$Y[[4]][,1])),digits=2)
round(cor(cbind(ptr_tds=resPLS$Y[[1]][,2],tds=resPLS$Y[[2]][,2],ptr_tcata=resPLS$Y[[3]][,2],tcata=resPLS$Y[[4]][,2])),digits=2)
library(FactoMineR)
coeffRV(resPLS$Y[[1]][,1:2],resPLS$Y[[2]][,1:2])$rv
coeffRV(resPLS$Y[[3]][,1:2],resPLS$Y[[4]][,1:2])$rv
coeffRV(resPLS$Y[[2]][,1:2],resPLS$Y[[4]][,1:2])$rv
coeffRV(resPLS$Y[[1]][,1:2],resPLS$Y[[3]][,1:2])$rv

# Clusters (table 1)
reshclust_ptrtds=hclust(dist(resPLS$Y[[1]][,1:2]),method="ward.D2")
res_classif_ptrtds=cutree(reshclust_ptrtds,k=3)
summary(factor(substr(names(res_classif_ptrtds[res_classif_ptrtds==1]),6,6)))
summary(factor(substr(names(res_classif_ptrtds[res_classif_ptrtds==2]),6,6)))
summary(factor(substr(names(res_classif_ptrtds[res_classif_ptrtds==3]),6,6)))

reshclust_tds=hclust(dist(resPLS$Y[[2]][,1:2]),method="ward.D2")
res_classif_tds=cutree(reshclust_tds,k=3)
summary(factor(substr(names(res_classif_tds[res_classif_tds==1]),6,6)))
summary(factor(substr(names(res_classif_tds[res_classif_tds==2]),6,6)))
summary(factor(substr(names(res_classif_tds[res_classif_tds==3]),6,6)))

reshclust_ptrtcata=hclust(dist(resPLS$Y[[3]][,1:2]),method="ward.D2")
res_classif_ptrtcata=cutree(reshclust_ptrtcata,k=3)
summary(factor(substr(names(res_classif_ptrtcata[res_classif_ptrtcata==1]),6,6)))
summary(factor(substr(names(res_classif_ptrtcata[res_classif_ptrtcata==2]),6,6)))
summary(factor(substr(names(res_classif_ptrtcata[res_classif_ptrtcata==3]),6,6)))

reshclust_tcata=hclust(dist(resPLS$Y[[4]][,1:2]),method="ward.D2")
res_classif_tcata=cutree(reshclust_tcata,k=3)
summary(factor(substr(names(res_classif_tcata[res_classif_tcata==1]),6,6)))
summary(factor(substr(names(res_classif_tcata[res_classif_tcata==2]),6,6)))
summary(factor(substr(names(res_classif_tcata[res_classif_tcata==3]),6,6)))

# Getting bootstrap results
#============================
bootres=bootstrap(resPLS,n_boot=1000)
nmark=25
p_boot1=plot(bootres,block=1,n_mark=nmark)
p_boot2=plot(bootres,block=2,n_mark=nmark)
p_boot3=plot(bootres,block=3,n_mark=nmark)
p_boot4=plot(bootres,block=4,n_mark=nmark)
grid.arrange(p_boot1+ggtitle("PTR-TDS: Axis 1 \n (1000 runs)"),p_boot2+ggtitle("TDS: Axis 1 \n (1000 runs)"),p_boot3 +ggtitle("PTR-TCATA: Axis 1 \n (1000 runs)"),p_boot4+ggtitle("TCATA: Axis 1 \n (1000 runs)"),nrow=2,ncol=2)
p_boot1_2=plot(bootres,block=1,comp=2,n_mark=nmark)
p_boot2_2=plot(bootres,block=2,comp=2,n_mark=nmark)
p_boot3_2=plot(bootres,block=3,comp=2,n_mark=nmark)
p_boot4_2=plot(bootres,block=4,comp=2,n_mark=nmark)
grid.arrange( p_boot1_2+ggtitle("PTR-TDS: Axis 2 \n (1000 runs)"),p_boot2_2+ggtitle("TDS: Axis 2 \n (1000 runs)"),p_boot3_2 +ggtitle("PTR-TCATA: Axis 2 \n (1000 runs)"),p_boot4_2+ggtitle("TCATA: Axis 2 \n (1000 runs)"),     nrow=2,ncol=2)



# Supplementary analyses
#=============================

# Getting sensory curves (Figure 4)
library(chemosensR)
analysis(tds,type="Standardized dominance curves",selection='product!="WarmUp"')
analysis(tcata,type="Standardized dominance curves",selection='product!="WarmUp"')
p_1=analysis(tds,type="Standardized dominance curves",selection='product!="WarmUp"')
p_2=analysis(tcata,type="Standardized dominance curves",selection='product!="WarmUp"',draw="")
grid.arrange(p_1[[1]]$output$dominanceCurves+ggtitle("TDS curves"),p_2[[1]]$output$dominanceCurves+ggtitle("TCATA curves"),ncol=2)

# Getting anova results (Table 2 and 3)

anova_tcata=anovaTable( df=tcata_final,model= "intensity~product+(1|subject)+(1|product:subject)",  alpha = 0.05,  unit = "descriptor",
                        columns = c("G_Mean", "F_Product", "P_Product", "F_Product_Significance",
                                    "SE_Product", "R2", "P_ShapiroWilks", "P_Levene_Product", "Mean_Product",
                                    "Group_Product"))

anova_tds=anovaTable( df=tds_final,model= "intensity~product+(1|subject)+(1|product:subject)",  alpha = 0.05,  unit = "descriptor",
                      columns = c("G_Mean", "F_Product", "P_Product", "F_Product_Significance",
                                  "SE_Product", "R2", "P_ShapiroWilks", "P_Levene_Product", "Mean_Product",
                                  "Group_Product"))

antds=anova_tds[,c("descriptor","F_Product","F_Product_Significance","Group_A","Group_B","Group_C")];colnames(antds)=c("descriptor","TDS F","TDS sig.","A_tds","B_tds","C_tds")
antcata=anova_tcata[,c("descriptor","F_Product","F_Product_Significance","Group_A","Group_B","Group_C")];colnames(antcata)=c("descriptor","TCATA F","TCATA sig.","A_tcata","B_tcata","C_tcata")

resanovasenso=merge(antcata,antds,by.x="descriptor",by.y="descriptor")
write.table(resanovasenso,"anovaTablesenso.csv",row.names=F,sep=";")


anova_ptrtcata=anovaTable( df=auc_tcata_final,model= "intensity~product+(1|subject)+(1|product:subject)",  alpha = 0.05,  unit = "ion",
                           columns = c("G_Mean", "F_Product", "P_Product", "F_Product_Significance",
                                       "SE_Product", "R2", "P_ShapiroWilks", "P_Levene_Product", "Mean_Product",
                                       "Group_Product"))

anova_ptrtds=anovaTable( df=auc_tds_final,model= "intensity~product+(1|subject)+(1|product:subject)",  alpha = 0.05,  unit = "ion",
                         columns = c("G_Mean", "F_Product", "P_Product", "F_Product_Significance",
                                     "SE_Product", "R2", "P_ShapiroWilks", "P_Levene_Product", "Mean_Product",
                                     "Group_Product"))

anptrtcata=anova_ptrtcata[,c("ion","F_Product","F_Product_Significance","Group_A","Group_B","Group_C")];colnames(anptrtcata)=c("ion","TCATA F","TCATA sig.","A_tcata","B_tcata","C_tcata")
anptrtds=anova_ptrtds[,c("ion","F_Product","F_Product_Significance","Group_A","Group_B","Group_C")];colnames(anptrtds)=c("ion","TDS F","TDS sig.","A_tds","B_tds","C_tds")
resanovaptr=merge(anptrtcata,anptrtds,by.x="ion",by.y="ion")
#write.table(resanovaptr,"anovaTableptr.csv",row.names=F,sep=";")

