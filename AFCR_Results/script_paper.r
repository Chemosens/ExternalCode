library(devtools)
library(openxlsx)
library(plotly)
library(gridExtra)
library(openxlsx)
library(MASS)
library(AFCR)
# To run only for updating the package
#remove.packages("AFCR")
#install_github("https://github.com/ChemoSens/AFCR")

geom_avg=function(conc)
{
  ga=rep(NA,length(conc)-1)
  for(i in 2:length(conc))
  {
    ga[i-1]=sqrt(conc[i-1]*conc[i])
  }
  return(ga)
}
get_groupBET=function(dfT,nb_total=193)
{
  bet_group=1
  for(i in 1:length(dfT[,"t"]))
  {
    nb_i=dfT[i,"t"]*nb_total
    bet_group=bet_group*dfT[i,"BET"]^nb_i
  }
  bet_group=bet_group^(1/nb_total)
  return(bet_group)
}


# Simulations: application 1, Figure 1
#===============
par(mfrow=c(2,2))
pS=rep(0,7);pS[4]=1
barplot(200*pS,ylim=c(0,200),names=paste0("S_c",0:6,"-S",1:7), main="Theoretical distribution for Tr=3",col="red")
matTassumingS=AFCR:::getMatrixTsachantS(2/3,6)
pT_estimated=pS%*%matTassumingS
pT=as.vector(pT_estimated)
dfT=data.frame(t=pT,names=paste0("]c",0:6,";c",1:7,"]"))
dfS=data.frame(t=pS,names=paste0("]c",0:6,";c",1:7,""))
real=ggplot(dfS,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("a. Distribution of the true threshold")+theme_bw()+ylab("Probability")+xlab("Threshold")
obs=ggplot(dfT,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("b. Distribution of BET threshold")+theme_bw()+xlab("Threshold")+ylab("Probability")
grid.arrange(real,obs,nrow=1)
barplot(200*pT_estimated,ylim=c(0,200),names=paste0("T_c",0:6,"-c",1:7), main="Theoretical results of 6 successive 3-AFC tests \n assuming Tr=3 for everyone",col="light blue")

# Real applications
#=======================
#setwd("C:/Users/capeltier/Nextcloud/MyDrive/DataAnalysis/2023_ThèseAlix/AFC")
triangular=read.xlsx("SEANCE 22 09 2023.xlsx",sheet="TDS")
triangular=triangular[,c("Répétition","Panéliste","Produit","Descripteur","Score","Temps")]
colnames(triangular)=c("replicate","subject","product","descriptor","score","time")
res_klo=keepLastOccurence(triangular,subjectName="subject",productName="product",descriptorName="descriptor",timeName="time")
threshold_data=getThreshold(res=res_klo,rata=NULL,
                                 decreasingConcentrations=c("C8","C7","C6","C5","C4","C3","C2","C1"))
head(threshold_data)

# Plotting BET distribution
distributionThreshold=summary(as.factor(threshold_data[,"threshold"]))
probaT=distributionThreshold/sum(distributionThreshold)
pT=as.vector(probaT)
dfT=data.frame(t=pT,names=paste0("]c",0:8,";c",1:9,"]"))
obs=ggplot(dfT,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("a. Observed BET threshold distribution")+theme_bw()+xlab("Threshold")+ylab("Probability")

# Calculating corrected BET distribution
conditional_proba=AFCR:::getMatrixTsachantS(2/3,8)
round(conditional_proba,digits=2)
probaS=t(ginv(conditional_proba))%*%probaT
summary(sum(probaS)-1) # Is a probability
probaS=as.vector(probaS)
names(probaS)=paste0("]c_", 0:8,";c_",1:9,"]")

# Plotting corrected BET distribution
pS=probaS
dfS=data.frame(t=pS,names=paste0("]c",0:8,";c",1:9,"]"))
real=ggplot(dfS,aes(x=names,y=t))+geom_col()+ylim(0,1)+theme_bw()+ylab("Probability")+xlab("Threshold")+ggtitle("b. True threshold distribution after correction")
grid.arrange(obs+ylim(0,0.35),real+ylim(0,0.35),nrow=1)

# comparing group threshold
decreasingNumConcentrations0=c(1.001,0.550,0.302,0.166,0.091,0.05,0.027,0.015)
pas=decreasingNumConcentrations0[1]/decreasingNumConcentrations0[2]
minConc=decreasingNumConcentrations0[length(decreasingNumConcentrations0)]/pas
maxConc=decreasingNumConcentrations0[1]*pas
concentrations=c(minConc,rev(decreasingNumConcentrations0),maxConc)
BETvalues=geom_avg(concentrations)
dfT=cbind(dfT,BET=BETvalues)
dfS=cbind(dfS,BET=BETvalues)
geom_avg(c(0.0015,1.001))
get_groupBET(dfT,nb_total=193)
get_groupBET(dfS,nb_total=193)

# Other threshold calculations from literature
afcMatrix=getTriangularMatrix(triangular, subjectName="subject", productName="product",descriptorName="descriptor",timeName="time", increasingConcentrations=paste0("C",1:8))
p_correct=apply(afcMatrix,2,sum)/nrow(afcMatrix)
p_discr=(p_correct-1/3)/(2/3)
plot(p_discr,pch=16)

# Getting lawless results
res_lawless=lawlessCorrection(numericConcentrations=rev(exp(decreasingNumConcentrations)),
                              afcMatrix)
exp(res_lawless) 

# Getting hough correction
houghTriangular=houghCorrection(afcMatrix,p_chance=1/3)

p_dis=apply(houghTriangular,2,sum)/nrow(houghTriangular)
plot(rev(decreasingNumConcentrations),p_dis, pch=16,xlab="Log-concentrations",ylab="Proportion of discriminators",main="Hough correction")
coeffLm=lm(p_dis~rev(decreasingNumConcentrations))$coef
abline(a=coeffLm[1],b=coeffLm[2])
abline(h=0.5)
logConc=(0.5-coeffLm[1])/coeffLm[2]
exp(logConc)
summary(factor(seuils0_ini[,"threshold"]))
pd_our=cumsum(summary(factor(seuils0_ini[,"threshold"])))/dim(seuils0_ini)[1]
plot( sort(unique(seuils0_ini[,"thresholdNum"])),pd_our, pch=16)



# Autre exemple sur TOM
#============================
library(openxlsx)
setwd("C:/Users/capeltier/Nextcloud/MyDrive/DataAnalysis/2023_ThèseAlix/AFC")
tom0=read.xlsx("Donnée_sensibilite_TOM.xlsx",sheet="Score sensibilité (0 à 6)")
# Remarque: dans tom, les seuils sont de 0 si aucun test réussi, 6 si tous les tests sont réussis....
# C'est l'exact inverse de ce qu'on a considéré
# Du coup, je transforme les données pour obtenir la même chose que nous:
tom=tom0
tom[,-1]=6-tom0[,-1]
K=6
saveur="amer"
saveur="sucre"
if(saveur=="sucre"){
  resTom=tom[,2:4];resTomM=tom[,2:4]-apply(tom[,2:4],1,median)}
if(saveur=="sale"){resTom=tom[,5:7];resTomM=tom[,5:7]-apply(tom[,5:7],1,median)}
if(saveur=="acide"){resTom=tom[,8:10];resTomM=tom[,8:10]-apply(tom[,8:10],1,median)}
if(saveur=="amer"){resTom=tom[,11:13];resTomM=tom[,11:13]-apply(tom[,11:13],1,median)}
if(saveur=="umami"){resTom=tom[,14:16];resTomM=tom[,14:16]-apply(tom[,14:16],1,median)}
resTomTotalCol=cbind(tom[,2:4,drop=F],tom[,5:7],tom[,8:10],tom[,11:13],tom[,14:16])
resTomTotalCol

names(tom)[2:4] <- names(tom)[5:7] <- names(tom)[8:10] <- names(tom)[11:13] <- names(tom)[14:16]
resTomTotal=rbind(tom[,2:4,drop=F],tom[,5:7],tom[,8:10],tom[,11:13],tom[,14:16])

distributionSeuilSucre=summary(factor(resTomTotalCol[,"Sucre-1"]))
distributionSeuilSale=summary(factor(resTomTotalCol[,"Sale-1"]))
distributionSeuilAcide=summary(factor(resTomTotalCol[,"Acide-1"]))
distributionSeuilAmer=summary(factor(resTomTotalCol[,"Amer1"]))
distributionSeuilUmami=summary(factor(resTomTotalCol[,"Umami1"]))
decreasingConcentrationsSucre=c(35.4,16.1,7.1,4.1,2.7,2.1)
decreasingConcentrationsSale=c(12.6,8.8,6.7,5.6,5.2,4.9)
decreasingConcentrationsAcide=c(24.7,17.1,9.5,5.6,4.1,3.6)
decreasingConcentrationsAmer=c(133.5,98.3,48.8,19.2,9.5,2.9)
decreasingConcentrationsUmami=c(29.6,19.6,3.7,9,1.6,1.1)

decreasingConcentrationsSucre[-length(decreasingConcentrationsSucre)]/decreasingConcentrationsSucre[-1]
decreasingConcentrationsSale[-length(decreasingConcentrationsSale)]/decreasingConcentrationsSale[-1]
decreasingConcentrationsAcide[-length(decreasingConcentrationsAcide)]/decreasingConcentrationsAcide[-1]
decreasingConcentrationsAmer[-length(decreasingConcentrationsAmer)]/decreasingConcentrationsAmer[-1]
decreasingConcentrationsUmami[-length(decreasingConcentrationsUmami)]/decreasingConcentrationsUmami[-1]

distributionSeuils=distributionSeuilUmami # change the name for each descriptor
decreasingNumConcentrations0=decreasingConcentrationsUmami

probaT=distributionSeuils/sum(distributionSeuils)
probaS=t(ginv(AFCR:::getMatrixTsachantS(2/3,6)))%*%probaT
probaS=as.vector(probaS)
names(probaS)=paste0("S", 0:6)
probaS*100

pas=decreasingNumConcentrations0[1]/decreasingNumConcentrations0[2]
minConc=decreasingNumConcentrations0[length(decreasingNumConcentrations0)]/pas
maxConc=decreasingNumConcentrations0[1]*pas
concentrations=c(minConc,rev(decreasingNumConcentrations0),maxConc)
BETvalues=geom_avg(concentrations)
dfT=cbind(t=distributionSeuils/100,BET=BETvalues)
dfS=cbind(t=probaS,BET=BETvalues)
get_groupBET(dfT,nb_total=100)
get_groupBET(dfS,nb_total=100)

# Other threshold calculations from literature
afcMatrix=getTriangularMatrix(triangular, subjectName="subject", productName="product",descriptorName="descriptor",timeName="time", increasingConcentrations=paste0("C",1:8))
p_correct=apply(afcMatrix,2,sum)/nrow(afcMatrix)
p_discr=(p_correct-1/3)/(2/3)
plot(p_discr,pch=16)

# Getting lawless results
res_lawless=lawlessCorrection(numericConcentrations=rev(exp(decreasingNumConcentrations)),
                              afcMatrix)
exp(res_lawless) 

# Getting hough correction
houghTriangular=houghCorrection(afcMatrix,p_chance=1/3)
