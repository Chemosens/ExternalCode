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

pS=rep(0,7);pS[4]=1
matTassumingS=AFCR:::getMatrixTsachantS(2/3,6)
pT_estimated=pS%*%matTassumingS
pT=as.vector(pT_estimated)
dfT=data.frame(t=pT,names=paste0("]c",0:6,";c",1:7,"]"))

dfS=data.frame(t=pS,names=paste0("]c",0:6,";c",1:7,"]"))
dfT[,"t_lab"]=round(dfT[,"t"],digits=2)
dfS[,"t_lab"]=round(dfS[,"t"],digits=2)
real=ggplot(dfS,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("a. Distribution of the true threshold")+theme_bw()+ylab("Probability")+xlab("Threshold")

obs=ggplot(dfT,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("b. Distribution of BET threshold")+theme_bw()+xlab("Threshold")+ylab("Probability")+ geom_text(aes(label = t_lab), vjust = -0.5)


pS_2=rep(0,7);pS_2[1]=pS_2[7]=0.04;pS_2[2]=pS_2[6]=0.09;pS_2[3]=pS_2[5]=0.22;pS_2[4]=0.3;
pT_estimated_2=pS_2%*%matTassumingS
pT_2=as.vector(pT_estimated_2)
dfT_2=data.frame(t=pT_2,names=paste0("]c",0:6,";c",1:7,"]"))
dfS_2=data.frame(t=pS_2,names=paste0("]c",0:6,";c",1:7,"]"))
dfT_2[,"t_lab"]=round(dfT_2[,"t"],digits=2)
dfS_2[,"t_lab"]=round(dfS_2[,"t"],digits=2)
real_2=ggplot(dfS_2,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("c. Distribution of the true threshold (2)")+theme_bw()+ylab("Probability")+xlab("Threshold") + ylim(0,0.4)+ geom_text(aes(label = t_lab), vjust = -0.5)
obs_2=ggplot(dfT_2,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("d. Distribution of BET threshold (2)")+theme_bw()+xlab("Threshold")+ylab("Probability")+ ylim(0,0.4)+ geom_text(aes(label = t_lab), vjust = -0.5)
grid.arrange(real,obs,real_2,obs_2,nrow=2)

# Other simulations for choosing the concentrations
#====================================================
# Here we imagine that the threshold concentration in the population follows a log-normal distribution.
m1=4.5
sd1=1
abscissa=seq(0,10000,length.out=100000)
probaNorm=dnorm(log(abscissa),mean=m1,sd=sd1)
probaLNorm=dlnorm(abscissa,meanlog=m1,sdlog=sd1)

# This example shows how the impact of the choice of the concentrations in the distribution of threshold
# We have chosen 6 successive AFC tests
K=6
cmax=exp(8)
cmin=exp(1)
# First, regular concentration are selected
log_concentrations_reg=seq(log(cmin),log(cmax),length.out=K+2)

# Then, concentrations are selected such as the threshold distribution is equiprobable
log_concentrations_equi=c(max(log(cmin),0),qnorm(mean=m1,sd=sd1,p=c((1:K)/(K+1))),log(cmax))
concentrations_reg=exp(log_concentrations_reg)
concentrations_equi=exp(log_concentrations_equi)

# Finally concentrations are selected with ASTM with 3 different minimal concentrations
cmin1=cmin;
cmin2=exp(3);
cmin3=exp(2)
concentrations_astm1=c(cmin1,2*cmin1,4*cmin1,8*cmin1,16*cmin1,32*cmin1,64*cmin1,128*cmin1)
concentrations_astm2=c(cmin2,2*cmin2,4*cmin2,8*cmin2,16*cmin2,32*cmin2,64*cmin2,128*cmin2)
concentrations_astm3=c(cmin3,2*cmin3,4*cmin3,8*cmin3,16*cmin3,32*cmin3,64*cmin3,128*cmin3)

log_concentrations_astm1=log(concentrations_astm1)
log_concentrations_astm2=log(concentrations_astm2)
log_concentrations_astm3=log(concentrations_astm3)

length(concentrations_equi)==length(concentrations_reg)
length(concentrations_equi)==length(concentrations_astm1)
length(concentrations_equi)==K+2

rbind(reg=concentrations_reg,equi=concentrations_equi,astm1=concentrations_astm1,astm2=concentrations_astm2,astm3=concentrations_astm3)
# Probability to be in these concentrations
quantile_equi=pnorm(log_concentrations_equi,mean=m1,sd=sd1)
quantile_reg=pnorm(log_concentrations_reg,mean=m1,sd=sd1)
quantile_astm1=pnorm(log_concentrations_astm1,mean=m1,sd=sd1)
quantile_astm2=pnorm(log_concentrations_astm2,mean=m1,sd=sd1)
quantile_astm3=pnorm(log_concentrations_astm3,mean=m1,sd=sd1)
quantile_reg[K+2]=quantile_equi[K+2]=quantile_astm1[K+2]=quantile_astm2[K+2]=quantile_astm3[K+2]=1

proba_equi=diff(quantile_equi)
proba_reg=diff(quantile_reg)
proba_astm1=diff(quantile_astm1)
proba_astm2=diff(quantile_astm2)
proba_astm3=diff(quantile_astm3)

t_equi=AFCR:::probaT(probaS=proba_equi,p=2/3)
t_reg=AFCR:::probaT(probaS=proba_reg,p=2/3)
t_astm1=AFCR:::probaT(probaS=proba_astm1,p=2/3)
t_astm2=AFCR:::probaT(probaS=proba_astm2,p=2/3)
t_astm3=AFCR:::probaT(probaS=proba_astm3,p=2/3)
# Total graph
colors=c("blue","lightgrey","cyan","deeppink","orange")


#The following graph represents the selected concentrations for a given normal distribution
plot(abscissa,probaLNorm,type="l",main="Threshold distribution according to concentrations",xlim=c(0,700),xlab="Concentrations")
abline(v=concentrations_equi,col=colors[1],lty=2)
abline(v=concentrations_reg,col=colors[2],lty=3,lwd=2)
abline(v=concentrations_astm1,col=colors[3],lty=1)
abline(v=concentrations_astm2,col=colors[4],lty=4)
abline(v=concentrations_astm3,col=colors[5],lty=4)

# Same with the log-normal distribution
plot(log(abscissa),probaNorm,type="l",main="Threshold distribution according to log-concentrations",xlim=c(0,9))
abline(v=log_concentrations_equi,col=colors[1],lty=2)
abline(v=log_concentrations_reg,col=colors[2],lty=3,lwd=2)
abline(v=log_concentrations_astm1,col=colors[3],lty=1)
abline(v=log_concentrations_astm2,col=colors[4],lty=4)
abline(v=log_concentrations_astm3,col=colors[5],lty=4)
points(log(cmin1),0,col=colors[3],pch=16)
points(log(cmin2),0,col=colors[4],pch=16)
points(log(cmin3),0,col=colors[5],pch=16)

barplot(proba_equi,names.arg=paste0("S_c",0:(K),"-c",1:(K+1)),col=colors[1], main="P(S=...) \n supposing a quantile repartition of log-conc.")
barplot(t_equi,names.arg=paste0("T_c",0:K,"-c",1:(K+1)),col=colors[1], main="P(T=...) \n supposing a quantile repartition of log-conc.")
barplot(proba_reg,names.arg=paste0("S_c",0:K,"-c",1:(K+1)), col=colors[2],main="P(S=...) \n supposing a regular repartition of log-conc.")
barplot(t_reg,names.arg=paste0("T_c",0:K,"-c",1:(K+1)), col=colors[2],main="P(T=...) \n supposing a regular repartition of log-conc.")
barplot(proba_astm1,names.arg=paste0("S_c",0:K,"-c",1:(K+1)),col=colors[3], main=paste0("P(S=...)  with an ASTM repartition \n beginning to ",round(cmin1,digits=1)))
barplot(t_astm1,names.arg=paste0("T_c",0:K,"-c",1:(K+1)),col=colors[3], main=paste0("P(T=...)  with an ASTM repartition \n beginning to ",round(cmin1,digits=1)))
barplot(proba_astm2,names.arg=paste0("S_c",0:K,"-c",1:(K+1)),col=colors[4], main=paste0("P(S=...) \n with an ASTM repartition \n beginning to ",round(cmin2,digits=1)))
barplot(t_astm2,names.arg=paste0("T_c",0:K,"-c",1:(K+1)),col=colors[4], main=paste0("P(T=...) \n with an ASTM repartition \n beginning to ",round(cmin2,digits=1)))
barplot(proba_astm3,names.arg=paste0("S_c",0:K,"-c",1:(K+1)),col=colors[5], main=paste0("P(S=...) \n with an ASTM repartition \n beginning to ",round(cmin3,digits=1)))
barplot(t_astm3,names.arg=paste0("T_c",0:K,"-c",1:(K+1)),col=colors[5], main=paste0("P(T=...) \n with an ASTM repartition \n beginning to ",round(cmin3,digits=1)))

rmse_equi=1/(K+2)*sum(100*(t_equi-proba_equi)^2)
rmse_reg=1/(K+2)*sum(100*(t_reg-proba_reg)^2)
rmse_astm1=1/(K+2)*sum(100*(t_astm1-proba_astm1)^2)
rmse_astm2=1/(K+2)*sum(100*(t_astm2-proba_astm2)^2)
rmse_astm3=1/(K+2)*sum(100*(t_astm3-proba_astm3)^2)

rmse_equi
rmse_reg
rmse_astm1
rmse_astm2
rmse_astm3

# The equiprobable repartition returns the smallest RMSE. Among ASTM, the 3rd one is the better choice. 



#======================================================
# Code for the real applications presented in the paper
#======================================================

setwd("C:/Users/capeltier/Nextcloud/MyDrive/DataAnalysis/2023_ThèseAlix/AFC")
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
dfT=data.frame(t=pT,names=paste0("]c",0:8,";c",1:9,"]"),t_lab=round(pT,digits=2))
obs=ggplot(dfT,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("a. Observed BET threshold distribution")+theme_bw()+xlab("Threshold")+ylab("Probability")+ geom_text(aes(label = t_lab), vjust = -0.5)

# Calculating corrected BET distribution
conditional_proba=AFCR:::getMatrixTsachantS(2/3,8)
round(conditional_proba,digits=2)
probaS=t(ginv(conditional_proba))%*%probaT
summary(sum(probaS)-1) # Is a probability
probaS=as.vector(probaS)
names(probaS)=paste0("]c_", 0:8,";c_",1:9,"]")

# Plotting corrected BET distribution
pS=probaS
dfS=data.frame(t=pS,names=paste0("]c",0:8,";c",1:9,"]"),t_lab=round(pS,digits=2))
real=ggplot(dfS,aes(x=names,y=t))+geom_col()+ylim(0,1)+theme_bw()+ylab("Probability")+xlab("Threshold")+ggtitle("b. Estimated corrected threshold distribution")+ geom_text(aes(label = t_lab), vjust = -0.5)
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
res_lawless=lawlessCorrection(numericConcentrations=rev(decreasingNumConcentrations0),
                              afcMatrix)
10^(res_lawless$p50) 

# Getting hough correction
houghTriangular=houghCorrection(afcMatrix,p_chance=1/3)
p_dis=apply(houghTriangular$triangularMatrix,2,sum)/nrow(houghTriangular$triangularMatrix)
plot(rev(log(decreasingNumConcentrations0,base=10)),p_dis, pch=16,xlab="Log-concentrations",ylab="Proportion of discriminators",main="Hough correction")
coeffLm=lm(p_dis~rev(log(decreasingNumConcentrations0,base=10)))$coef
abline(a=coeffLm[1],b=coeffLm[2])
abline(h=0.5)
logConc=(0.5-coeffLm[1])/coeffLm[2]
10^logConc


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
