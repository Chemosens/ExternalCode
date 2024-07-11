library(devtools)
library(openxlsx)
library(plotly)
library(gridExtra)
library(openxlsx)
library(AFCR)
library(MASS)
# To run only for updating the package
#remove.packages("AFCR")
#install_github("https://github.com/ChemoSens/AFCR")


# Simulations: application 1
#===========================
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

# Real applications
#=======================
#setwd("C:/Users/capeltier/Nextcloud/MyDrive/DataAnalysis/2023_ThèseAlix/AFC")
# Charger les fichiers de timesens
triangular=read.xlsx("SEANCE 22 09 2023.xlsx",sheet="TDS")
triangular=triangular[,c("Répétition","Panéliste","Produit","Descripteur","Score","Temps")]
colnames(triangular)=c("replicate","subject","product","descriptor","score","time")

# Getting initial BET distribution
res_klo=keepLastOccurence(triangular,subjectName="subject",productName="product",descriptorName="descriptor",timeName="time")
threshold_data=getThreshold(res=res_klo,rata=NULL,
                            decreasingConcentrations=c("C8","C7","C6","C5","C4","C3","C2","C1"))
distributionThreshold=summary(as.factor(threshold_data[,"threshold"]))
probaT=distributionThreshold/sum(distributionThreshold)
pT=as.vector(probaT)
dfT=data.frame(t=pT,names=paste0("]c",0:8,";c",1:9,"]"))
obs=ggplot(dfT,aes(x=names,y=t))+geom_col()+ylim(0,1)+ggtitle("a. Observed BET threshold distribution")+theme_bw()+xlab("Threshold")+ylab("Probability")

# Calculating corrected BET distribution
conditional_proba=AFCR:::getMatrixTsachantS(2/3,8)
probaS=t(ginv(conditional_proba))%*%probaT
summary(sum(probaS)-1) # Is a probability
probaS=as.vector(probaS)
names(probaS)=paste0("]c_", 0:8,";c_",1:9,"]")
pS=probaS
dfS=data.frame(t=pS,names=paste0("]c",0:8,";c",1:9,"]"))
real=ggplot(dfS,aes(x=names,y=t))+geom_col()+ylim(0,1)+theme_bw()+ylab("Probability")+xlab("Threshold")+ggtitle("b. True threshold distribution after correction")
grid.arrange(obs+ylim(0,0.35),real+ylim(0,0.35),nrow=1)

# writing the results
#df_ts=cbind(dfT[,c(2)],round(100*dfT[,1],digits=2),round(100*dfS[,1],digits=2))
#write.table(df_ts,"test.csv",sep=";")

# Group threshold calculation
#-------------------------------
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
groupBET(dfT,nb_total=193)
groupBET(dfS,nb_total=193)

# Other threshold calculations from literature
triangularMatrix=getTriangularMatrix(triangular, subjectName="subject",descriptorName="descriptor", productName="product", timeName="time",increasingConcentrations=paste0("C",1:8))
p_correct=apply(triangularMatrix,2,sum)/nrow(triangularMatrix)
p_discr=(p_correct-1/3)/(2/3)
plot(p_discr,pch=16)

# Getting lawless results
res_lawless=lawlessCorrection(numericConcentrations=rev(decreasingNumConcentrations0),
                              triangularMatrix)
exp(res_lawless) 

# Getting hough correction
sim=rep(NA,100)
for(k in 1: 100)
{
  print(k)
  houghTriangular=houghCorrection(triangularMatrix,p_chance=1/3,numericConcentrations=rev(decreasingNumConcentrations0))
  p_dis=apply(houghTriangular$triangularMatrix,2,sum)/nrow(houghTriangular$triangularMatrix)
  plot(rev(log(decreasingNumConcentrations0)),p_dis, pch=16,xlab="Log-concentrations",ylab="Proportion of discriminators",main="Hough correction")
 # spl=splinefun(p_dis,rev(log(decreasingNumConcentrations0)))
  coeffLm=lm(p_dis[5:8]~rev(log(decreasingNumConcentrations0))[5:8])$coef
  abline(a=coeffLm[1],b=coeffLm[2])
  abline(h=0.5)
  logConc=(0.5-coeffLm[1])/coeffLm[2]
  sim[k]=exp(logConc)
  
}
summary(sim)
plot(rev(log(decreasingNumConcentrations0)),p_dis, pch=16,xlab="Log-concentrations",ylab="Proportion of discriminators",main="Hough correction")
points(rev(log(decreasingNumConcentrations0)),p_discr, pch=16,col="blue",xlab="Log-concentrations",ylab="Proportion of discriminators",main="Hough correction")
points(rev(log(decreasingNumConcentrations0)),cumsum(probaS), pch=16,col="red",xlab="Log-concentrations",ylab="Proportion of discriminators",main="BET distribution correction")

# ff
res=selectExtreme(N=193,n=15,lower=TRUE,p=1/3,probaS,K=8)
res$gg+ggtitle("Distribution of true threshold in the extreme group")
