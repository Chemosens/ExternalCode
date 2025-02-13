#######################################################################################################
# Statistical modeling of categorical trajectories with multivariate functional principal components. 
# H. Cardot, C. Peltier
#
# TCATA Application
# Version of the 2025/02/13
#######################################################################################################

# This scripts allows the different results of the paper regarding TCATA data to be reproduced



library(cfda)
library(ggplot2)
library(openxlsx)
library(gridExtra)
library(plotly)
library(reshape2)
library(funData)
library(MFPCA)
# The following functions are in the git repository https://github.com/Chemosens/ExternalCode/MFPCAWithCategoricalTrajectories.
source("tcataRead.r")
source("plotSignal.r")
source("getIndicators.r")
# datapaper.xlsx is obtained on https://data.mendeley.com/datasets/3j6h7mrxnf/1 and should be cited as
# Visalli, Michel; Béno, Noëlle (2023), “A dataset of consumer perceptions of gustometer-controlled stimuli measured with three temporal sensory evaluation methods”, Paris-Saclay University, V1, doi: 10.17632/3j6h7mrxnf.1


# Reading data  (choosing a  working directory where )
#==================
setwd(wd)
signalSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="Stimuli")
tdsSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="TDS")
tcataSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="TCATA")

colors=c(Acid="deeppink1",Basil="darkgreen",Bitter="dodgerblue",Lemon="yellow",Licorice="brown",
         Mint="lightgreen",Salty="lightblue",stop="white",Sweet="red")

colors_nb <- c(Acid="gray5",Basil="gray20",Bitter="gray30",Lemon="gray40",Licorice="gray50",
         Mint="gray60",Salty="gray80",stop="gray80",Sweet="gray90")

line_type_nb <- c(Acid="solid",Basil="dotted",Bitter="longdash",Lemon="solid",Licorice="longdash",
               Mint="twodash",Salty="longdash",stop="blank",Sweet="solid")


# Analyzing TCATA
#==================
# This returns normalized tcata data with score 0 for 0 and 1
tcataAll1=tcataSheet[,1:6]
colnames(tcataAll1)=c("rep","subject","product","time","descriptor","score")
tcataA=tcataRead(tcataAll1)


# Pre-processing for MFPCA
#==============================
signal="S06"
signal2="S07"
signal3="S04"
df_norm1=tcataA$df[tcataA$df[,"product"]%in%c(signal,signal2,signal3),]
df_norm1[,"state"]=paste0(df_norm1[,"descriptor"],"_",df_norm1[,"score"])
nom_ind <- levels(as.factor( df_norm1$id))
n_ind <- length(nom_ind) ### 150
etats <- levels(as.factor(df_norm1$state)) 
df_norm1$id <- as.character(df_norm1$id) 
df_norm1$state <- as.character(df_norm1$state) 
n_state <- length(unique(df_norm1$state)) 
pos_state <- 1:n_state
temps.out <- seq(0,1,length=100)
n.points <- length(temps.out)

df_mfpca <- list()
for (i in 1:n_state){df_mfpca[[i]] <- matrix(0,nrow=n_ind, ncol=n.points);
rownames(df_mfpca[[i]])=nom_ind;}
names(df_mfpca) <- etats

descriptors=unique(df_norm1[,"descriptor"])
for (i in 1:n_ind) # Here, unlike TDS, it must be filled according to the score and the descriptor
{
  for(descriptor in descriptors)
  {
    df_temp <- df_norm1[df_norm1$id==nom_ind[i]&df_norm1$descriptor==descriptor,]
    mat_temp <- getIndicators(df_temp$state,df_temp$time,times.out=temps.out)
   state_temp <- unique(df_temp$state) 
   for (j in 1:length(state_temp)){
    nom_etat <- state_temp[j]
    pos_liste <- pos_state[etats == nom_etat]
    df_mfpca[[pos_liste]][i,] = mat_temp[j,]
  }  
  }
}
all(df_mfpca[["Acid_0"]]+df_mfpca[["Acid_1"]]==1)
all(df_mfpca[["Sweet_0"]]+df_mfpca[["Sweet_1"]]==1)
df_fundata <- list()
descriptors_1=paste0(descriptors,"_1")
for (desc in descriptors_1){
  df_fundata[[desc]] <- funData(argvals=temps.out,X=df_mfpca[[desc]])
}
names(df_fundata)=descriptors_1
df_multiFun <- multiFunData(df_fundata,values=names(df_fundata))
summary(df_multiFun)
plot(df_multiFun[["Acid_1"]],main="Acid",xlab="Time")

# Empirical proabability (Figure 11)
#====================================================
pr_descripteur <- unlist(lapply(df_mfpca,mean))
pdf("TCATAmean_pr_NB.pdf",height=7,width=10)
meanFuncMFPCA <- resMFPCA$meanFunction
plot(meanFuncMFPCA [["Acid_1"]], ylim=c(0,0.5),col=colors["Acid"],lty=line_type_nb["Acid"],lwd=3, xlab="Time",ylab="Empirical probability") ## Acid
plot(meanFuncMFPCA [["Sweet_1"]],add=TRUE,col=colors["Sweet"],lty=line_type_nb["Sweet"],lwd=3) ## Sweet
plot(meanFuncMFPCA [["Lemon_1"]],add=TRUE,col=colors["Lemon"],lty=line_type_nb["Lemon"],lwd=3) ## 
plot(meanFuncMFPCA [["Salty_1"]],add=TRUE,col=colors["Salty"],lty=line_type_nb["Salty"],lwd=3) ## Salty
plot(meanFuncMFPCA [["Basil_1"]],add=TRUE,col=colors["Basil"],lty=line_type_nb["Basil"],lwd=3) ## 
plot(meanFuncMFPCA [["Bitter_1"]],add=TRUE,col=colors["Bitter"],lty=line_type_nb["Bitter"],lwd=3) ## Bitter
#plot(meanFuncMFPCA [["Mint_1"]],add=TRUE,col=colors["Mint"],lty=line_type_nb["Bitter"],lwd=3) ## Bitter
#plot(meanFuncMFPCA [["Licorice_1"]],add=TRUE,col=colors["Licorice"],lty=line_type_nb["Bitter"],lwd=3) ## Bitter
title(main="Empirical mean trajectories")
interestingAttributes=c("Acid", "Sweet","Lemon","Salty","Basil","Bitter")
legend(0,0.5,interestingAttributes,col =colors[interestingAttributes],lty=line_type_nb[interestingAttributes],lwd=2)
dev.off()


# Use the MFPCA with a basis of k=10 splines
#============================================
k=10
resMFPCA <- MFPCA(df_multiFun, M = 75, uniExpansions = list(
  list(type = "splines1Dpen", k = k),
  list(type = "splines1Dpen", k = k),
  list(type = "splines1Dpen", k = k),
  list(type = "splines1Dpen", k = k),
  list(type = "splines1Dpen", k = k),
  list(type = "splines1Dpen", k = k),
  list(type = "splines1Dpen", k = k),
  list(type = "splines1Dpen", k = k)
),fit = TRUE) # calculate reconstruction, too
# outputs of MFPCA
#summary(resMFPCA)

# Eigenvalues of MFPCA with TCATA (Figure 12)
#=====================================
plot(resMFPCA$values)
round(resMFPCA$values/sum(resMFPCA$values),digits=2)
pdf("TCATAEigenvalues.pdf",height=7,width=10)
plot(resMFPCA$values/sum(resMFPCA$values),ylim=c(0,0.2),pch=19,main="Eigenvalues (%), equal weights",xlab="Dimension",ylab=expression(lambda))
grid()
dev.off()
#scoreplot(resMFPCA,main="scoreplot")

# Select only some attributes in the meanFunction
#plot(resMFPCA$meanFunction,dim=which(names(resMFPCA$meanFunction)%in%c("Acid","Basil"))) ##il faut zoomer
pc_mfpca <- resMFPCA$scores


# Principal components (individual graph) (Figure 13)
#==============================================
colnames(pc_mfpca) <- paste0("pc",1:ncol(pc_mfpca))
pc <- as.data.frame(pc_mfpca)
pc[,"Signal"] <- substr(nom_ind,10,12)
pc[,"ind"] <- nom_ind
pdf("TCATAcp_MFPCA_S04-S06-S07.pdf",height=7,width=10)
gg_ind <- ggplot(data=pc,aes(x=pc1,y=pc2,color=Signal))+geom_point(size=3)+theme_bw()+ggtitle("Principal components scores ", paste0(signal,"-",signal2,"-",signal3))
gg_ind + geom_hline(yintercept=0,col="gray80") + geom_vline(xintercept=0,col="gray80")+ labs(x="PC1 (19%)", y="PC2 (12%)")+ scale_color_grey()
dev.off()



# Importance of variable in axes
#==============================================
eigenfunMFPCA <- resMFPCA$functions
names(eigenfunMFPCA)
contribAx1=data.frame()
for (nom in names(eigenfunMFPCA))
{
  for(i in 1:3)
  {
    id=paste0(nom,i)
    contribAx1[id,"state"]=nom
    contribAx1[id,"dim"]=paste0("dim.",i)
    contribAx1[id,'value']=sum(eigenfunMFPCA[[nom]][i,]@X^2)/length(eigenfunMFPCA[[nom]][i,]@X)
  }
}
colorsForPlot=colors
names(colorsForPlot)=paste0(names(colors),"_1")
ggplot(contribAx1,aes(x=dim,y=value,fill=state))+geom_col()+scale_fill_manual(values=colorsForPlot)+theme_bw()


apply(resMFPCA$scores[,1:3],2,quantile)
apply(abs(resMFPCA$scores[,1:3]),2,quantile)


# First eigenfunction (Figure 14)
#==============================================
i_vp <- 1 #### premier vecteur propre

pdf("TCATAvp1_sweetsalty.pdf",height=7,width=10)
par(mfrow=c(2,2), mar = c(2, 2, 2, 1), oma = c(2, 2, 2, 2))
perception <- "Acid_1"
perception2=result <- sub("_.*", "", perception)
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception2)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Lemon_1"
perception2<- sub("_.*", "", perception)
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Salty_1"
perception2= sub("_.*", "", perception)
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception2)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Sweet_1"
perception2=result <- sub("_.*", "", perception)
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception2)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)
dev.off()       

# second eigenfunction (Figure 14)
#==============================================
i_vp <- 2
pdf("TCATAvp2_sweetsalty.pdf",height=7,width=10)
par(mfrow=c(2,2), mar = c(2, 2, 2, 1), oma = c(2, 2, 2, 2))

perception <- "Acid_1"
perception2=result <- sub("_.*", "", perception)
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception2)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Lemon_1"
perception2=result <- sub("_.*", "", perception)
plot(meanFuncMFPCA[[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception2)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Salty_1"
perception2=result <- sub("_.*", "", perception)
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception2)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Sweet_1"
perception2=result <- sub("_.*", "", perception)
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception2],lty=line_type_nb[perception2],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception2)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception2],"gray1","gray50"), lty=c(1,5,2),lwd=2)
dev.off()

