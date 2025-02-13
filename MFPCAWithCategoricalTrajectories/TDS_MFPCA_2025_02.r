#######################################################################################################
# Statistical modeling of categorical trajectories with multivariate functional principal components
# H. Cardot, C. Peltier
#
# TDS application
# Version of the 2025/02/13
#######################################################################################################

# This scripts allows the different results of the paper regarding TDS data to be reproduced


# Loading libraries and data
#====================
library(cfda)
library(ggplot2)
library(openxlsx)
library(gridExtra)
library(plotly)
library(reshape2)
library(xtable)
library(funData)
library(MFPCA)
# The following functions are in the git repository https://github.com/Chemosens/ExternalCode/MFPCAWithCategoricalTrajectories.
source("convertForCfda.r")
source("plotSignal.r")
source("getIndicatrices.r")
source("tdsRead.r")
source("cleanDf.r")
source("splitInPeriods.r")# Reading data

# datapaper.xlsx is obtained on https://data.mendeley.com/datasets/3j6h7mrxnf/1
# It should be cited as Visalli, Michel; Béno, Noëlle (2023), “A dataset of consumer perceptions of gustometer-controlled stimuli measured with three temporal sensory evaluation methods”, Paris-Saclay University, V1, doi: 10.17632/3j6h7mrxnf.1

signalSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="Stimuli")
tdsSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="TDS")


colors=c(Acid="deeppink1",Basil="darkgreen",Bitter="dodgerblue",Lemon="yellow",Licorice="brown",
         Mint="lightgreen",Salty="lightblue",stop="white",Sweet="red")
colors_nb <- c(Acid="gray5",Basil="gray20",Bitter="gray30",Lemon="gray40",Licorice="gray50",
         Mint="gray60",Salty="gray80",stop="gray80",Sweet="gray90")
line_type_nb <- c(Acid="twodash",Basil="dotted",Bitter="longdash",Lemon="solid",Licorice="longdash",
               Mint="twodash",Salty="longdash",stop="blank",Sweet="solid")


# Analyzing TDS
#==================
# Renaming columns and reshaping data
tdsAll1=tdsSheet[,1:6]
colnames(tdsAll1)=c("rep","subject","product","time","descriptor","score")
tdsA=tdsRead(df=tdsAll1,discretization=1) # had to choose discretization=1 for some signals

# Discriminating 3 signals
#==============================
signal="S06"
signal2="S07"
signal3="S04"
prod1=tdsA$df[tdsA$df[,"product"]==signal,]
prod2=tdsA$df[tdsA$df[,"product"]==signal2,]
prod3=tdsA$df[tdsA$df[,"product"]==signal3,]
prod1[,"subject"]=paste0("A",substr(prod1[,"subject"],5,6))
prod2[,"subject"]=paste0("B",substr(prod2[,"subject"],5,6))
prod3[,"subject"]=paste0("C",substr(prod3[,"subject"],5,6))
all_dat0S1=convertForCfda(tds=NULL,tdsdf=prod1)
all_dat0S2=convertForCfda(tds=NULL,tdsdf=prod2)
all_dat0S3=convertForCfda(tds=NULL,tdsdf=prod3)

# Simplifying the name of subjects
all_dat0S1[,"id"]=substr(all_dat0S1[,"id"],4,6)
all_dat0S2[,"id"]=substr(all_dat0S2[,"id"],4,6)
all_dat0S3[,"id"]=substr(all_dat0S3[,"id"],4,6)

# Raw data: Getting Graph TDS (Fig. 1)  
#=================================
p5ter=plotData(all_dat0S1,col=colors)
p6ter=plotData(all_dat0S2,col=colors)
p7ter=plotData(all_dat0S3,col=colors)
#pdf("TDSS060704_col.pdf",width=10,height=7) 
grid.arrange(p5ter+ggtitle(signal),p6ter+ggtitle(signal2),p7ter+ggtitle(signal3),nrow=1)
#dev.off()

# Visualizing the signals (Fig. 2)
#=========================
signals=names(summary(factor(signalSheet[,"Product"])))
signal=signals[1]
### graphe IWFOS 2025 couleurs
p10 <- plotSignal(signalSheet,signals[10],colors=colors,linetype=line_type_nb)+theme_bw()
p11 <- plotSignal(signalSheet,signals[11],colors=colors,linetype=line_type_nb)+theme_bw()
p8 <- plotSignal(signalSheet,signals[8],colors=colors,linetype=line_type_nb)+theme_bw()
#pdf("SignalS06S07S04_col.pdf",width=10,height=7)
grid.arrange(p10,p11,p8,nrow=3)
#dev.off()



# Preprocessing for MFPCA No replicates, normalisation
all_dat_pls <- convertForCfda(tds=NULL,tdsdf=rbind(prod1,prod2,prod3))
nom_ind <- levels(as.factor(all_dat_pls$id))
n_ind <- length(nom_ind) ### 150
ind_stop <- all_dat_pls$state == "stop"
df_norm1 <- all_dat_pls[!ind_stop,] # removing STOP
etats <- levels(as.factor(df_norm1$state)) ## levels of state modality
df_norm1$id <- as.character(df_norm1$id) ## names of individuals
df_norm1$state <- as.character(df_norm1$state) ## modality of the states
n_state <- length(unique(df_norm1$state)) 
pos_state <- 1:n_state
temps.out <- seq(0,1,length=100)
n.points <- length(temps.out)

# Getting the list required for MFPCA
df_mfpca <- list()
for (i in 1:n_state){df_mfpca[[i]] <- matrix(0,nrow=n_ind, ncol=n.points)}
names(df_mfpca) <- etats
for (i in 1:n_ind)
{
  df_temp <- df_norm1[df_norm1$id==nom_ind[i],]
  mat_temp <- getIndicatrices(df_temp$state,df_temp$time,times.out=temps.out)
  state_temp <- unique(df_temp$state) 
  for (j in 1:length(state_temp)){
    nom_etat <- state_temp[j]
    pos_liste <- pos_state[etats == nom_etat]
    df_mfpca[[pos_liste]][i,] = mat_temp[j,]
  }
}

# Getting the final object to be passed in MFPCA
df_fundata <- list()
for (i in 1:n_state){df_fundata[[i]] <- funData(argvals=temps.out,X=df_mfpca[[i]])}
names(df_fundata)=names(df_mfpca)
df_multiFun <- multiFunData(df_fundata,values=names(df_fundata))
summary(df_multiFun)


# Exploration of the final object with MFPCA package and example of binary trajectories
p1 <- plotData(all_dat0S3[1:14,],col=colors) +theme_bw()
p1
plot(df_multiFun[["Bitter"]][1],xlab="Time",col=colors["Bitter"],lty=line_type_nb["Bitter"],lwd=3)
grid()
plot(df_multiFun[["Lemon"]][1],xlab="Time",col=colors["Lemon"],lty=line_type_nb["Lemon"],lwd=3, add=TRUE)
plot(df_multiFun[["Acid"]][1],col=colors["Acid"],lty=line_type_nb["Acid"],lwd=3, add=TRUE)
plot(df_multiFun[["Sweet"]][1],col=colors["Sweet"],lty=line_type_nb["Sweet"],lwd=3, add=TRUE)

CAt <- c("Bitter", "Lemon", "Acid", "Sweet")
legend(0.6,0.5,CAt,col=colors[CAt],lty=line_type_nb[CAt],lwd=2)

# MFPCA with equal weights ( Use the MFPCA with a basis of k=10 splines)
k=10 ## number of basis functions
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
summary(resMFPCA) # outputs of MFPCA

# Graphs of Eigenvalues (Figure 3a.)
#=====================================
plot(resMFPCA$values)
round(resMFPCA$values/sum(resMFPCA$values),digits=2)



# Plot of mean functions (Figure 4)
#=============================
pr_descripteur <- unlist(lapply(df_mfpca,mean))
nom <- names(pr_descripteur)
round(pr_descripteur,digits=2)
#Acid    Basil   Bitter    Lemon Licorice     Mint    Salty    Sweet 
#0.20     0.06     0.07     0.16     0.01     0.01     0.17     0.32 

#pdf("mean_pr_col.pdf",height=7,width=10)
par(mfrow=c(1,1))
meanFuncMFPCA <- resMFPCA$meanFunction
plot(meanFuncMFPCA [["Acid"]], ylim=c(0,0.5),col=colors["Acid"],lty=line_type_nb["Acid"],lwd=3, xlab="Time",ylab="Empirical probability") ## Acid
plot(meanFuncMFPCA [["Sweet"]],add=TRUE,col=colors["Sweet"],lty=line_type_nb["Sweet"],lwd=3) ## Sweet
plot(meanFuncMFPCA [["Lemon"]],add=TRUE,col=colors["Lemon"],lty=line_type_nb["Lemon"],lwd=3) ## 
plot(meanFuncMFPCA [["Salty"]],add=TRUE,col=colors["Salty"],lty=line_type_nb["Salty"],lwd=3) ## Salty
plot(meanFuncMFPCA [["Basil"]],add=TRUE,col=colors["Basil"],lty=line_type_nb["Basil"],lwd=3) ## 
plot(meanFuncMFPCA [["Bitter"]],add=TRUE,col=colors["Bitter"],lty=line_type_nb["Bitter"],lwd=3) ## Bitter
grid()
title(main="Empirical mean trajectories")
interestingAttributes=c("Acid", "Sweet","Lemon","Salty","Basil","Bitter") # Selection of most frequent states
legend(0.8,0.25,interestingAttributes,col=colors[interestingAttributes],lty=line_type_nb[interestingAttributes],lwd=2)
#dev.off()

# Getting Principal components (Figure 5)
#==========================================
pc_mfpca <- resMFPCA$scores
colnames(pc_mfpca) <- paste0("pc",1:ncol(pc_mfpca))
pc <- as.data.frame(pc_mfpca)
pc[,"Signal"] <- substr(nom_ind,1,3)
pc[,"ind"] <- nom_ind
#pdf("cp_MFPCA_S04-S06-S07_col.pdf",height=7,width=10)
gg_ind <- ggplot(data=pc,aes(x=pc1,y=pc2,color=Signal))+geom_point(size=3)+theme_classic()+ggtitle("Principal components scores ", paste0(signal,"-",signal2,"-",signal3))
gg_ind + geom_hline(yintercept=0,col="gray80") + geom_vline(xintercept=0,col="gray80")+ labs(x="PC1 (24%)", y="PC2 (12%)") + scale_color_grey() + theme_classic()#+ scale_color_brewer()
#dev.off()

# Indicators of variable importance
#===================================
eigenfunMFPCA <- resMFPCA$functions
#### Indicator of variable importance
Poids_var <- matrix(NA,nrow=8,ncol=5)
for (i in 1:8){
  Poids_var[i,] <- apply(eigenfunMFPCA[[i]]@X[1:5,]^2,1,mean)
}
colSums(Poids_var)
colnames(Poids_var) <- c("dim 1", "dim 2", "dim 3", "dim 4", "dim 5")
rownames(Poids_var) <- names(eigenfunMFPCA)
xtable(Poids_var)
round(Poids_var,digits=2)
#dim 1 dim 2 dim 3 dim 4 dim 5
#Acid      0.08  0.26  0.41  0.27  0.32
#Basil     0.04  0.07  0.00  0.04  0.03
#Bitter    0.00  0.02  0.02  0.02  0.02
#Lemon     0.10  0.02  0.48  0.14  0.27
#Licorice  0.00  0.00  0.00  0.00  0.00
#Mint      0.00  0.00  0.00  0.00  0.00
#Salty     0.22  0.30  0.02  0.16  0.09
#Sweet     0.56  0.34  0.06  0.37  0.27


for (nom in interestingAttributes){cat(nom, round(norm(eigenfunMFPCA[[nom]])[1:3],digits=2), "\n")}


# AXE 1 : Acid, Lemon, Salty, Sweet
# AXE 2 : Acid, Basil, Salty, Sweet
## AXE3 très interessant mais pas de place pour IWFOS

apply(resMFPCA$scores[,1:3],2,quantile)
#[,1]        [,2]         [,3]
#0%   -0.624761897 -0.56902618 -0.491894362
#25%  -0.358808950 -0.14250653 -0.123703365
#50%  -0.008880389  0.02906378 -0.008489644
#75%   0.364805922  0.19429476  0.115700641
#100%  0.619853730  0.53978355  0.524401579
 
apply(abs(resMFPCA$scores[,1:3]),2,quantile)
#[,1]        [,2]        [,3]
#0%   0.006397905 0.005370227 0.003087872
#25%  0.195083721 0.084990357 0.068325886
#50%  0.364045922 0.174990792 0.122062586
#75%  0.482712265 0.320824239 0.253602502
#100% 0.624761897 0.569026183 0.524401579

# Variation of the first eigenfunction around the average (Figure 6)
#==========================================
i_vp <- 1 #### premier vecteur propre
#pdf("vp1_sweetsalty_col.pdf",height=7,width=10)
par(mfrow=c(2,2), mar = c(2, 2, 2, 1), oma = c(2, 2, 2, 2))
perception <- "Acid"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Lemon"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Salty"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Sweet"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.3*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.3 ", phi[1])),expression(paste(hat(p), " - 0.3 ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)
#dev.off()       

# Variation of second eigenfunction around the average (Figure 7)
#==========================================
i_vp <- 2 #### premier vecteur propre
#pdf("vp2_sweetsalty_col.pdf",height=7,width=10)
par(mfrow=c(2,2), mar = c(2, 2, 2, 1), oma = c(2, 2, 2, 2))

perception <- "Acid"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Basil"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Salty"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Sweet"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.78),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-0.2*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + 0.2 ", phi[2])),expression(paste(hat(p), " - 0.2 ", phi[2]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)
#dev.off()

#Comparison of the weight values for the different weighted MFPCA
#=============================================================
ndesc <- length(pr_descripteur)
pr_descripteur <- unlist(lapply(df_mfpca,mean))
W.mat <- matrix(NA,nrow=3,ncol=ndesc)
W.mat[1,] <- 1/ndesc 
#W.mat[2,] <- wt_descripteur/sum(wt_descripteur)

Proba_emp <- matrix(NA,ncol=8,nrow=300)
W.mat[1,] <- 1/ndesc
for (i in 1:8){Proba_emp[,i] <- resMFPCA$meanFunction[[i]]@X}
w_temp <- 1/colMeans(Proba_emp*(1-Proba_emp))
W.mat[2,] <- w_temp/sum(w_temp)

wtp_descripteur <- (pr_descripteur)^(-1)
wtp_descripteur <- 1/colMeans(Proba_emp) 
W.mat[3,] <- wtp_descripteur/sum(wtp_descripteur)
colnames(W.mat) <- names(pr_descripteur)

round(100*W.mat,digits=1)
library(xtable)
xtable(W.mat)

# Here, the same work was conducted for weighted MFPCA
#=================================================
#  Weighted MFPCA w = 1/(p(1-p))
#==================================================
pr_descripteur <- unlist(lapply(df_mfpca,mean))
wt_descripteur <- (pr_descripteur*(1-pr_descripteur))^(-1)
resMFPCA_wt <- MFPCA(df_multiFun, weights = wt_descripteur,M = 75, uniExpansions = list(
  list(type = "splines1D", k = k),
  list(type = "splines1D", k = k),
  list(type = "splines1D", k = k),
  list(type = "splines1D", k = k),
  list(type = "splines1D", k = k),
  list(type = "splines1D", k = k),
  list(type = "splines1D", k = k),
  list(type = "splines1D", k = k)
),fit = TRUE)

# Figure 3.B
#==========
#pdf("valp_MFPCA.pdf",height=7,width=10)
par(mfrow=c(1,2))
plot(resMFPCA$values/sum(resMFPCA$values),ylim=c(0,0.24),main="Eigenvalues (%), equal weights",xlab="Dimension",ylab=expression(lambda),pch=19)
grid()
plot(resMFPCA_wt$values/sum(resMFPCA_wt$values),ylim=c(0,0.24),main="Eigenvalues (%), weights 1/(p(1-p))",xlab="Dimension",ylab=expression(lambda),pch=19 ) ### voir code apres
grid()
#dev.off()

# outputs of MFPCA
summary(resMFPCA_wt)
plot(resMFPCA_wt$values)
plot(resMFPCA_wt$values/sum(resMFPCA_wt$values))
round(resMFPCA_wt$values[1:5]/sum(resMFPCA_wt$values),digits=2)
#[1] 0.14 0.09 0.07 0.06 0.05

# Plot of the principal component scores (Figure 8)
#===================================================
pc_mfpca <- resMFPCA_wt$scores
colnames(pc_mfpca) <- paste0("pc",1:ncol(pc_mfpca))
pc <- as.data.frame(pc_mfpca)
pc[,"Signal"] <- substr(nom_ind,1,3)
pc[,"ind"] <- nom_ind

# Principal components
#pdf("cp_MFPCAW_S04-S06-S07.pdf",height=7,width=10)
#gg_ind <- ggplot(data=pc,aes(x=pc1,y=pc4,color=product,shape=ind_speciaux))+geom_point(size=3)+theme_bw()+ggtitle("MFPCA ", paste0(signal,"-",signal2,"-",signal3))
gg_ind <- ggplot(data=pc,aes(x=pc1,y=pc2,color=Signal))+geom_point(size=3)+theme_classic()+ggtitle("Principal components scores ", paste0(signal,"-",signal2,"-",signal3))
gg_ind + geom_hline(yintercept=0,col="gray80") + geom_vline(xintercept=0,col="gray80")+ labs(x="PC1 (14%)", y="PC2 (09%)") + scale_color_grey() + theme_bw()#+ scale_color_brewer()
#dev.off()

####### Eigenfunctions
### Importance of the States
eigenfunMFPCA_wt <- resMFPCA_wt$functions
Poids_varW <- matrix(NA,nrow=8,ncol=5)
for (i in 1:8){
  Poids_varW[i,] <- wt_descripteur[i]*apply(eigenfunMFPCA_wt[[i]]@X[1:5,]^2,1,mean)
}
colSums(Poids_varW)
colnames(Poids_varW) <- c("dim 1", "dim 2", "dim 3", "dim 4", "dim 5")
rownames(Poids_varW) <- names(eigenfunMFPCA_wt)

xtable(Poids_varW)
round(Poids_varW,digits=2)

#dim 1 dim 2 dim 3 dim 4 dim 5
#Acid      0.08  0.05  0.04  0.20  0.14
#Basil     0.21  0.20  0.09  0.03  0.02
#Bitter    0.00  0.02  0.14  0.07  0.19
#Lemon     0.15  0.04  0.06  0.10  0.38
#Licorice  0.01  0.02  0.01  0.56  0.14
#Mint      0.00  0.22  0.56  0.00  0.03
#Salty     0.19  0.36  0.05  0.00  0.04
#Sweet     0.36  0.08  0.05  0.03  0.05


## Axe 1 : Sweet, Salty, Lemon, Basil
## Axe 2 : Salty, Mint, Basil, 

# Variation of first eigenfunction around the average, weighted case (Figure 9)
#==========================================
meanFuncMFPCA <- resMFPCA_wt$meanFunction
eigenfunMFPCA <- resMFPCA_wt$functions

i_vp <- 1 #### premier vecteur propre
#pdf("vp1_weight_col.pdf",height=7,width=10)
par(mfrow=c(2,2), mar = c(2, 2, 2, 1), oma = c(2, 2, 2, 2))
alpha <- 1
perception <- "Basil"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.9),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.8,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Lemon"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.9),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.71,0.8,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Salty"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.9),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.8,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Sweet"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.9),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-alpha*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.8,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

#dev.off()

# Variation of first eigenfunction around the average, weighted case (Figure 10)
#==========================================
i_vp <- 2 #### premier vecteur propre
#pdf("vp2_weight_col.pdf",height=7,width=10)
par(mfrow=c(2,2), mar = c(2, 2, 2, 1), oma = c(2, 2, 2, 2))

perception <- "Basil"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.8),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.71,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Mint"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.8),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Salty"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.8),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

perception <- "Sweet"
plot(meanFuncMFPCA [[perception]], ylim=c(0,0.8),col=colors[perception],lty=line_type_nb[perception],lwd=3, xlab="Time",ylab="Empirical probability")
plot(meanFuncMFPCA [[perception]]+1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=5, col="gray1" ,lwd=2 )
plot(meanFuncMFPCA [[perception]]-1*eigenfunMFPCA[[perception]][i_vp],add=T,lty=2,col="gray50", lwd=2 )
title(perception)
legend(0.01,0.78,legend= c(expression(hat(p)),expression(paste(hat(p), " + ", phi[1])),expression(paste(hat(p), " - ", phi[1]))),col=c(colors[perception],"gray1","gray50"), lty=c(1,5,2),lwd=2)

#dev.off() 

