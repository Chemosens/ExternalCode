# Code for CFDA-TDS paper
#----- Caroline Peltier
#----- 2023 04 14
#========================================

# This code is separated in three parts:
# -the first one for simulations (Material and Method part),
# - the second one for the Result part
# -and the third one for supplementary data material


#libraries to load
#==================
library(scales)
library(cfda)
library(gridExtra)
library(ggplot2)
library(openxlsx)
library(mixOmics)
library(gridExtra)
library(ggrepel)
library(CSUtils)
library(FactoMineR)
library(chemosensR)
source("C:/Users/capeltier/Desktop/DataAnalysis/2022_QualitativeFDA/FunctionsCFDA_forPaper2.R")
reposave="C:/Users/capeltier/Desktop/Communications/Articles/CFDATDS/Review220308/"
repodata="C:/Users/capeltier/Desktop/DataAnalysis/2022_QualitativeFDA"
#==========================================
# Part 1: Material and methods: simulations (Fig 1-5)
#=========================================
# Plotting bases (Figure 1)

# Generating toy datasets
#=======================
# Generation of toy dataset 1
all_dat_P2S1=data.frame(state=c("A","B","C","C"),time=c(0,0.3,0.71,1),id=rep("P2S1",4))
all_dat_P2S3=data.frame(state=c("A","B","C","C"),time=c(0,0.35,0.72,1),id=rep("P2S3",4))
all_dat_P2S2=data.frame(state=c("A","B","C","C"),time=c(0,0.25,0.73,1),id=rep("P2S2",4))
all_dat_P1S1=data.frame(state=c("B","A","C","C"),time=c(0,0.3,0.71,1),id=rep("P1S1",4))
all_dat_P1S2=data.frame(state=c("B","A","C","C"),time=c(0,0.35,0.72,1),id=rep("P1S2",4))
all_dat_P1S3=data.frame(state=c("B","A","C","C"),time=c(0,0.25,0.73,1),id=rep("P1S3",4))
all_dat1=rbind(all_dat_P2S1,all_dat_P2S2,all_dat_P2S3, all_dat_P1S1,all_dat_P1S2,all_dat_P1S3)
plotData(all_dat1)
# Generation of toy dataset 2
all_dat_P1S1=data.frame(state=c("F","B","F","F"),time=c(0,0.4,0.6,1),id=rep("P1S1",4))
all_dat_P1S2=data.frame(state=c("D","B","D","D"),time=c(0,0.45,0.65,1),id=rep("P1S2",4))
all_dat_P1S3=data.frame(state=c("E","B","E","E"),time=c(0,0.5,0.7,1),id=rep("P1S3",4))
all_dat_P2S1=data.frame(state=c("F","C","F","F"),time=c(0,0.401,0.601,1),id=rep("P2S1",4))
all_dat_P2S2=data.frame(state=c("D","C","D","D"),time=c(0,0.451,0.651,1),id=rep("P2S2",4))
all_dat_P2S3=data.frame(state=c("E","C","E","E"),time=c(0,0.501,0.701,1),id=rep("P2S3",4))
all_dat_P3S1=data.frame(state=c("F","A","F","F"),time=c(0,0.399,0.599,1),id=rep("P3S1",4))
all_dat_P3S2=data.frame(state=c("D","A","D","D"),time=c(0,0.449,0.649,1),id=rep("P3S2",4))
all_dat_P3S3=data.frame(state=c("E","A","E","E"),time=c(0,0.499,0.699,1),id=rep("P3S3",4))
all_b1=data.frame(state=c("E","B","D","B"),time=c(0,0.41,0.61,1),id=rep("P1S4",4))
all_b2=data.frame(state=c("E","C","D","B"),time=c(0,0.71,0.81,1),id=rep("P2S4",4))
all_b3=data.frame(state=c("E","A","D","B"),time=c(0,0.51,0.61,1),id=rep("P3S4",4))
all_dat2=rbind(all_dat_P1S1,all_dat_P1S2,all_dat_P1S3,#all_b1,
               all_dat_P2S1,all_dat_P2S2,all_dat_P2S3,#all_b2,
               all_dat_P3S1,all_dat_P3S2,all_dat_P3S3#,all_b3
)
#
colorsEx1=c(A="red",B="blue",C="green")
colorsEx2=rainbow(6);names(colorsEx2)=c("A","B","C","D","E","F")
#

# CFDA for toy datasets
#=======================
#Analysis of toy dataset 1 (Figure 3)
setwd(reposave)
enc_all1=getOptEnc(all_dat1,nbasis=8,norder=1,colors=colorsEx1)
p_harm1=plotHarm(enc_all1,harm=1,sizeLine="none",colors=colorsEx1,maxSize=4)
p_harm2=plotHarm(enc_all1,harm=2,sizeLine="none",colors=colorsEx1,maxSize=4)

#pdf(file=paste0(reposave,"/Figure2.pdf"),width=9,height=10)
grid.arrange(
  plotData(all_dat1,col=colorsEx1)+ggtitle("a. Toy dataset 1"),
             enc_all1$p_comp+xlab("Component 1")+ylab("Component 2")+ggtitle("b. Individual map")+xlim(c(-1,1))+ylim(c(-1,1)),
             p_harm1$p_exp+ggtitle("c. Harmonic 1"),
             p_harm2$p_exp+ggtitle("d. Harmonic 2"),
  enc_all1$p_eig+ggtitle("e. Cumulative eigenvalues")+theme(panel.grid = element_line(color = "grey92", linetype = 2), panel.background = element_blank())
)
#dev.off()

#Analysis of toy dataset 2 (Figure 3)
enc_all2=getOptEnc(all_dat2,nbasis=8,norder=1,colors=colorsEx2,repel=T)
p_harm1Ex2=plotHarm(enc_all2,harm=1,sizeLine="none",colors=colorsEx2,maxSize=4)
p_harm2Ex2=plotHarm(enc_all2,harm=2,sizeLine="none",colors=colorsEx2,maxSize=4)

#pdf(file=paste0(reposave,"/Figure3.pdf"),width=9,height=10)
grid.arrange(      plotData(all_dat2,col=rainbow(6))+ggtitle("b. Toy dataset 2"),

             enc_all2$p_comp+xlab("Component 1")+ylab("Component 2")+ggtitle("b. Individual map"),
             p_harm1Ex2$p_exp+ggtitle("c. Harmonic 1"),
             p_harm2Ex2$p_exp+ggtitle("d. Harmonic 2"),
             enc_all2$p_eig+ggtitle("e. Cumulative eigenvalues")+theme(panel.grid = element_line(color = "grey92", linetype = 2), panel.background = element_blank())
             )
#dev.off()

# SPLSDA on toy dataset 2
#=======================
nchar_prod=2
groups=as.factor(substr(rownames(enc_all2$fmca$pc),1,nchar_prod))
names(groups)=rownames(enc_all2$fmca$pc)
X=enc_all2$fmca$pc
colnames(X)=paste0("V",1:ncol(X))
nbComp=8
list.keepX=c(1:nbComp)
set.seed(seed=2)
tune.splsda.srbct <- tune.splsda(X=X[,1:nbComp], Y=groups, ncomp = 2, # calculate for first 4 components
                                 validation = "loo",
                                  test.keepX = list.keepX,
                                 cpus = 1) # allow for paralleliation to decrease runtime

ressplsda0=splsda(X=X[,1:nbComp],Y=groups,keepX=1)
df_ind=data.frame(ressplsda0$variates$X,label=rownames(ressplsda0$variates$X))
p_ind=ggplot(data=df_ind,aes(x=comp1,y=comp2,label=label))+geom_vline(xintercept=0,col="grey")+geom_hline(yintercept=0,col="grey")+theme_bw()+geom_text_repel()
p_h1=trajectoiresPropres(enc_all2,vecteurPropre=ressplsda0$loadings$X[,1],colors=rainbow(6),lwd=1)$p_exp
p_h2=trajectoiresPropres(enc_all2,vecteurPropre=ressplsda0$loadings$X[,2],colors=rainbow(6),lwd=1)$p_exp

# Figure 4
#pdf(file=paste0(reposave,"/Figure4.pdf"),width=9,height=7)
grid.arrange(plot(tune.splsda.srbct)+ggtitle("a. Tuning results"),
             p_ind+ggtitle("b. Individual map"),
             p_h1+ggtitle("c. Harmonic 1"),
             p_h2+ggtitle("d. Harmonic 2"),nrow=2)
#dev.off()
#=========================================
# Part 2:  Analysis on real dataset (Results, Fig 6-10)
#=========================================

#tdsAll=read.xlsx("C:/Users/capeltier/Desktop/DataAnalysis/2022_QualitativeFDA/datapaper.xlsx",sheet="TDS")
#tdsAll2=tdsAll[,-c(1,3)]
# tdsAll3=tdsAll2[tdsAll[,"ProductFamily"]=="Crisp",]
# colnames(tdsAll3)=c("SubjectCode","ProductCode","Time","AttributeCode","Score")
# write.xlsx(tdsAll3,file="tdsCrisps.xlsx")
#colnames(tdsAll2)=c("subject","product","time","descriptor","score")
#tdsAll2[,"rep"]=1
#tdsA=tdsRead(df=tdsAll2)
#save(tdsA,file="tdsA.RData")



# Loading the data
setwd(repodata)
load("tdsA.Rdata")
# Chips
products=c("C1","C2","C3","C4")

# restricting TDSa to only some products
tdsA$durations=tdsA$duration[tdsA$durations[,"product"]%in%products,]
tdsA$durations[,"product"]=as.factor(as.character(tdsA$durations[,"product"]))
tdsA$stdDominances=tdsA$stdDominances[tdsA$stdDominances[,"product"]%in%products,]
tdsA$stdDominances[,"product"]=as.factor(as.character(tdsA$stdDominances[,"product"]))
tdsA$stdDominances[,"descriptor"]=as.factor(as.character(tdsA$stdDominances[,"descriptor"]))
tdsA$df=tdsA$df[tdsA$df[,"product"]%in%products,]
tdsA$df[,"product"]=as.factor(as.character(tdsA$df[,"product"]))


# Colors with stop
colors=c("dodgerblue","dodgerblue4","forestgreen","red","darkorange1","darkolivegreen1","brown","darkorchid1","hotpink","white")
names(colors)=c("Bland","Crackly_Hard","Crispy","Fat","Melting","Potato","Roasted","Salty","Sticky_Pasty","stop")
# Colors without stop
colors2=c("dodgerblue","dodgerblue4","forestgreen","red","darkorange1","darkolivegreen1","brown","darkorchid1","hotpink")
names(colors2)=c("Bland","Crackly_Hard","Crispy","Fat","Melting","Potato","Roasted","Salty","Sticky_Pasty")

prod1=tdsA$df[tdsA$df[,"product"]==products[1],]
prod2=tdsA$df[tdsA$df[,"product"]==products[2],]
prod3=tdsA$df[tdsA$df[,"product"]==products[3],]
prod4=tdsA$df[tdsA$df[,"product"]==products[4],]
all_prod=rbind(prod1,prod2,prod3,prod4)#,prod5)

# Plotting the legend
plot(NULL)
legend(x="center",fill=colors2,names(colors2))
summary(as.factor(as.character(tdsA$df[tdsA$df[,"product"]%in%products,"descriptor"])))

# Outputs without CFDA (PCA and curves)
#===================================
# TDS curves (Figure 5)
dom=tdsA$stdDominances[tdsA$stdDominances[,"product"]%in% products,]
dom[,"product"]=as.factor(as.character(dom[,"product"]))
dom=dom[dom[,"descriptor"]%in%names(colors2),]
dom[,"descriptor"]=as.factor(as.character(dom[,"descriptor"]))
summary(dom[,"descriptor"])
p=tdsDominanceCurves(dom)$gg
p+scale_color_manual(values=colors[1:9])
p+geom_line(size=1)
# PCA with PCAgg (Figure 5)
pca=PCAgg(tdsA$durations[tdsA$durations[,"product"]%in%products,],representation="DistanceBiplot",expandBiplot=50,option="Covariance")
p3=plotPCAgg(pca,type="biplot",indSup=c("points"))
p3

# Result of Hotelling tests (FactoMineR)
library(reshape2)
dur0=tdsA$durations[tdsA$durations[,"product"]%in%products &tdsA$durations[,"descriptor"]%in%names(colors2),]
dur2=dcast(dur0,"product~descriptor",value.var="score", fun.aggregate=mean)
durSuj=dcast(dur0,"product+subject+rep~descriptor",value.var="score", fun.aggregate=mean)
durSuj2=durSuj[,-c(1:3)]
rownames(durSuj2)=paste0(durSuj[,1],durSuj[,2],durSuj[,3])
dur3=dur2[,-1];rownames(dur3)=dur2[,1]
dur4=rbind(dur3,durSuj2)
pcaFacto=FactoMineR::PCA(dur4,scale.unit=F,ind.sup=5:nrow(dur4))
colPca=rainbow(4); names(colPca)=c("C1","C2","C3","C4")
colPca[substr(rownames(durSuj2),1,2)]
colIndSup=colPca[substr(rownames(durSuj2),1,2)];names(colIndSup)=rownames(durSuj2)
summary(manova(pcaFacto$ind.sup$coord[,1:2]~as.factor(substr(rownames(pcaFacto$ind.sup$coord),1,2)),test="Wilks"))
h12=HotellingTest(pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C1",1:2],
              pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C2",1:2])
h23=HotellingTest(pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C3",1:2],
                  pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C2",1:2])
h13=HotellingTest(pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C1",1:2],
                  pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C3",1:2])
h14=HotellingTest(pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C1",1:2],
                  pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C4",1:2])
h24=HotellingTest(pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C4",1:2],
                  pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C2",1:2])
h34=HotellingTest(pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C4",1:2],
                  pcaFacto$ind.sup$coord[substr(rownames(pcaFacto$ind.sup$coord),1,2)=="C3",1:2])

print(paste0(as.character(round(h12$p.value,digits=4)),paste0(" (H = ",round(h12$stat,digits=2),")")))
print(paste0(as.character(round(h13$p.value,digits=4)),paste0(" (H = ",round(h13$stat,digits=2),")")))
print(paste0(as.character(round(h14$p.value,digits=4)),paste0(" (H = ",round(h14$stat,digits=2),")")))
print(paste0(as.character(round(h23$p.value,digits=4)),paste0(" (H = ",round(h23$stat,digits=2),")")))
print(paste0(as.character(round(h24$p.value,digits=4)),paste0(" (H = ",round(h24$stat,digits=2),")")))
print(paste0(as.character(round(h34$p.value,digits=4)),paste0(" (H = ",round(h34$stat,digits=2),")")))


# Converting TDS for CFDA
all_dat=convertForCfda(tds=NULL,tdsdf=all_prod)
# Visualzing the data
pt_evol1 <- estimate_pt(all_dat[substr(all_dat[,"id"],1,nchar(products)[1])==products[1],])
pt_evol2 <- estimate_pt(all_dat[substr(all_dat[,"id"],1,nchar(products)[2])==products[2],])
pt_evol3 <- estimate_pt(all_dat[substr(all_dat[,"id"],1,nchar(products)[3])==products[3],])
pt_evol4 <- estimate_pt(all_dat[substr(all_dat[,"id"],1,nchar(products)[4])==products[4],])

pt_evol_all <- estimate_pt(all_dat)
p1=plot(pt_evol1,col=colors,ribbon = TRUE)+ggtitle(paste0("P(X(t)=x) for ",products[1]))+scale_color_discrete(rainbow(n=5))
p2=plot(pt_evol2,col=colors,ribbon = TRUE)+ggtitle(paste0("P(X(t)=x) for ",products[2]))
p3=plot(pt_evol3,col=colors,ribbon = TRUE)+ggtitle(paste0("P(X(t)=x) for ",products[3]))
p4=plot(pt_evol4,col=colors,ribbon = TRUE)+ggtitle(paste0("P(X(t)=x) for ",products[4]))
p6=plot(pt_evol_all,col=colors[names(colors)!="stop"],ribbon=TRUE)+ggtitle("e. P(X(t)=x for all products")
grid.arrange(p1,p2,p3,p4,p6)

# Running CFDA on real TDS dataset (Fig. 6)
#===================================
all_dat_wr=all_dat
enc_all=getOptEnc(all_dat_wr,nbasis=8,norder=3)
pc=enc_all$fmca$pc
colnames(pc)=paste0("pc",1:ncol(pc))
pc=as.data.frame(pc)
pc[,"product"]=substr(rownames(pc),1,nchar_prod)
gg_ind=ggplot(data=pc,aes(x=pc1,y=pc2,color=product))+theme(panel.grid = element_line(color = "grey92",
                                                                                       size = 0.01,
                                                                                       linetype = 2),
                                                             panel.background = element_blank())+geom_vline(xintercept=0,color="grey")+geom_hline(yintercept=0,color="grey")+geom_vline(xintercept=0.5,color="darkgoldenrod1",linetype=2)+geom_hline(yintercept=0.5,color="darkgoldenrod1",linetype=2)+geom_point(aes(pch=product),size=2)
p_harm1=plotHarm(enc_all,harm=1,sizeLine="none",colors=colors2)
p_harm2=plotHarm(enc_all,harm=2,sizeLine="none",colors=colors2)
pt_evol_all <- estimate_pt(all_dat[all_dat[,"state"]!="stop",])
p5=plot(pt_evol_all,col=colors[names(colors)!="stop"],ribbon=TRUE)+ggtitle("e. P(X(t)=x for all products")
p6=plotData(all_dat[all_dat[,"id"]%in%rownames(pc[pc[,1]>0.5&pc[,2]>0.5,]),],col=colors2)+ggtitle("f. Raw data of the individuals \n in the upper right quadrant of 6.d")

mi=min(enc_all$fmca$pc[,1:2])
ma=max(enc_all$fmca$pc[,1:2])
#pdf(file=paste0(reposave,"/Figure6.pdf"),width=7,height=10)
grid.arrange(enc_all$p_eig+ggtitle("a. Cumulative eigenvalues"),
             gg_ind+ggtitle("b. Individual map")+xlab("Component 1")+ylab("Component 2")+xlim(mi,ma)+ylim(mi,ma),
             p_harm1$p_exp+theme(legend.position="none")+ggtitle("c. Harmonic 1"),
             p_harm2$p_exp+theme(legend.position="none")+ggtitle("d. Harmonic 2"),
             p5,p6,nrow=3)
#dev.off()


# Results of sPLSDA (Fig. 7)
#===================================
groups=as.factor(substr(rownames(enc_all$fmca$pc),1,nchar_prod))
names(groups)=rownames(enc_all$fmca$pc)
X=enc_all$fmca$pc
colnames(X)=paste0("V",1:ncol(X))

# Tuning sparsity parameter (/!\ can be long to run)
list.keepX=c(1:30)
set.seed(5)
tune.splsda.srbct <- tune.splsda(X, groups, ncomp = 2, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10,
                                 nrepeat = 1, # use repeated cross-validation
                                 test.keepX = list.keepX,
                                 cpus = 6) # allow for parallelisation to decrease runtime
# Results for Crisps # 7; 7; 3
p_tune=plot(tune.splsda.srbct)
plot(p_tune)
tune.splsda.srbct$choice.keepX

# Running sPLSDA
ressplsda=splsda(X=X,Y=groups,keepX=c(5,5),ncomp=2)

# Getting Hotelling results on the 1-2 map obtained
X2=ressplsda$variates[[1]]
hcfda12=HotellingTest(X2[substr(rownames(X2),1,2)=="C1",1:2],
                      X2[substr(rownames(X2),1,2)=="C2",1:2])
hcfda23=HotellingTest(X2[substr(rownames(X2),1,2)=="C3",1:2],
                      X2[substr(rownames(X2),1,2)=="C2",1:2])
hcfda13=HotellingTest(X2[substr(rownames(X2),1,2)=="C1",1:2],
                      X2[substr(rownames(X2),1,2)=="C3",1:2])
hcfda14=HotellingTest(X2[substr(rownames(X2),1,2)=="C1",1:2],
                      X2[substr(rownames(X2),1,2)=="C4",1:2])
hcfda24=HotellingTest(X2[substr(rownames(X2),1,2)=="C4",1:2],
                      X2[substr(rownames(X2),1,2)=="C2",1:2])
hcfda34=HotellingTest(X2[substr(rownames(X2),1,2)=="C3",1:2],
                      X2[substr(rownames(X2),1,2)=="C4",1:2])

print(paste0(as.character(round(hcfda12$p.value,digits=4)),paste0(" (H = ",round(hcfda12$stat,digits=2),")")))
print(paste0(as.character(round(hcfda13$p.value,digits=4)),paste0(" (H = ",round(hcfda13$stat,digits=2),")")))
print(paste0(as.character(round(hcfda14$p.value,digits=4)),paste0(" (H = ",round(hcfda14$stat,digits=2),")")))
print(paste0(as.character(round(hcfda23$p.value,digits=4)),paste0(" (H = ",round(hcfda23$stat,digits=2),")")))
print(paste0(as.character(round(hcfda24$p.value,digits=4)),paste0(" (H = ",round(hcfda24$stat,digits=2),")")))
print(paste0(as.character(round(hcfda34$p.value,digits=4)),paste0(" (H = ",round(hcfda34$stat,digits=2),")")))

# Results of sPLS (Fig. 7)
df_ind=as.data.frame(ressplsda$variates$X)
df_ind[,"product"]=groups
mi=min(df_ind[,1:2])
ma=max(df_ind[,1:2])
p_ind=ggplot(df_ind,aes(x=comp1,y=comp2,col=product))+geom_point()+theme_bw()+xlim(mi,ma)+ylim(mi,ma)+ggtitle("a. Individual plot")+geom_vline(xintercept=0,col="grey")+geom_hline(yintercept=0,col="grey")
p_h1=trajectoiresPropres(enc_all,vecteurPropre=ressplsda$loadings$X[,1],colors=colors,sizeLine="none")$p_exp+ggtitle("b. Harm. 1")+geom_hline(yintercept=0,col="grey")
p_h2=trajectoiresPropres(enc_all,vecteurPropre=ressplsda$loadings$X[,2],colors=colors,sizeLine="none")$p_exp+ggtitle("c. Harm. 2")+geom_hline(yintercept=0,col="grey")

setwd(reposave)
#pdf(file=paste0(reposave,"/Figure7.pdf"),width=7, height=6 )
grid.arrange(p_tune+ggtitle("a. Tuning results"),
             p_ind+ggtitle("b. Individual map")+theme(panel.grid = element_line(color = "grey92",
                                                                                size = 0.01,
                                                                                linetype = 2),
                                                      panel.background = element_blank()),
             p_h1+theme(legend.position="none")+ggtitle("c. Harmonic 1"),
             p_h2+theme(legend.position="none")+ggtitle("d. Harmonic 2"),
             nrow=2)
#graphics.off()
summary(manova(ressplsda$variates$X[,1:2]~groups),test="Hotelling-Lawley")
summary(manova(as.matrix(pca$indSup[,c("x","y")])~as.factor(pca$indSup[,"product"])),test="Hotelling-Lawley")

# Classification (Fig. 8)
#=====================
# Selecting only prod3=C3
prodC3=rbind(prod3)
prod_converted=convertForCfda(tds=NULL,tdsdf=prodC3)
plotData(prod_converted,col=colors)
# Running cfda
enc_allC3=getOptEnc(prod_converted,norder=3,nbasis=8)
pcC3=enc_allC3$fmca$pc
colnames(pcC3)=paste0("pc",1:ncol(pcC3))
pcC3=as.data.frame(pcC3)
pcC3[,"ind"]=substr(rownames(pcC3),3,8)
gg_ind=ggplot(data=pcC3,aes(x=pc1,y=pc2,label=ind))+geom_text()+theme_bw()
reshclust=hclust(dist(pcC3[,1:2]),method="ward.D2")
plot(100*(1-rev(reshclust$height)/sum(reshclust$height)),xlab="number of groups",ylab="Percent explained",main="Proportion of inertia explained by the agglomeration",pch=16);
groups=cutree(reshclust,k=3)
summary(factor(groups))
gp1= names(groups[groups==1])
gp2= names(groups[groups==2])
gp3= names(groups[groups==3])
sum(duplicated(c(substr(gp1[1:40],7,10),substr(gp1[41:76],3,7))))#20
sum(duplicated(c(substr(gp2[1:21],7,10),substr(gp2[21:44],3,7))))#5
sum(duplicated(c(substr(gp3[1:7],7,10),substr(gp3[8:17],3,7))))#2

pcC3[,"group"]=factor(groups)
summary(factor(groups))[1]
p01=gg_ind=ggplot(data=pcC3,aes(x=pc1,y=pc2,color=group))+geom_point()+theme(panel.grid = element_line(color = "grey92",
                                                                                                        size = 0.01,
                                                                                                        linetype = 2),
                                                                              panel.background = element_blank())+geom_vline(xintercept=0,col="darkgrey")+geom_hline(yintercept = 0,col="darkgrey")
p_harm1=plotHarm(enc_all=enc_allC3,harm=1,sizeLine="none",colors=colors,maxSize=4)
p_harm2=plotHarm(enc_all=enc_allC3,harm=2,sizeLine="none",colors=colors,maxSize=4)


p1=plotData(prod_converted[prod_converted[,"id"]%in%gp1,],col=colors)+ggtitle(paste0("d. Group 1 (n=",summary(factor(groups))[1],")"))+theme(legend.position="none")
p2=plotData(prod_converted[prod_converted[,"id"]%in%gp2,],col=colors)+ggtitle(paste0("e. Group 2 (n=",summary(factor(groups))[2],")"))+theme(legend.position="none")
p3=plotData(prod_converted[prod_converted[,"id"]%in%gp3,],col=colors)+ggtitle(paste0("f. Group 3 (n=",summary(factor(groups))[3],")"))+theme(legend.position="none")
setwd(reposave)
pdf(file=paste0(reposave,"/Figure8.pdf"),width=7, height=5 )
grid.arrange(p01+ggtitle("a. Individual map"),
             p_harm1$p_exp+ggtitle("b. Harmonic 1")+theme(legend.position="none"),
             p_harm2$p_exp+ggtitle("c. Harmonic 2")+theme(legend.position="none"),
             p1,p2,p3,
             nrow=2)
dev.off()
#======================================
# Supplementary data (Part 3)
#======================================

# Reconstructing a barplot of a subject
#======================================
Ncomp=70
res=list()
pctNcomp=rep(NA,Ncomp)
for(ncomp in 1:Ncomp)
{
  print(ncomp)
  res[[ncomp]]=reconstructBarplot(ncomp,enc_allC3,prod_converted)
  pctNcomp[ncomp]=sum(res[[ncomp]]$real_att_j==res[[ncomp]]$att_j,na.rm=T)/(dim(res[[ncomp]]$att_j)[1]*dim(res[[ncomp]]$att_j)[2])
}

res1=res[[1]]$df_res
res2=res[[2]]$df_res
res3=res[[3]]$df_res
res5=res[[5]]$df_res
res10=res[[10]]$df_res
res20=res[[20]]$df_res
res50=res[[50]]$df_res

p0=plotData(prod_converted,col=colors,addBorder=F)+ggtitle("h. 0riginal data")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p1=plotData(res1,col=colors,addBorder=F)+ggtitle("a. Estimation (p=1)")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p2=plotData(res2,col=colors,addBorder=F)+ggtitle("b. Estimation (p=2)")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p3=plotData(res3,col=colors,addBorder=F)+ggtitle("c. Estimation (p=3)")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p5=plotData(res5,col=colors,addBorder=F)+ggtitle("d. Estimation (p=5)")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p10=plotData(res10,col=colors,addBorder=F)+ggtitle("e. Estimation (p=10)")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p20=plotData(res20,col=colors,addBorder=F)+ggtitle("f. Estimation (p=20)")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p50=plotData(res50,col=colors,addBorder=F)+ggtitle("g. Estimation (p=50)")+theme(legend.position="none",axis.text.y = element_blank(),axis.ticks = element_blank())
p_prop=ggplot(data.frame(x=1:70,y=pctNcomp),aes(x=x,y=y))+geom_point()+theme_bw()+ggtitle("i. Quality of estimation")+ylab("Proportion of good estimates")+xlab("Number of components")
setwd(reposave)
tiff(filename="Figure11.tiff",width=1000, height=1000,res=100 )
grid.arrange(p1,p2,p3,p5,p10,p20,p50,p0,p_prop)
dev.off()


tds123=getTdsObjectFromReconstruction(res,c(1,5,50))
tds5to20=getTdsObjectFromReconstruction(res,c(2,10,60))
tds50=getTdsObjectFromReconstruction(res,c(3,20,70))
p=tdsDominanceCurves(dom[dom[,"product"]=="C3",])
p$gg+gg+scale_color_manual(values=colors[1:9])
tdsDominanceCurves(tds123$tds)$gg+ylim(c(0,70))+scale_color_manual(values=colors[1:9])
tdsDominanceCurves(tds5to20$tds)$gg+ylim(c(0,70))+scale_color_manual(values=colors[1:9])
tdsDominanceCurves(tds50$tds)$gg+ylim(c(0,70))+scale_color_manual(values=colors[1:9])

p=tdsDominanceCurves(dom[dom[,"product"]=="C3",])$gg
p+scale_color_manual(values=colors[1:9])+ylim(c(0,60))
