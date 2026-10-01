
library(ggplot2)
library(gridExtra)
library(openxlsx)
library(dplyr)
library(tidyr)
library(funData)
library(factoextra)
library(MFPCA)
# Reading data
#==================
# File downloaded from:
# Béno, Nicolle & Visalli (2023)
# https://doi.org/10.17632/3j6h7mrxnf.1... name it as 'datapaper.xlsx'

# Load the functions in 'Function repository'

lapply(
  list.files("Functions/", full.names = TRUE),
  source
)


signalSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="Stimuli")
tdsSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="TDS")
tcataSheet=openxlsx::read.xlsx("datapaper.xlsx",sheet="TCATA")

colors=c(Acid="deeppink1",Basil="darkgreen",Bitter="dodgerblue",Lemon="yellow",Licorice="brown",
         Mint="lightgreen",Salty="lightblue",stop="white",Sweet="red")


line_type_nb <- c(Acid="solid",Basil="dotted",Bitter="longdash",Lemon="solid",Licorice="longdash",
               Mint="twodash",Salty="longdash",stop="blank",Sweet="solid")


# Analyzing TCATA
#==================
# This returns normalized tcata data with score 0 for 0 and 1
tcataAll1=tcataSheet[,1:6]
colnames(tcataAll1)=c("rep","subject","product","time","descriptor","score")
# Adding all zeros for all descriptors at time =0 and time 1 for all descriptors
# Adds identifiers

# Selecting 3 signals with CFDA
#==============================
 signal="S06"
 signal2="S07"
 signal3="S04"
 df_spe=tcataAll1[tcataAll1[,"product"]%in%c(signal,signal2,signal3),]
df_spe[,"id"]=paste0(df_spe[,"subject"],"_",df_spe[,"product"])
dfnorm=tcataComplete(df_spe,norm=T)
#===========================
# Usual stats
 #=========================
# biplot
 durations <- dfnorm %>%
   arrange(id, descriptor, time) %>%
   group_by(id, descriptor) %>%
   mutate(
     next_time = lead(time),
     duration = ifelse(score == 1,
                       next_time - time,
                       0)
   ) %>%
   summarise(
     total_duration = sum(duration, na.rm = TRUE),
     .groups = "drop"
   )


 X <- durations %>%
   pivot_wider(
     names_from = descriptor,
     values_from = total_duration,
     values_fill = 0
   )
 X=as.data.frame(X)
 rownames(X) <- X$id
 product <- sub(".*_(S[0-9]+)_.*", "\\1", X$id)


res.pca <- prcomp(X[,-1], scale. = FALSE)
p_biplot=fviz_pca_biplot(
   res.pca,
   repel = TRUE,
   geom.ind = "point",
   habillage =product
 )
p_biplot

# curves

library(dplyr)
library(ggplot2)
df4=dfnorm[dfnorm[,"product"]=="S04",]
df6=dfnorm[dfnorm[,"product"]=="S06",]
df7=dfnorm[dfnorm[,"product"]=="S07",]
p_S04=tcataCurve(df4,times=seq(0, 1, by = 0.05),colors=colors)
p_S06=tcataCurve(df6,times=seq(0, 1, by = 0.05),colors=colors)
p_S07=tcataCurve(df7,times=seq(0, 1, by = 0.05),colors=colors)

grid.arrange(p_biplot+ggtitle("a. Biplot of durations"),
             p_S04+ggtitle("b. TCATA curve: P04"),
             p_S06+ggtitle("c. TCATA curve: P06"),
             p_S07+ggtitle("d. TCATA curve: P07"),nrow=2)
######### MFPCA ##########################################################

summary(dfnorm)
n.points=300
times.out=seq(0,max(dfnorm[,"time"]),length.out=  n.points)
df_multiFun=prepareForMFPCA(dfnorm,times.out )
plot(df_multiFun[["Acid"]],main="Acid",xlab="Time")

# Use the MFPCA with a basis of k=8 splines
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
),fit = TRUE)



# calculate reconstruction, too
# outputs of MFPCA
#summary(resMFPCA)
p_scree=screeplotMFPCA(resMFPCA,returnValue=FALSE)

eigen_values=screeplotMFPCA(resMFPCA,returnValue = TRUE)
cumsum(eigen_values[,"lambda"])
p_avg=plotAvgMFPCA(resMFPCA)
p_contrib=plotContribMFPCA(resMFPCA)

# plot PC
pc_mfpca <- resMFPCA$scores
colnames(pc_mfpca) <- paste0("pc",1:ncol(pc_mfpca))
pc <- as.data.frame(pc_mfpca)
nom_ind=rownames(pc_mfpca)
pc[,"Signal"] <- substr(nom_ind,10,12)
pc[,"ind"] <- nom_ind
gg_ind <- ggplot(data=pc,aes(x=pc1,y=pc2,color=Signal))+geom_point(size=3)+theme_bw()+ggtitle("Principal components scores ", paste0(signal,"-",signal2,"-",signal3))
expl_pc1=round(resMFPCA$values[1]/sum(resMFPCA$values)*100)
expl_pc2=round(resMFPCA$values[2]/sum(resMFPCA$values)*100)
p_ind=gg_ind + geom_hline(yintercept=0,col="gray80") + geom_vline(xintercept=0,col="gray80")+ labs(x=paste0("PC1 (",expl_pc1,"%)"), y=paste0("PC2 (",expl_pc2,"%)"))#+ scale_color_grey()

grid.arrange(p_avg+ggtitle("a. Empirical mean trajectories"),
             p_scree+ggtitle('b. Eigenvalue (proportion)'),
             p_contrib+ggtitle("c. Descriptor contributions"),
             p_ind+ggtitle("d. MFPCA scores"),nrow=2)

p_salty_1=plotEigenFunction(resMFPCA,perception="Salty",comp=1,colors=colors)
p_lemon_1=plotEigenFunction(resMFPCA,perception="Lemon",comp=1,colors=colors)
p_sweet_1=plotEigenFunction(resMFPCA,perception="Sweet",comp=1,colors=colors)
p_basil_1=plotEigenFunction(resMFPCA,perception="Basil",comp=1,colors=colors)

p_veigen_all=plotAllEigenFunction(resMFPCA,comp=1,colors=colors)

grid.arrange(p_lemon_1+ggtitle("a. Lemon"),p_salty_1+ggtitle("b. Salty"),p_sweet_1+ggtitle("c. Sweet"),p_basil_1+ggtitle("d. Basil"),nrow=2)
#### AXE 1
i_vp <- 2 #### second vecteur propre
#pdf("TCATAvp2_sweetsaltybasil.pdf",height=10,width=10)
p_acid_2=plotEigenFunction(resMFPCA,perception="Acid",comp=2,colors=colors)
p_lemon_2=plotEigenFunction(resMFPCA,perception="Lemon",comp=2,colors=colors)
p_salty_2=plotEigenFunction(resMFPCA,perception="Salty",comp=2,colors=colors)
p_sweet_2=plotEigenFunction(resMFPCA,perception="Sweet",comp=2,colors=colors)
p_basil_2=plotEigenFunction(resMFPCA,perception="Basil",comp=2,colors=colors)
grid.arrange(p_lemon_2+ggtitle("a. Lemon"),p_salty_2+ggtitle("b. Salty"),p_sweet_2+ggtitle("c. Sweet"),p_basil_2+ggtitle("d. Basil"),p_acid_2+ggtitle("e. Acid"),nrow=3)

i_vp <- 3 #### premier vecteur propre
#pdf("TCATAvp2_sweetsaltybasil.pdf",height=10,width=10)
p_acid_3=plotEigenFunction(resMFPCA,perception="Acid",comp=2,colors=colors)
p_lemon_3=plotEigenFunction(resMFPCA,perception="Lemon",comp=2,colors=colors)
p_salty_3=plotEigenFunction(resMFPCA,perception="Salty",comp=2,colors=colors)
p_sweet_3=plotEigenFunction(resMFPCA,perception="Sweet",comp=2,colors=colors)
p_basil_3=plotEigenFunction(resMFPCA,perception="Basil",comp=2,colors=colors)
grid.arrange(p_lemon_3+ggtitle("a. Lemon"),p_salty_3+ggtitle("b. Salty"),p_sweet_3+ggtitle("c. Sweet"),p_basil_3+ggtitle("d. Basil"),p_acid_3+ggtitle("e. Acid"),nrow=3)

#==================
# PLSDA and product discrimination
#==================
library(RGCCA)
cumsum(eigen_values[,"lambda"])
pc[,-which(colnames(pc)%in%c("Signal","ind"))]
n_comp=9
pc_pls=pc[,1:n_comp]
res_manova=manova(as.matrix(pc[,1:n_comp])~pc[,"Signal"])
library(Hotelling)
ht_46=hotelling.test(as.matrix(pc[pc[,"Signal"]=="S04",1:n_comp]),as.matrix(pc[pc[,"Signal"]=="S06",1:n_comp]))
ht_47=hotelling.test(as.matrix(pc[pc[,"Signal"]=="S04",1:n_comp]),as.matrix(pc[pc[,"Signal"]=="S07",1:n_comp]))
ht_67=hotelling.test(as.matrix(pc[pc[,"Signal"]=="S06",1:n_comp]),as.matrix(pc[pc[,"Signal"]=="S07",1:n_comp]))

summary(res_manova)
pc[,1:n_comp]
respls_opt=rgcca(blocks=list(pc=pc_pls,product=pc[,"Signal"]),method="pls",response=2,ncomp=2)
p_pls=plot(respls_opt,type="samples",block=1,show_sample_names=FALSE)+scale_color_manual(values=c("S06"="#00BA38","S07"="#619CFF","S04"="#F8766D"))+scale_shape_manual(values = c("S04"=16, "S06"=17, "S07"=18))
#plot(respls_opt,type="weights")
#plot(respls_opt,type="weights",comp=2)
p_var1_pls=plotPlsVar(respls_opt,resMFPCA,dim_sel=1)
p_var2_pls=plotPlsVar(respls_opt,resMFPCA,dim_sel=2)
pred_opt=rgcca_predict(respls_opt)$score
res_perm=plsPerm(pc_pls,product=pc[,"Signal"],nperm=1000)

grid.arrange(res_perm$p + ggtitle("a. Histogram of accuracies of permuted datasets"),
             p_pls+ggtitle("b. PLSDA: samples space"),
             p_var1_pls+ggtitle("c. PLSDA: First eigenfunction"),
             p_var2_pls+ggtitle("d. PLSDA: Second eigen function"),nrow=2)
res_perm$p
res_perm$pval
res_perm_S04S06=plsPerm(pc_pls,selec=c("S04","S06"),
                 product=pc[,"Signal"],nperm=100)
res_perm_S04S07=plsPerm(pc_pls,selec=c("S04","S07"),
                        product=pc[,"Signal"],nperm=100)
res_perm_S06S07=plsPerm(pc_pls,selec=c("S07","S06"),
                        product=pc[,"Signal"],nperm=100)
res_perm_S04S06$pval
res_perm_S06S07$pval
res_perm_S04S07$pval
p_pls46=plot(res_perm_S04S06$pls,type="samples",block=1,show_sample_names=FALSE)+scale_color_manual(values=c("S06"="#00BA38","S07"="#619CFF","S04"="#F8766D"))+scale_shape_manual(values = c("S04"=16, "S06"=17, "S07"=18))
p_pls47=plot(res_perm_S04S07$pls,type="samples",block=1,show_sample_names=FALSE)+scale_color_manual(values=c("S06"="#00BA38","S07"="#619CFF","S04"="#F8766D"))+scale_shape_manual(values = c("S04"=16, "S06"=17, "S07"=18))
p_pls67=plot(res_perm_S06S07$pls,type="samples",block=1,show_sample_names=FALSE)+scale_color_manual(values=c("S06"="#00BA38","S07"="#619CFF","S04"="#F8766D"))+scale_shape_manual(values = c("S04"=16, "S06"=17, "S07"=18))

p_var46=plotPlsVar(res_perm_S04S06$pls,resMFPCA,dim_sel=1)
p_var47=plotPlsVar(res_perm_S04S07$pls,resMFPCA,dim_sel=1)
p_var67=plotPlsVar(res_perm_S06S07$pls,resMFPCA,dim_sel=1)
grid.arrange(p_pls46 +ggtitle("PLS P04-P06"),
             p_pls47+ggtitle("PLS P04-P07"),
             p_pls67+ggtitle("PLS P06-P07"),
             p_var46+ggtitle("var P04-P06"),
             p_var47+ggtitle("var P04-P07"),
             p_var67+ggtitle("var P06-P07"),nrow=2)


#==================
# Clustering
#===================
# MFPCA classique
mat_for_dist=pc[,1:n_comp]
distance=dist(mat_for_dist)
reshclust=hclust(distance,method="ward.D2")
plot(reshclust)
gp=as.character(cutree(reshclust,3))
table(gp,pc[,"Signal"])


# Total MFPCA on the graph
#=============================
df_total=dfnorm
df_total[df_total[,"product"]=="S06","time"]=1+dfnorm[dfnorm[,"product"]=="S06","time"]
df_total[df_total[,"product"]=="S07","time"]=2+dfnorm[dfnorm[,"product"]=="S07","time"]
df_total[,"id"]=df_total[,"subject"]
df_total_o=df_total[order(df_total[,"time"]),]
k=30
n.points.tot=900
times.out.tot=seq(0,max(df_total[,"time"]),length.out=n.points.tot)

df_multiFun_tot=prepareForMFPCA(df_total_o,times.out=times.out.tot)
plot(df_multiFun_tot[["Acid"]],main="Acid",xlab="Time")

# Use the MFPCA with a basis of k=8 splines of order 0
knots=c(seq(0,3,length.out=30),1,2)
knots2=sort(knots)
knots3=knots2[!knots2%in%c(0,3)]
length(knots3)
basis_0=createMyBasis(argvals=times.out.tot,knots=knots3,order=0)
matplot(  times.out.tot,  t(basis_0@X),  type = "l")

knots_per0=c(seq(0,1,length.out=10))
knots_per=knots_per0[!knots_per0%in%c(0,1)]
times.out.per=times.out.tot[times.out.tot>=0&times.out.tot<=1]
basis_per=createMyBasis(argvals=times.out.per,knots=knots_per,order=2)
matplot(  times.out.per,  t(basis_per@X),  type = "l")
dim(basis_per@X)
n_basis_per=dim(basis_per@X)[1]
custom_mat=matrix(0,n_basis_per*3,length(times.out.tot))
custom_mat[1:n_basis_per,1:length(times.out.per)]=basis_per@X
custom_mat[(n_basis_per+1):(2*n_basis_per),length(times.out.per)+(1:length(times.out.per))]=basis_per@X
custom_mat[(2*n_basis_per+1):(3*n_basis_per),2*length(times.out.per)+(1:length(times.out.per))]=basis_per@X
myBasis=createMyBasis(argvals=times.out.tot,type="custom", custom_mat=t(custom_mat))

matplot(  times.out.tot,  t(myBasis@X),  type = "l")




hist(df_total_o[!df_total_o[,"time"]%in%c(0,1,2,3),"time"],xlim=c(0,3),breaks=100)
resMFPCA_tot <- MFPCA(df_multiFun_tot, M = 75, uniExpansions = list(
  list(type = "given", functions=basis_0),
  list(type = "given", functions=basis_0),
  list(type = "given", functions=basis_0),
  list(type = "given", functions=basis_0),
  list(type = "given", functions=basis_0),
  list(type = "given", functions=basis_0),
  list(type = "given", functions=basis_0),
  list(type = "given",functions=basis_0)
),fit = TRUE)

resMFPCA_tot <- MFPCA(df_multiFun_tot, M = 75, uniExpansions = list(
  list(type = "given", functions=myBasis),
  list(type = "given", functions=myBasis),
  list(type = "given", functions=myBasis),
  list(type = "given", functions=myBasis),
  list(type = "given", functions=myBasis),
  list(type = "given", functions=myBasis),
  list(type = "given", functions=myBasis),
  list(type = "given",functions=myBasis)
),fit = TRUE)

# calculate reconstruction, too
pc_tot=resMFPCA_tot[["scores"]]
pc_tot=as.data.frame(pc_tot)

colnames(pc_tot)=paste0("pc",1:ncol(pc_tot))
pc_tot[,"ind"]=substr(rownames(pc_tot),7,8)
gg_ind_tot <- ggplot(data=pc_tot,aes(x=pc1,y=pc2,label=ind))+geom_text()+theme_bw()+ggtitle("Principal components scores ", paste0(signal,"-",signal2,"-",signal3))
gg_ind_tot
p_screeplot_tot=screeplotMFPCA(resMFPCA_tot,returnValue=F)

scree_tot=screeplotMFPCA(resMFPCA_tot,returnValue=T)
cumsum(scree_tot[,"lambda"]) # 18
p_avg_tot=plotAvgMFPCA(resMFPCA_tot)+ylim(0,1)
p_contrib_tot=plotContribMFPCA(resMFPCA_tot)
p_v1_tot=plotAllEigenFunction(resMFPCA_tot,colors=colors)
p_v2_tot=plotAllEigenFunction(resMFPCA_tot,comp=2,colors=colors)

grid.arrange(gg_ind_tot+ggtitle("a. MFPCA on entiere sequence"),p_v1_tot+ggtitle("b. First eigenfunction"),p_avg_tot+ggtitle("c. Frequence of citations"),p_screeplot_tot+ggtitle("d. Screeplot")+ylim(0,0.12))
# outliers: 42, 24, 27,37,35
outlier="42"
ggplot(df_total[df_total[,"subject"]==paste0("TCATA_",outlier),],aes(x=time,y=score,color=descriptor))+geom_point()+scale_color_manual(values=colors)+theme_bw()
#
n_comp_tot=7
mat_for_dist_tot=pc_tot[,1:n_comp_tot]
rownames(mat_for_dist_tot)=substr(rownames(pc_tot),7,8)
distance_tot=dist(mat_for_dist_tot,method="euclidea")
reshclust_tot=hclust(distance_tot,method="ward.D2")
plot(reshclust_tot)
groups_tot=as.character(cutree(reshclust_tot,3))
names(groups_tot)=rownames(pc_tot)
gp1=names(groups_tot[groups_tot==1])
gp2=names(groups_tot[groups_tot==2])
gp3=names(groups_tot[groups_tot==3])

df4_1=dfnorm[dfnorm[,"product"]=="S04"&dfnorm[,"subject"]%in%gp1,]
df6_1=dfnorm[dfnorm[,"product"]=="S06"&dfnorm[,"subject"]%in%gp1,]
df7_1=dfnorm[dfnorm[,"product"]=="S07"&dfnorm[,"subject"]%in%gp1,]
df4_2=dfnorm[dfnorm[,"product"]=="S04"&dfnorm[,"subject"]%in%gp2,]
df6_2=dfnorm[dfnorm[,"product"]=="S06"&dfnorm[,"subject"]%in%gp2,]
df7_2=dfnorm[dfnorm[,"product"]=="S07"&dfnorm[,"subject"]%in%gp2,]
df4_3=dfnorm[dfnorm[,"product"]=="S04"&dfnorm[,"subject"]%in%gp3,]
df6_3=dfnorm[dfnorm[,"product"]=="S06"&dfnorm[,"subject"]%in%gp3,]
df7_3=dfnorm[dfnorm[,"product"]=="S07"&dfnorm[,"subject"]%in%gp3,]
df4_4=dfnorm[dfnorm[,"product"]=="S04"&dfnorm[,"subject"]%in%gp4,]
df6_4=dfnorm[dfnorm[,"product"]=="S06"&dfnorm[,"subject"]%in%gp4,]
df7_4=dfnorm[dfnorm[,"product"]=="S07"&dfnorm[,"subject"]%in%gp4,]

p_S04_1=tcataCurve(df4_1,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S06_1=tcataCurve(df6_1,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S07_1=tcataCurve(df7_1,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S04_2=tcataCurve(df4_2,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S06_2=tcataCurve(df6_2,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S07_2=tcataCurve(df7_2,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S04_3=tcataCurve(df4_3,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S06_3=tcataCurve(df6_3,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")
p_S07_3=tcataCurve(df7_3,times=seq(0, 1, by = 0.05),colors=colors)+theme(legend.position="none")

grid.arrange(
             p_S04_1+ggtitle("b. Group 1: P04"),
             p_S06_1+ggtitle("c. Group 1: P06"),
             p_S07_1+ggtitle("d. Group 1: P07"),
             p_S04_2+ggtitle("b. Group 2: P04"),
             p_S06_2+ggtitle("c. Group 2: P06"),
             p_S07_2+ggtitle("d. Group 2: P07"),
             p_S04_3+ggtitle("b. Group 3: P04"),
             p_S06_3+ggtitle("c. Group 3: P06"),
             p_S07_3+ggtitle("d. Group 3: P07"),
             nrow=3)
summary(factor(groups_tot))
groups_tot[groups_tot==3]
# using oclust
library(oclust)
library(isotree)
# using isofree
model <- isolation.forest(pc_tot[,1:n_comp_tot], ndim=1, ntrees=1000, nthreads=1)
predict_if <- predict(model, pc_tot, type="avg_depth")# scores ou "avg_depth" (plus petites valeurs)
hist(predict_if)

outliers_predict <- predict_if  < 6
outliers=names(outliers_predict[outliers_predict])
predict_if[outliers_predict]

df_total_o[df_total_o[,"id"]==outliers[1]&!df_total_o[,"time"]%in%c(0,1,2,3),]
df_total_o[df_total_o[,"id"]==outliers[2]&!df_total_o[,"time"]%in%c(0,1,2,3),]
df_total_o[df_total_o[,"id"]==outliers[3]&!df_total_o[,"time"]%in%c(0,1,2,3),]
df_total_o[df_total_o[,"id"]==outliers[4]&!df_total_o[,"time"]%in%c(0,1,2,3),]
df_total_o[df_total_o[,"id"]==outliers[5]&!df_total_o[,"time"]%in%c(0,1,2,3),]
df_total_o[df_total_o[,"id"]==outliers[6]&!df_total_o[,"time"]%in%c(0,1,2,3),]
df_total_o[df_total_o[,"id"]==outliers[7]&!df_total_o[,"time"]%in%c(0,1,2,3),]


tcataA$df[tcataA$df[,"id"]%in%c("TCATA_42_S06_TCATA")&!tcataA$df[,"time"]%in%c(0,1),]


library(dplyr)

library(ggplot2)

scores_pca_duration <- as.data.frame(res.pca$x)
scores_pca_duration$product=product2
ggplot(scores_pca_duration, aes(PC1, PC2, colour = product)) +
  geom_point(size = 3) + ggtitle("PCA of durations")+theme_bw()
scores[which.max(scores[,1]),]
tcataA$df[tcataA$df[,"id"]%in%c(scores_pca_duration)&!tcataA$df[,"time"]%in%c(0,1),]
model_pca_dur <- isolation.forest(scores_pca_duration, ndim=1, ntrees=1000, nthreads=1)
predict_if_pca_dur <- predict(model_pca_dur,scores_pca_duration, type="avg_depth")# scores ou "avg_depth" (plus petites valeurs)
plot(predict_if_pca_dur)
hist(predict_if_pca_dur)
outliers_predict_pca_dur <- predict_if_pca_dur  < 9
outliers=names(outliers_predict_pca_dur[outliers_predict_pca_dur])
predict_if_pca_dur[outliers_predict_pca_dur]
outlier_pca=rownames(scores[scores[,1]>2,])
