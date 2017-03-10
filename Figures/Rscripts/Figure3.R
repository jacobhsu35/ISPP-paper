setwd("/home/jacobhsu/")
tiff("./ISPP_Paper/Figures/Figure3.tiff", width = 86, height = 180, units = 'mm',pointsize = 10, res = 900)
par(mfrow=c(3,1),family="mono")
par(mai=c(0.6,0.52,0.25,0.1)) #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
library(ROCR);                                                                                        

#### Path
##########################    CGD_AD  ####################################################
CGD_AD_pos <- unique(read.table("./ISPP_Paper/Data/AD/AD_genelist.txt")[,1]);
CGD_AD_neg <- unique(read.table("./ISPP_Paper/Data/AD/control_geneset.txt")[,1]);
###############################################################################
CONS <- read.table("./ISPP_Paper/Data/For_ROC/Approved_fordist_constraint_scores.txt",header=T,sep="\t",fill=T);
CONS$IsCase <- as.numeric(CONS$Gene %in% CGD_AD_pos)
predCONS <- prediction(CONS$score*1, CONS$IsCase)
perfCONS <- performance(predCONS, 'tpr','fpr')
plot(perfCONS,col=2,main = "CGD_AD (876)",lwd=0.5)
auc_CONS <- performance(predCONS, "auc")
auc.area <- slot(auc_CONS, "y.values")[[1]]
cat(auc.area)
###############################################################################
INTOR <- read.table("./ISPP_Paper/Data/For_ROC/Approved_Genic_score.txt",header=T,fill=T);
INTOR$IsCase <- as.numeric(INTOR$Gene %in% CGD_AD_pos)
predINTOR <- prediction(INTOR$score*-1, INTOR$IsCase)
perfINTOR <- performance(predINTOR, 'tpr','fpr')
plot(perfINTOR,col=3,add=TRUE,lwd=0.5)
auc_INTOR <- performance(predINTOR, "auc")
auc.area <- slot(auc_INTOR, "y.values")[[1]]
cat(auc.area)
##############################################################################
RECE <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_rece_score_forROC.txt",header=T,sep="\t",fill=T);
RECE$IsCase <- as.numeric(RECE$Gene %in% CGD_AD_pos)
predRECE <- prediction(RECE$score*1, RECE$IsCase)
perfRECE <- performance(predRECE, 'tpr','fpr')  
plot(perfRECE,col=4,add=TRUE,lwd=0.5) 
auc_RECE <- performance(predRECE, "auc")
auc.area <- slot(auc_RECE, "y.values")[[1]]
cat(auc.area)
##############################################################################
HI <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_gene_HI_forROC.txt",header=T,sep="\t",fill=T);
HI$IsCase <- as.numeric(HI$Gene %in% CGD_AD_pos)
predHI <- prediction(HI$score*1, HI$IsCase)
perfHI <- performance(predHI, 'tpr','fpr')
plot(perfHI,col=5,add=TRUE,lwd=0.5)
auc_HI <- performance(predHI, "auc")
auc.area_HI <- slot(auc_HI, "y.values")[[1]]
cat(auc.area_HI)
##############################################################################
GDI <- read.table("./ISPP_Paper/Data/For_ROC/GDI_gene_score.txt",header=T,sep="\t",fill=T);
GDI$IsCase <- as.numeric(GDI$Gene %in% CGD_AD_pos)
predGDI <- prediction(GDI$GDI_Phred*1, GDI$IsCase)
perfGDI <- performance(predGDI, 'tpr','fpr')
plot(perfGDI,col=8,add=TRUE,lwd=0.5)
auc_GDI <- performance(predGDI, "auc")
auc.area_GDI <- slot(auc_GDI, "y.values")[[1]]
cat(auc.area_GDI)
##############################################################################
Net <- read.table("./ISPP_Paper/Data/For_ROC/network_gene_score.txt",header=T,sep="\t",fill=T);
Net$IsCase <- as.numeric(Net$GENE_NAME %in% CGD_AD_pos)
predNet <- prediction(Net$PREDICTED_INDISPENSABILITY_SCORE*1, Net$IsCase)
perfNet <- performance(predNet, 'tpr','fpr')
plot(perfNet,col=9,add=TRUE,lwd=0.5)
auc_Net <- performance(predNet, "auc")
auc.area_Net <- slot(auc_Net, "y.values")[[1]]
cat(auc.area_Net)
##############################################################################
CGD_AD_M <- read.csv("./ISPP_Paper/Data/AD/merge_AD_forAUC.csv",header=F);
lines(CGD_AD_M$V2,CGD_AD_M$V1,col=6,lwd=0.5)
std <- read.csv("./ISPP_Paper/Data/For_ROC/standard.csv",header=F)
lines(std$V1,std$V2,col=1,lwd=0.5)
legend("bottomright",pch=15,col=c(5,4,3,9,2,8,6),legend = c("HI:      0.69",
                                                            "REC:     0.73",
                                                            "RVIS:    0.64",
                                                            "NET:     0.70",
                                                            "CONS:    0.65",
                                                            "GDI:     0.49",
                                                            "ISPP_AD: 0.75"))
##########################    CGD_AR  ####################################################
CGD_AR_pos <- unique(read.table("./ISPP_Paper/Data/AR/AR_genelist.txt")[,1]);
CGD_AR_neg <- unique(read.table("./ISPP_Paper/Data/AR/control_geneset.txt")[,1]);
###############################################################################
CONS <- read.table("./ISPP_Paper/Data/For_ROC/Approved_fordist_constraint_scores.txt",header=T,sep="\t",fill=T);
CONS$IsCase <- as.numeric(CONS$Gene %in% CGD_AR_pos)
predCONS <- prediction(CONS$score*1, CONS$IsCase)
perfCONS <- performance(predCONS, 'tpr','fpr')
plot(perfCONS,col=2,main = "CGD_AR (1500)",lwd=0.5)
auc_CONS <- performance(predCONS, "auc")
auc.area <- slot(auc_CONS, "y.values")[[1]]
cat(auc.area)
###############################################################################
INTOR <- read.table("./ISPP_Paper/Data/For_ROC/Approved_Genic_score.txt",header=T,fill=T);
INTOR$IsCase <- as.numeric(INTOR$Gene %in% CGD_AR_pos)
predINTOR <- prediction(INTOR$score*-1, INTOR$IsCase)
perfINTOR <- performance(predINTOR, 'tpr','fpr')
plot(perfINTOR,col=3,add=TRUE,lwd=0.5)
auc_INTOR <- performance(predINTOR, "auc")
auc.area <- slot(auc_INTOR, "y.values")[[1]]
cat(auc.area)
##############################################################################
RECE <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_rece_score_forROC.txt",header=T,sep="\t",fill=T);
RECE$IsCase <- as.numeric(RECE$Gene %in% CGD_AR_pos)
predRECE <- prediction(RECE$score*1, RECE$IsCase)
perfRECE <- performance(predRECE, 'tpr','fpr')  
plot(perfRECE,col=4,add=TRUE,lwd=0.5) 
auc_RECE <- performance(predRECE, "auc")
auc.area <- slot(auc_RECE, "y.values")[[1]]
cat(auc.area)
##############################################################################
HI <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_gene_HI_forROC.txt",header=T,sep="\t",fill=T);
HI$IsCase <- as.numeric(HI$Gene %in% CGD_AR_pos)
predHI <- prediction(HI$score*1, HI$IsCase)
perfHI <- performance(predHI, 'tpr','fpr')
plot(perfHI,col=5,add=TRUE,lwd=0.5)
auc_HI <- performance(predHI, "auc")
auc.area_HI <- slot(auc_HI, "y.values")[[1]]
cat(auc.area_HI)
##############################################################################
GDI <- read.table("./ISPP_Paper/Data/For_ROC/GDI_gene_score.txt",header=T,sep="\t",fill=T);
GDI$IsCase <- as.numeric(GDI$Gene %in% CGD_AR_pos)
predGDI <- prediction(GDI$GDI_Phred*1, GDI$IsCase)
perfGDI <- performance(predGDI, 'tpr','fpr')
plot(perfGDI,col=8,add=TRUE,lwd=0.5)
auc_GDI <- performance(predGDI, "auc")
auc.area_GDI <- slot(auc_GDI, "y.values")[[1]]
cat(auc.area_GDI)
##############################################################################
Net <- read.table("./ISPP_Paper/Data/For_ROC/network_gene_score.txt",header=T,sep="\t",fill=T);
Net$IsCase <- as.numeric(Net$GENE_NAME %in% CGD_AR_pos)
predNet <- prediction(Net$PREDICTED_INDISPENSABILITY_SCORE*1, Net$IsCase)
perfNet <- performance(predNet, 'tpr','fpr')
plot(perfNet,col=9,add=TRUE,lwd=0.5)
auc_Net <- performance(predNet, "auc")
auc.area_Net <- slot(auc_Net, "y.values")[[1]]
cat(auc.area_Net)
##############################################################################
CGD_AR_M <- read.csv("./ISPP_Paper/Data/AR/merge_AR_forAUC.csv",header=F);
lines(CGD_AR_M$V2,CGD_AR_M$V1,col=6,lwd=0.5)
std <- read.csv("./ISPP_Paper/Data/For_ROC/standard.csv",header=F)
lines(std$V1,std$V2,col=1,lwd=0.5)
legend("bottomright",pch=15,col=c(5,4,3,9,2,8,6),legend = c("HI:      0.49",
                                                            "REC:     0.65",
                                                            "RVIS:    0.54",
                                                            "NET:     0.62",
                                                            "CONS:    0.43",
                                                            "GDI:     0.58",
                                                            "ISPP_AR: 0.73"))
########################  CGD_XL                   ######################################################
CGD_XL_pos <- unique(read.table("./ISPP_Paper/Data/XL/XL_genelist.txt")[,1]);
CGD_XL_negonly <- unique(read.table("./ISPP_Paper/Data/XL/control_geneset.txt")[,1]);
###############################################################################
CONS <- read.table("./ISPP_Paper/Data/XL/X_chr_Cons_forAUC.txt",header=T,sep="\t",fill=T);
CONS$IsCase <- as.numeric(CONS$Gene %in% CGD_XL_pos)
predCONS <- prediction(CONS$score*1, CONS$IsCase)
perfCONS <- performance(predCONS, 'tpr','fpr')
plot(perfCONS,col=2,main="CGD_XL (182)",lwd=0.5)
auc_CONS <- performance(predCONS, "auc")
auc.area <- slot(auc_CONS, "y.values")[[1]]
cat(auc.area)
###############################################################################
INTOR <- read.table("./ISPP_Paper/Data/XL/X_chr_Genic_forAUC.txt",header=T,fill=T);
INTOR$IsCase <- as.numeric(INTOR$Gene %in% CGD_XL_pos)
predINTOR <- prediction(INTOR$score*-1, INTOR$IsCase)
perfINTOR <- performance(predINTOR, 'tpr','fpr')
plot(perfINTOR,col=3,add=TRUE,lwd=0.5)
auc_INTOR <- performance(predINTOR, "auc")
auc.area <- slot(auc_INTOR, "y.values")[[1]]
cat(auc.area)
##############################################################################
RECE <- read.table("./ISPP_Paper/Data/XL/X_chr_Rece_forAUC.txt",header=T,sep="\t",fill=T);
RECE$IsCase <- as.numeric(RECE$Gene %in% CGD_XL_pos)
predRECE <- prediction(RECE$score*1, RECE$IsCase)
perfRECE <- performance(predRECE, 'tpr','fpr')  
plot(perfRECE,col=4,add=TRUE,lwd=0.5) 
auc_RECE <- performance(predRECE, "auc")
auc.area <- slot(auc_RECE, "y.values")[[1]]
cat(auc.area)
##############################################################################
HI <- read.table("./ISPP_Paper/Data/XL/X_chr_HI_forAUC.txt",header=T,sep="\t",fill=T);
HI$IsCase <- as.numeric(HI$Gene %in% CGD_XL_pos)
predHI <- prediction(HI$score*1, HI$IsCase)
perfHI <- performance(predHI, 'tpr','fpr')
plot(perfHI,col=5,add=TRUE,lwd=0.5)
auc_HI <- performance(predHI, "auc")
auc.area_HI <- slot(auc_HI, "y.values")[[1]]
cat(auc.area_HI)
##############################################################################
GDI <- read.table("./ISPP_Paper/Data/XL/X_chr_GDI_forAUC.txt",header=T,sep="\t",fill=T);
GDI$IsCase <- as.numeric(GDI$Gene %in% CGD_XL_pos)
predGDI <- prediction(GDI$GDI_Phred*1, GDI$IsCase)
perfGDI <- performance(predGDI, 'tpr','fpr')
plot(perfGDI,col=8,add=TRUE,lwd=0.5)
auc_GDI <- performance(predGDI, "auc")
auc.area_GDI <- slot(auc_GDI, "y.values")[[1]]
cat(auc.area_GDI)
##############################################################################
Net <- read.table("./ISPP_Paper/Data/XL/X_chr_Net_forAUC.txt",header=T,sep="\t",fill=T);
Net$IsCase <- as.numeric(Net$GENE_NAME %in% CGD_XL_pos)
predNet <- prediction(Net$PREDICTED_INDISPENSABILITY_SCORE*1, Net$IsCase)
perfNet <- performance(predNet, 'tpr','fpr')
plot(perfNet,col=9,add=TRUE,lwd=0.5)
auc_Net <- performance(predNet, "auc")
auc.area_Net <- slot(auc_Net, "y.values")[[1]]
cat(auc.area_Net)
##############################################################################
CGD_XL_M <- read.csv("./ISPP_Paper/Data/XL/CGD_XL_auc_onlyX_chr_genes.csv",header=F);
lines(CGD_XL_M$V2,CGD_XL_M$V1,col=6,lwd=0.5)
std <- read.csv("./ISPP_Paper/Data/For_ROC/standard.csv",header=F)
lines(std$V1,std$V2,col=1,lwd=0.5)
legend("bottomright",pch=15,col=c(5,4,3,9,2,8,6),legend = c("HI:      0.69",
                                                            "REC:     0.76",
                                                            "RVIS:    0.60",
                                                            "NET:     0.66",
                                                            "CONS:    0.64",
                                                            "GDI:     0.49",
                                                            "ISPP_XL: 0.85"))
#############################################################################
dev.off()
