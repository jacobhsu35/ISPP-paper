setwd("/home/jacobhsu/")
tiff("./ISPP_Paper/Figures/S3_Figure.tiff", width = 200, height = 150, units = 'mm',pointsize = 10, res = 100)
par(mfrow=c(1,2),family="mono")
par(mai=c(1.22,1.02,0.82,0.42))
library(ROCR);                                                                                        

##########################    hOMIM_earlyOnset  ####################################################
h_earlyOnset_AD_pos <- unique(read.table("./ISPP_Paper/Data/hOMIM_early/Approved_hOMIM_earlyOnset_list_forROC.txt")[,1]);
h_earlyOnset_AD_neg <- unique(read.table("./ISPP_Paper/Data/hOMIM_early/control_geneset.txt")[,1]);
###############################################################################
CONS <- read.table("./ISPP_Paper/Data/For_ROC/Approved_fordist_constraint_scores.txt",header=T,sep="\t",fill=T);
CONS$IsCase <- as.numeric(CONS$Gene %in% h_earlyOnset_AD_pos)
predCONS <- prediction(CONS$score*1, CONS$IsCase)
perfCONS <- performance(predCONS, 'tpr','fpr')
plot(perfCONS,col=2,main = "hOMIM_earlyOnset (646)")
auc_CONS <- performance(predCONS, "auc")
auc.area <- slot(auc_CONS, "y.values")[[1]]
cat(auc.area)
###############################################################################
INTOR <- read.table("./ISPP_Paper/Data/For_ROC/Approved_Genic_score.txt",header=T,fill=T);
INTOR$IsCase <- as.numeric(INTOR$Gene %in% h_earlyOnset_AD_pos)
predINTOR <- prediction(INTOR$score*-1, INTOR$IsCase)
perfINTOR <- performance(predINTOR, 'tpr','fpr')
plot(perfINTOR,col=3,add=TRUE)
auc_INTOR <- performance(predINTOR, "auc")
auc.area <- slot(auc_INTOR, "y.values")[[1]]
cat(auc.area)
##############################################################################
RECE <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_rece_score_forROC.txt",header=T,sep="\t",fill=T);
RECE$IsCase <- as.numeric(RECE$Gene %in% h_earlyOnset_AD_pos)
predRECE <- prediction(RECE$score*1, RECE$IsCase)
perfRECE <- performance(predRECE, 'tpr','fpr')  
plot(perfRECE,col=4,add=TRUE) 
auc_RECE <- performance(predRECE, "auc")
auc.area <- slot(auc_RECE, "y.values")[[1]]
cat(auc.area)
##############################################################################
HI <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_gene_HI_forROC.txt",header=T,sep="\t",fill=T);
HI$IsCase <- as.numeric(HI$Gene %in% h_earlyOnset_AD_pos)
predHI <- prediction(HI$score*1, HI$IsCase)
perfHI <- performance(predHI, 'tpr','fpr')
plot(perfHI,col=5,add=TRUE)
auc_HI <- performance(predHI, "auc")
auc.area_HI <- slot(auc_HI, "y.values")[[1]]
cat(auc.area_HI)
##############################################################################
GDI <- read.table("./ISPP_Paper/Data/For_ROC/GDI_gene_score.txt",header=T,sep="\t",fill=T);
GDI$IsCase <- as.numeric(GDI$Gene %in% h_earlyOnset_AD_pos)
predGDI <- prediction(GDI$GDI_Phred*1, GDI$IsCase)
perfGDI <- performance(predGDI, 'tpr','fpr')
plot(perfGDI,col=8,add=TRUE,lwd=0.5)
auc_GDI <- performance(predGDI, "auc")
auc.area_GDI <- slot(auc_GDI, "y.values")[[1]]
cat(auc.area_GDI)
##############################################################################
Net <- read.table("./ISPP_Paper/Data/For_ROC/network_gene_score.txt",header=T,sep="\t",fill=T);
Net$IsCase <- as.numeric(Net$GENE_NAME %in% h_earlyOnset_AD_pos)
predNet <- prediction(Net$PREDICTED_INDISPENSABILITY_SCORE*1, Net$IsCase)
perfNet <- performance(predNet, 'tpr','fpr')
plot(perfNet,col=9,add=TRUE,lwd=0.5)
auc_Net <- performance(predNet, "auc")
auc.area_Net <- slot(auc_Net, "y.values")[[1]]
cat(auc.area_Net)
##############################################################################
h_earlyOnset_AD_M <- read.csv("./ISPP_Paper/Data/hOMIM_early/merge_hOMIM_early_forAUC.csv",header=F);
lines(h_earlyOnset_AD_M$V2,h_earlyOnset_AD_M$V1,col=6)
std <- read.csv("./ISPP_Paper/Data/For_ROC/standard.csv",header=F)
lines(std$V1,std$V2,col=1)
legend("bottomright",pch=15,col=c(5,4,3,9,2,8,6),legend = c("HI:     0.59",
                                                            "REC:    0.81",
                                                            "RVIS:   0.58",
                                                            "NET:    0.69",
                                                            "CONS:   0.51",
                                                            "GDI:    0.57",
                                                            "ISPP:   0.82"))
############################    CGD_PAE          ##################################################
CGD_PAE_pos <- unique(read.table("./ISPP_Paper/Data/PAE/Pediatric_genelist.txt")[,1]);
CGD_PAE_neg <- unique(read.table("./ISPP_Paper/Data/PAE/control_geneset.txt")[,1]);
###############################################################################
CONS <- read.table("./ISPP_Paper/Data/For_ROC/Approved_fordist_constraint_scores.txt",header=T,sep="\t",fill=T);
CONS$IsCase <- as.numeric(CONS$Gene %in% CGD_PAE_pos)
predCONS <- prediction(CONS$score*1, CONS$IsCase)
perfCONS <- performance(predCONS, 'tpr','fpr')
plot(perfCONS,col=2,main="CGD_Paediatric (1401)")
auc_CONS <- performance(predCONS, "auc")
auc.area <- slot(auc_CONS, "y.values")[[1]]
cat(auc.area)
###############################################################################
INTOR <- read.table("./ISPP_Paper/Data/For_ROC/Approved_Genic_score.txt",header=T,fill=T);
INTOR$IsCase <- as.numeric(INTOR$Gene %in% CGD_PAE_pos)
predINTOR <- prediction(INTOR$score*-1, INTOR$IsCase)
perfINTOR <- performance(predINTOR, 'tpr','fpr')
plot(perfINTOR,col=3,add=TRUE)
auc_INTOR <- performance(predINTOR, "auc")
auc.area <- slot(auc_INTOR, "y.values")[[1]]
cat(auc.area)
##############################################################################
RECE <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_rece_score_forROC.txt",header=T,sep="\t",fill=T);
RECE$IsCase <- as.numeric(RECE$Gene %in% CGD_PAE_pos)
predRECE <- prediction(RECE$score*1, RECE$IsCase)
perfRECE <- performance(predRECE, 'tpr','fpr')  
plot(perfRECE,col=4,add=TRUE) 
auc_RECE <- performance(predRECE, "auc")
auc.area <- slot(auc_RECE, "y.values")[[1]]
cat(auc.area)
##############################################################################
HI <- read.table("./ISPP_Paper/Data/For_ROC/Approved_dbNSFP2.4_gene_HI_forROC.txt",header=T,sep="\t",fill=T);
HI$IsCase <- as.numeric(HI$Gene %in% CGD_PAE_pos)
predHI <- prediction(HI$score*1, HI$IsCase)
perfHI <- performance(predHI, 'tpr','fpr')
plot(perfHI,col=5,add=TRUE)
auc_HI <- performance(predHI, "auc")
auc.area_HI <- slot(auc_HI, "y.values")[[1]]
cat(auc.area_HI)
##############################################################################
GDI <- read.table("./ISPP_Paper/Data/For_ROC/GDI_gene_score.txt",header=T,sep="\t",fill=T);
GDI$IsCase <- as.numeric(GDI$Gene %in% CGD_PAE_pos)
predGDI <- prediction(GDI$GDI_Phred*1, GDI$IsCase)
perfGDI <- performance(predGDI, 'tpr','fpr')
plot(perfGDI,col=8,add=TRUE,lwd=0.5)
auc_GDI <- performance(predGDI, "auc")
auc.area_GDI <- slot(auc_GDI, "y.values")[[1]]
cat(auc.area_GDI)
##############################################################################
Net <- read.table("./ISPP_Paper/Data/For_ROC/network_gene_score.txt",header=T,sep="\t",fill=T);
Net$IsCase <- as.numeric(Net$GENE_NAME %in% CGD_PAE_pos)
predNet <- prediction(Net$PREDICTED_INDISPENSABILITY_SCORE*1, Net$IsCase)
perfNet <- performance(predNet, 'tpr','fpr')
plot(perfNet,col=9,add=TRUE,lwd=0.5)
auc_Net <- performance(predNet, "auc")
auc.area_Net <- slot(auc_Net, "y.values")[[1]]
cat(auc.area_Net)
##############################################################################
CGD_PAE_M <- read.csv("./ISPP_Paper/Data/PAE/merge_PAE_forAUC.csv",header=F);
lines(CGD_PAE_M$V2,CGD_PAE_M$V1,col=6)
std <- read.csv("./ISPP_Paper/Data//For_ROC/standard.csv",header=F)
lines(std$V1,std$V2,col=1)
legend("bottomright",pch=15,col=c(5,4,3,9,2,8,6),legend = c("HI:     0.57",
                                                            "REC:    0.76",
                                                            "RVIS:   0.56",
                                                            "NET:    0.66",
                                                            "CONS:   0.51",
                                                            "GDI:    0.54",
                                                            "ISPP:   0.76"))
dev.off()

