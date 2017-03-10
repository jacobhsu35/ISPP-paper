CGD_AD <- read.csv("/home/jacobhsu/ISPP_Paper/Data/REC_compare/CGD_AD_trainingset_all_scores_forR.csv", header = T)
CGD_AR <- read.csv("/home/jacobhsu/ISPP_Paper/Data/REC_compare/CGD_AR_trainingset_all_scores_forR.csv", header = T)
h_AD <- read.csv("/home/jacobhsu/ISPP_Paper/Data/REC_compare/hOMIM_AD_all_scores_forR.csv", header = T)
h_AR <- read.csv("/home/jacobhsu/ISPP_Paper/Data/REC_compare/hOMIM_AR_all_scores_forR.csv", header = T)
Genic_AD <- read.csv("/home/jacobhsu/ISPP_Paper/Data/REC_compare/Genic_AD_trainingset_all_scores_forR.csv",header = T)
Genic_AR <- read.csv("/home/jacobhsu/ISPP_Paper/Data/REC_compare/Genic_AR_trainingset_all_scores_forR.csv",header = T)

par(mfrow=c(3,1))
par(mai=c(1.02,0.82,0.82,0.42))

wilcox.test(CGD_AD$rec,CGD_AR$rec)
wilcox.test(h_AD$rec,h_AR$rec)
wilcox.test(Genic_AD$rec,Genic_AR$rec)

boxplot(CGD_AD$rec,CGD_AR$rec, main = "Recessive score (REC) on CGD genelists, (p-value = 1.554e-07, Wilcoxon rank sum test)", names=c("CGD_AD (867)", "CGD_AR (1480)"))
boxplot(h_AD$rec,h_AR$rec, main = "Recessive score (REC) on hOMIM genelists, (p-value = 0.06134, Wilcoxon rank sum test)", names=c("h_AD (418)", "h_AR (569)"))
boxplot(Genic_AD$rec,Genic_AR$rec,main = "Recessive score (REC) on Genic Intolerance genelists, (p-value = 7.321e-10, Wilcoxon rank sum test)", names=c("Genic_AD (362)", "Genic_AR (817)"))

tiff("/home/jacobhsu/ISPP_Paper/Figures/S2Figure.tiff", width = 178, height = 190, units = 'mm',pointsize = 13.5, res = 150)
par(mfrow=c(3,1))
#par(mai=c(1.02,0.82,0.82,0.42))
boxplot(CGD_AD$rec,CGD_AR$rec, main = "Recessive score (REC) on CGD gene lists, (p-value = 1.554e-07, Wilcoxon rank sum test)", names=c("CGD_AD (867)", "CGD_AR (1480)"))
boxplot(h_AD$rec,h_AR$rec, main = "Recessive score (REC) on hOMIM gene lists, (p-value = 0.06134, Wilcoxon rank sum test)", names=c("h_AD (418)", "h_AR (569)"))
boxplot(Genic_AD$rec,Genic_AR$rec,main = "Recessive score (REC) on RVIS gene lists, (p-value = 7.321e-10, Wilcoxon rank sum test)", names=c("Genic_AD (362)", "Genic_AR (817)"))

dev.off()  
