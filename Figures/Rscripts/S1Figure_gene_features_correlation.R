setwd("/home/jacobhsu")
data <- read.csv("./ISPP_Paper/Data/S1.csv",header=T)
data$gene = NULL
library(corrplot)
tiff("./ISPP_Paper/Figures/S1Figure.tiff", width = 178, height = 200, units = 'mm',pointsize = 8, res = 100)
par(mai=c(10.6,1.52,1.25,2.1))
cor_data=cor(data,use="complete.obs",method="spearman")
corrplot(cor_data, method="number", tl.pos="lt", type="upper",tl.col="black", tl.cex=1.0, tl.srt=45,addCoef.col="black", addCoefasPercent = TRUE,p.mat = 1-abs(cor_data), sig.level=0.7, insig = "blank")  
write.csv(cor_data, file = "./ISPP_Paper/Data/S2_test.csv")

dev.off()
#cor(data, use="complete.obs", method="kendall")
#library(corrgram)
#corrgram(data, order=TRUE, lower.panel=panel.shade,upper.panel=panel.pie, text.panel=panel.txt)
#corrgram(data, order=TRUE, lower.panel=panel.ellipse,upper.panel=panel.pts, text.panel=panel.txt,diag.panel=panel.minmax)
#corrgram(data, order=NULL, lower.panel=panel.shade,upper.panel=NULL, text.panel=panel.txt,main="Features_matrix")
