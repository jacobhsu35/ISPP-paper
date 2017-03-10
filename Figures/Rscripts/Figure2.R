

#### Path
#wd <- getwd()
setwd("/home/jacobhsu/")
tiff("./ISPP_Paper/Figures/Figure2.tiff", width = 178, height = 190, units = 'mm',pointsize = 13.5, res = 300)
par(mai=c(0.6,0.52,0.25,0.1))
#setwd("/home/jacobhsu")
library(gplots)
library(RColorBrewer)

data <- read.csv("./ISPP_Paper/Data/S3.csv")
genelist <- data[,1]  # name of each row
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- genelist
my_palette <- colorRampPalette(c("darkgreen", "white", "darkred"))(n = 20)
col_breaks = c(seq(-35,-3.9011,length=10), # for green z-score < -3.874514 (0.05/(18*26))  Bonferroni correction
           +   seq(-3.901,3.9011,length=1), # for yellow
           +   seq(3.90111,35,length=10)) # for red z-score > 3.874514 (0.05/(18*26))  Bonferroni correctio

### New color key
lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(1.5,4)
lhei = c(0.8,6,1.5)
par(cex.main=0.8)
heatmap.2(mat_data,
          #cellnote = mat_data,  # same data set for cell labels
          main = " Z score (P < 9.615E-5, Bonferroni correction) ", # heat map title
          key = TRUE,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(7,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",            # turn off column clustering
          lmat = lmat, lwid = lwid, lhei = lhei)
dev.off()

