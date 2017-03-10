library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)

#par(mfrow=c(3,1))
par(mai=c(1.02,0.82,0.82,0.42))

overlap <- draw.quintuple.venn(area1=876, area2=1500, area3=182, area4=1401, area5=643, n12=0, n13=0, n14=348, n15=137,
                               n23=0, n24=736, n25=305, n34=74, n35=51, n45=410, n123=0, n124=0, n125=0, n134=0,
                               n135=0, n145=77, n234=0, n235=0, n245=206, n345=33, n1234=0, n1235=0,
                               n1245=0, n1345=0, n2345=0, n12345=0, category = c("CGD_AD","CGD_AR","CGD_XL","CGD_PAE","hOMIM_early"),
                               lty = "blank", fill = c("skyblue", "pink1", "mediumorchid","orange","green"),euler.d=TRUE,scaled = TRUE,margin = 0.05)




library(gridExtra)
#grid.arrange(gTree(children=hOMIM_AR), gTree(children=REC_AR),gTree(children=GENIC_AR), ncol=3 )
#grid.arrange(gTree(children=hOMIM_AR), gTree(children=REC_AR),gTree(children=GENIC_AR), ncol=3, top= "The overlapping between combined models (AD,AR,XL) and respective recessive genes in three published studies")
setwd("/home/jacobhsu/")
tiff("./ISPP_Paper/Figures/S5_150dpi.tiff", width = 86, height = 70, units = 'mm',pointsize = 3, res = 150)
#grid.arrange(gTree(children=hOMIM_AR), gTree(children=REC_AR),gTree(children=GENIC_AR), ncol=3 )
overlap <- draw.quintuple.venn(area1=876, area2=1500, area3=182, area4=1401, area5=643, n12=0, n13=0, n14=348, n15=137,
                               n23=0, n24=736, n25=305, n34=74, n35=51, n45=410, n123=0, n124=0, n125=0, n134=0,
                               n135=0, n145=77, n234=0, n235=0, n245=206, n345=33, n1234=0, n1235=0,
                               n1245=0, n1345=0, n2345=0, n12345=0, category = c("CGD_AD","CGD_AR","CGD_XL","CGD_PAE","hOMIM_early"),
                               lty = "blank", fill = c("skyblue", "pink1", "mediumorchid","orange","green"),euler.d=TRUE,scaled = TRUE,margin = 0.05)

# Write the file  
dev.off()  


