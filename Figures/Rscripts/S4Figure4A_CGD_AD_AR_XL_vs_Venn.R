library(VennDiagram)
library(grid)
library(gridBase)
library(lattice)

par(mfrow=c(3,1))
par(mai=c(1.02,0.82,0.82,0.42))

#grid.newpage()
hOMIM_AR <- draw.quad.venn(area1 = 876,area2 = 569,area3 = 1500,area4 = 182,n12 = 21,n13 = 0,n14 = 0,n23 = 419,n24 = 1,n34 = 0,n123 = 0,n124 = 0,n134 = 0,n234 = 0,n1234 = 0,category = c("CGD_AD","hOMIM_AR","CGD_AR","CGD_XL"),lty = "blank", fill = c("skyblue", "pink1", "mediumorchid","orange"))
#grid.newpage()
REC_AR <- draw.quad.venn(area1 = 876,area2 = 684,area3 = 1500,area4 = 182,n12 = 16,n13 = 0,n14 = 0,n23 = 513,n24 = 48,n34 = 0,n123 = 0,n124 = 0,n134 = 0,n234 = 0,n1234 = 0,category = c("CGD_AD","REC \nrecessive","CGD_AR","CGD_XL"),lty = "blank", fill = c("skyblue", "pink1", "mediumorchid","orange"))
#grid.newpage()
GENIC_AR <- draw.quad.venn(area1 = 876,area2 = 817,area3 = 1500,area4 = 182,n12 = 77,n13 = 0,n14 = 0,n23 = 491,n24 = 43,n34 = 0,n123 = 0,n124 = 0,n134 = 0,n234 = 0,n1234 = 0,category = c("CGD_AD","GENIC \nrecessive","CGD_AR","CGD_XL"),lty = "blank", fill = c("skyblue", "pink1", "mediumorchid","orange"))


library(gridExtra)
grid.arrange(gTree(children=hOMIM_AR), gTree(children=REC_AR),gTree(children=GENIC_AR), ncol=3 )
#grid.arrange(gTree(children=hOMIM_AR), gTree(children=REC_AR),gTree(children=GENIC_AR), ncol=3, top= "The overlapping between combined models (AD,AR,XL) and respective recessive genes in three published studies")
pdf(file="/home/jacobhsu/ISPP_Paper/Figures/S4Figure_a4size.pdf",height =5, width =45, onefile = TRUE, paper ="a4r",pointsize = 8)
grid.arrange(gTree(children=hOMIM_AR), gTree(children=REC_AR),gTree(children=GENIC_AR), ncol=3 )
# Write the file  
dev.off()  


