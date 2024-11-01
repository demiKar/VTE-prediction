
##### Import Packages

library("readxl")
library("ggplot2")
library("ComplexHeatmap")


##### Initialization

source("/src/pheatmap_fixed.R")

normalized_values <- read_xlsx("/data/Figure2A.xlsx")

##### Plot

figure2A <- pheatmap_fixed(mat = normalized_values,cluster_rows=F,cluster_cols =T, gaps_row = 36,angle_col = "45",  
                    fontsize_row = 5, fontsize_col =18,
                    show_rownames = FALSE,
                    cellwidth=25, cellheight=4,name="Normalized Expr Values",fontsize =20,color=colorRampPalette(c("navy", "white", "red"))(50),
                    heatmap_legend_param = list(title_gp=gpar(fontsize=20,fontface="bold"),
                                                labels_gp = gpar(fontsize = 20)))


png(filename="Figure2A.png",width=9,height=13,res=300, units = "in")

print(figure2A)

dev.off()
