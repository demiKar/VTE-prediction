# R version 4.3.2 was used to generate the plots.


# ===== Loading libraries =====

library("data.table")
library(ggplot2)
library("ggrepel")
library("cowplot")
library(grid)
library(gridExtra)
library(ggpubr)
library(scater)
library(ggplotify)
library(RColorBrewer)
library(viridis)
library("rjson")
library(MASS)

# ===== Initialization =====

source("/src/pheatmap_fixed.R")

raw <- readxl::read_xlsx("/data/HYPERCAN.xlsx")
raw <- as.data.table(raw)

protein_features <- c("PatientID","VTE","gal-8","LAT","DAB2","CD200R1","CNDP1","CD84","SPON2","CDH6","FR-gamma","SERPINB8","MARCO")

quantile_normalized_values <-raw

# ===== Figure2a =====

quantile_normalized_values <- quantile_normalized_values[,c(protein_features), with=F]
setnames(quantile_normalized_values,"gal-8","LGALS8")
quantile_normalized_values[VTE == "non.CAT"]$VTE <- "non-CAT"
quantile_normalized_values[VTE == "CAT"]$VTE <- "CAT"

samples <- quantile_normalized_values
annRow <- setDT(samples[,c("PatientID","VTE"), with=F])
annRow <- annRow[order(annRow$VTE)]
annRow <- setDF(annRow)
rownames(annRow) <- annRow$PatientID
annRow$`PatientID` <- as.factor(annRow$`PatientID`)
rownames(annRow) <- factor(rownames(annRow), levels = rownames(annRow))
annRow$PatientID <- NULL

setDF(quantile_normalized_values)
rownames(quantile_normalized_values) <- quantile_normalized_values$PatientID
quantile_normalized_values$PatientID <- NULL
quantile_normalized_values$VTE <- NULL
quantile_normalized_values <- quantile_normalized_values[rownames(annRow),]

getPalette = colorRampPalette(brewer.pal(n = 8, name = "Spectral"))
mycolors <- c("yellowgreen","black")
names(mycolors) <- unique(annRow$`VTE`)
mycolors <- list(`VTE` = mycolors)


scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}
quantile_normalized_values <- scale_mat(quantile_normalized_values, "column")

colors <-  c("#000004FF","#280B54FF","#65156EFF","#9F2A63FF","#D44842FF","#F57D15FF","#FAC127FF","#FCFFA4FF")


figureC <- pheatmap_fixed(mat = quantile_normalized_values,cluster_rows=T,cluster_cols =T, gaps_row = 35,
                          angle_col = "45",  
                    right_annotation = rowAnnotation(df = annRow, 
                                                     col = mycolors, 
                                                     annotation_name_gp = gpar(fontsize = 20, fontface = "bold"),  # Font for annotation labels
                                                     annotation_legend_param = list(
                                                       title_gp = gpar(fontsize = 20, fontface = "bold"),  # Font for legend title
                                                       labels_gp = gpar(fontsize = 20)  # Font for legend labels
                                                     )),
                    fontsize_row = 5, fontsize_col =18,
                    show_rownames = FALSE,
                    cellwidth=25, cellheight=4,name="Normalized Expr Values",fontsize =20,color=colorRampPalette(c("navy", "white", "red"))(50),
                    heatmap_legend_param = list(title_gp=gpar(fontsize=20,fontface="bold"),
                                                labels_gp = gpar(fontsize = 20)))



print(figureC)



# ===== Figure S1 =====
quantile_normalized_values <-raw

quantile_normalized_values <- quantile_normalized_values[,c(protein_features), with=F]
setnames(quantile_normalized_values,"gal-8","LGALS8")
quantile_normalized_values[VTE == "non.CAT"]$VTE <- "non-CAT"
quantile_normalized_values[VTE == "CAT"]$VTE <- "CAT"

samples <- quantile_normalized_values
annRow <- setDT(samples[,c("PatientID","VTE"), with=F])
annRow <- annRow[order(annRow$VTE)]
annRow <- setDF(annRow)
rownames(annRow) <- annRow$PatientID
annRow$`PatientID` <- as.factor(annRow$`PatientID`)
rownames(annRow) <- factor(rownames(annRow), levels = rownames(annRow))
annRow$PatientID <- NULL

setDF(quantile_normalized_values)
rownames(quantile_normalized_values) <- quantile_normalized_values$PatientID
quantile_normalized_values$PatientID <- NULL
quantile_normalized_values$VTE <- NULL
quantile_normalized_values <- quantile_normalized_values[rownames(annRow),]


annRow$sample <- rownames(annRow)
quantile_normalized_values <- as.data.frame(quantile_normalized_values)
all.pct <- quantile_normalized_values


PCA1 <- prcomp(all.pct, center = TRUE)
summary(PCA1)
pca_data_perc=round(100*(PCA1$sdev^2/sum(PCA1$sdev^2)),1)
df_pca_data=data.frame(PC1 = PCA1$x[,1], PC2 = PCA1$x[,2], sample = rownames(all.pct))

df_pca_data <- merge(df_pca_data, annRow, by.x="sample", by.y="sample")
df_pca_data$status <- df_pca_data[,"VTE"]


a<-1
b<-2
p1 <- ggplot(df_pca_data, aes(PC1, PC2, color = status)) +
  geom_point(size = 0.7, alpha = 0.8) +  
  stat_ellipse(aes(color = status), type = "norm", level = 0.7, linetype = "dashed") +
  scale_fill_manual(values = c( "limegreen","grey27")) +
  scale_color_manual(values = c("limegreen","grey27")) +
  labs(
    x = paste0("PC", a, " (", pca_data_perc[a], "%)"),
    y = paste0("PC", b, " (", pca_data_perc[b], "%)")
  ) +
  theme(
    axis.line = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_line(size = 0.1, colour = "grey"),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size = 8, colour = "black", face = "bold"),
    axis.text.x = element_text(size = 8, hjust = 1, color = "black", face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 8),
    legend.key.size = unit(0.3, "cm"),
    plot.title = element_text(size = 10, hjust = 0.5, lineheight = .7, face = "bold"),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
  )



print(p1)

