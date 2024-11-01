
##### Import Packages

library("pROC")
library("data.table")
library(ggplot2)

##### Functions

createROCObject <- function(dfSet){
  predictions <- dfSet$CAT
  labels <- ifelse(dfSet$obs=="CAT", 1, 0)
  return(list('roc'=pROC::roc(labels, predictions, ci=TRUE), 
              'auc'=trunc(pROC::auc(labels, predictions) * 10^3)/10^3))
}

##### Initialization
khorana_predictions <- fread("/data/Khorana_predictions.txt")
model <- readRDS("/data/Model.rds")

##### Model ROC

roclist = list()

for(i in unique(model$pred$Resample)) {
  foldSet <- filter(model$pred, Resample==i)
  tempROCAUC <- createROCObject(foldSet)
  rocFold <- tempROCAUC[['roc']]
  aucFold <- tempROCAUC[['auc']]
  roclist[[sprintf('%s (AUC = %g)', i, aucFold)]] <- rocFold
}

meanROCAUC <- createROCObject(model$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (AUC = %g)', aucAvg)]] <- rocAvg


##### Khorana ROC
khorana_predictions.roc=pROC::roc(khorana_predictions$VTE, khorana_predictions$Khorana_score, ci=TRUE)
khorana_predictions.auc=trunc(pROC::auc(khorana_predictions$VTE, khorana_predictions$Khorana_score) * 10^3)/10^3
roclist[[sprintf('Khorana ROC (AUC = %g)', khorana_predictions.auc)]] <- khorana_predictions.roc

#####  Graph ROC curves

optimal.threshold.train <- pROC::coords(rocAvg, x="best", best.method="youden", transpose = T)
  
rocplot_title <- sprintf(paste0('\nSpec:%g, Sens:%g'),
                         trunc(optimal.threshold.train[['specificity']]*10^3)/10^3,
                         trunc(optimal.threshold.train[['sensitivity']]*10^3)/10^3)
  
  
rocplot <- ggroc(roclist, aes='colour', legacy.axes = F, aes(size=Folds)) +
  scale_size_manual(values = c(0.25, 0.25, 0.25, 0.25, 0.25, 1, 1)) +
  geom_segment(aes(x=1, xend=0, y=0, yend=1), color='red', linetype='dashed') +
  geom_point(data=data.frame(as.list(optimal.threshold.train)),mapping=aes(x=specificity,y=sensitivity), colour='red') +
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.title = element_text(size=11,color="black",face="bold"), 
          axis.title=element_text(size=11,color="black",face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size =9),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.5, "cm"), 
          legend.position = c(0.75,0.25), 
          axis.line = element_line(size = 0.5, colour = "black"),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=10, face="bold"),
          axis.text.x = element_text(size=10, face="bold"))


png(filename=paste0("Figure1B.png"),width=5,height=4,res=300, units = "in")
  
print(rocplot)  

dev.off()
