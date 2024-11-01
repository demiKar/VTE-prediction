##### Import Packages

library("readxl")
library("data.table")
library("pROC")
library("ggplot2")

##### Functions

get.coords.for.ggplot <- function(roc) {
  df <- coords(roc, "all", transpose = FALSE)
  return(df[rev(seq(nrow(df))),])
}
  
createROCObject <- function(dfSet){
  predictions <- dfSet$CAT
  labels <- ifelse(dfSet$obs=="CAT", 1, 0)
  return(list('roc'=pROC::roc(labels, predictions, ci=TRUE), 
              'auc'=trunc(pROC::auc(labels, predictions) * 10^3)/10^3))
}

##### Initialization

biological_data <- read_xlsx("/data/biological_data_Figure1C.xlsx")
biological_data <- as.data.table(biological_data)

model <- readRDS("/data/Model.rds")


##### ROC curves

test_ROC_plot_colors = c("Chance"="red")
roclist <- list()

###### TOP model

meanROCAUC <- createROCObject(model$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (AUC = %g)', aucAvg)]] <- rocAvg

model <- pROC::coords(rocAvg)
test_ROC_plot_colors[[sprintf('Model (AUC = %g)', aucAvg)]] <- "blue"

###### Peak TF1
peakTF1 <- biological_data[,c("VTE","Peak TF1"), with=F]
peakTF1 <- peakTF1[!is.na(peakTF1$`Peak TF1`)]
obs = ifelse(peakTF1$VTE=="VTE", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, peakTF1$`Peak TF1`)
testAUC <- trunc(pROC::auc(obs, peakTF1$`Peak TF1`) * 10^3)/10^3
peakTF1 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('Peak TF1 (AUC = %g)', testAUC)]] <- "magenta"
  
###### ETP TF1
etpTF1 <- biological_data[,c("VTE","ETP TF1"), with=F]
etpTF1 <- etpTF1[!is.na(etpTF1$`ETP TF1`)]
obs = ifelse(etpTF1$VTE=="VTE", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, etpTF1$`ETP TF1`)
testAUC <- trunc(pROC::auc(obs, etpTF1$`ETP TF1`) * 10^3)/10^3
etpTF1 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('ETP TF1 (AUC = %g)', testAUC)]] <- "orange"
  
###### ETP TF5
etpTF5 <- biological_data[,c("VTE","ETP TF5"), with=F]
etpTF5 <- etpTF5[!is.na(etpTF5$`ETP TF5`)]
obs = ifelse(etpTF5$VTE=="VTE", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, etpTF5$`ETP TF5`)
testAUC <- trunc(pROC::auc(obs, etpTF5$`ETP TF5`) * 10^3)/10^3
etpTF5 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('ETP TF5 (AUC = %g)', testAUC)]] <- "purple"

###### Peak TF5
peakTF5 <- biological_data[,c("VTE","Peak TF5"), with=F]
peakTF5 <- peakTF5[!is.na(peakTF5$`Peak TF5`)]
obs = ifelse(peakTF5$VTE=="VTE", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, peakTF5$`Peak TF5`)
testAUC <- trunc(pROC::auc(obs, peakTF5$`Peak TF5`) * 10^3)/10^3
peakTF5 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('Peak TF5 (AUC = %g)', testAUC)]] <- "darkgreen"
  
###### F1+2
F12 <- biological_data[,c("VTE","F1+2"), with=F]
F12 <- F12[!is.na(F12$`F1+2`)]
obs = ifelse(F12$VTE=="VTE", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, F12$`F1+2`)
testAUC <- trunc(pROC::auc(obs, F12$`F1+2`) * 10^3)/10^3
F12 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('F1+2 (AUC = %g)', testAUC)]] <- "grey"
  
###### Fibrinogen (mg/dl)
Fibrinogen <- biological_data[,c("VTE","Fibrinogen (mg/dl)"), with=F]
Fibrinogen <- Fibrinogen[!is.na(Fibrinogen$`Fibrinogen (mg/dl)`)]
obs = ifelse(Fibrinogen$VTE=="VTE", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, Fibrinogen$`Fibrinogen (mg/dl)`)
testAUC <- trunc(pROC::auc(obs, Fibrinogen$`Fibrinogen (mg/dl)`) * 10^3)/10^3
Fibrinogen <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('Fibrinogen (AUC = %g)', testAUC)]] <- "green"
  
###### ROC plot

rocplot <- ggplot2::ggplot() + 
  ggplot2::geom_line(model, mapping=aes(x=specificity, y=sensitivity, colour="Model (AUC = 0.866)")) +
  ggplot2::geom_line(peakTF1, mapping=aes(x=specificity, y=sensitivity, colour="Peak TF1 (AUC = 0.513)")) +
  ggplot2::geom_line(etpTF1, mapping=aes(x=specificity, y=sensitivity, colour="ETP TF1 (AUC = 0.565)")) +
  ggplot2::geom_line(etpTF5, mapping=aes(x=specificity, y=sensitivity, colour="ETP TF5 (AUC = 0.599)")) +
  ggplot2::geom_line(peakTF5, mapping=aes(x=specificity, y=sensitivity, colour="Peak TF5 (AUC = 0.469)")) +
  ggplot2::geom_line(F12, mapping=aes(x=specificity, y=sensitivity, colour="F1+2 (AUC = 0.692)")) +
  ggplot2::geom_line(Fibrinogen, mapping=aes(x=specificity, y=sensitivity, colour="Fibrinogen (AUC = 0.628)")) +
  ggplot2::scale_x_reverse(lim=c(1, 0)) +
  geom_segment(aes(x=1, xend=0, y=0, yend=1), color='red', linetype='dashed') +
  scale_color_manual(name="Legend", values=test_ROC_plot_colors, breaks = names(test_ROC_plot_colors)) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =9),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = c(0.72,0.24), 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=10, face="bold"),
        axis.text.x = element_text(size=10, face="bold")) 


png(filename="Figure1C.png",width=5,height=4.5,res=300, units = "in")

print(rocplot)

dev.off()
