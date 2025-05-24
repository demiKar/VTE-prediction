# R version 4.3.2 was used to generate the plots.


# ===== Loading libraries =====
library("readxl")
library("data.table")
library("pROC")
library("ggplot2")

# ===== Functions =====

get.coords.for.ggplot <- function(roc) {
  df <- coords(roc, "all", transpose = FALSE)
  return(df[rev(seq(nrow(df))),])
}
  
createROCObject <- function(dfSet){
  predictions <- dfSet$CAT
  labels <- ifelse(dfSet$obs == "CAT", 1, 0)
  
  roc_obj <- pROC::roc(labels, predictions, ci = TRUE)
  auc_value <- round(pROC::auc(labels, predictions), 2) # round to 2 decimals clearly
  
  return(list('roc' = roc_obj, 'auc' = auc_value))
}


# ===== Initialization =====

biological_data <- readxl::read_xlsx("/data/HYPERCAN.xlsx")
biological_data <- as.data.table(biological_data)

model <- readRDS("/data/TOP.rds")


##### ROC curves

test_ROC_plot_colors = c("Chance"="red")
roclist <- list()

###### TOP model

meanROCAUC <- createROCObject(model$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (c-statistic = %g)', aucAvg)]] <- rocAvg

model <- pROC::coords(rocAvg)
test_ROC_plot_colors[[sprintf('Model (c-statistic = %g)', aucAvg)]] <- "blue"

###### Peak TF1
peakTF1 <- biological_data[,c("VTE","Peak TF1"), with=F]
peakTF1 <- peakTF1[!is.na(peakTF1$`Peak TF1`)]
obs = ifelse(peakTF1$VTE=="CAT", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, peakTF1$`Peak TF1`)
testAUC <- round(pROC::auc(obs, peakTF1$`Peak TF1`), 2)
peakTF1 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('Peak TF1 (c-statistic = %g)', testAUC)]] <- "magenta"
  
###### ETP TF1
etpTF1 <- biological_data[,c("VTE","ETP TF1"), with=F]
etpTF1 <- etpTF1[!is.na(etpTF1$`ETP TF1`)]
obs = ifelse(etpTF1$VTE=="CAT", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, etpTF1$`ETP TF1`)
testAUC <- round(pROC::auc(obs, etpTF1$`ETP TF1`), 2)
etpTF1 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('ETP TF1 (c-statistic = %g)', testAUC)]] <- "orange"
  
###### ETP TF5
etpTF5 <- biological_data[,c("VTE","ETP TF5"), with=F]
etpTF5 <- etpTF5[!is.na(etpTF5$`ETP TF5`)]
obs = ifelse(etpTF5$VTE=="CAT", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, etpTF5$`ETP TF5`)
testAUC <- round(pROC::auc(obs, etpTF5$`ETP TF5`), 2)
etpTF5 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('ETP TF5 (AUC = %g)', testAUC)]] <- "purple"

###### Peak TF5
peakTF5 <- biological_data[,c("VTE","Peak TF5"), with=F]
peakTF5 <- peakTF5[!is.na(peakTF5$`Peak TF5`)]
obs = ifelse(peakTF5$VTE=="CAT", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, peakTF5$`Peak TF5`)
testAUC <- round(pROC::auc(obs, peakTF5$`Peak TF5`) , 2)
peakTF5 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('Peak TF5 (c-statistic = %g)', testAUC)]] <- "darkgreen"
  
###### F1+2
F12 <- biological_data[,c("VTE","F1+2"), with=F]
F12 <- F12[!is.na(F12$`F1+2`)]
obs = ifelse(F12$VTE=="CAT", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, F12$`F1+2`)
testAUC <- round(pROC::auc(obs, F12$`F1+2`) , 2)
F12 <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('F1+2 (c-statistic = %g)', testAUC)]] <- "grey"
  
###### Fibrinogen (mg/dl)
Fibrinogen <- biological_data[,c("VTE","Fibrinogen (mg/dl)"), with=F]
Fibrinogen <- Fibrinogen[!is.na(Fibrinogen$`Fibrinogen (mg/dl)`)]
obs = ifelse(Fibrinogen$VTE=="CAT", "CAT", "non.CAT")
  
testROC <- pROC::roc(obs, Fibrinogen$`Fibrinogen (mg/dl)`)
testAUC <- round(pROC::auc(obs, Fibrinogen$`Fibrinogen (mg/dl)`), 2)
Fibrinogen <- get.coords.for.ggplot(testROC)
test_ROC_plot_colors[[sprintf('Fibrinogen (c-statistic = %g)', testAUC)]] <- "green"
  
###### ROC plot

rocplot <- ggplot2::ggplot() + 
  ggplot2::geom_line(model, mapping=aes(x=specificity, y=sensitivity, colour="Model (c-statistic = 0.84)")) +
  ggplot2::geom_line(peakTF1, mapping=aes(x=specificity, y=sensitivity, colour="Peak TF1 (c-statistic = 0.53)")) +
  ggplot2::geom_line(etpTF1, mapping=aes(x=specificity, y=sensitivity, colour="ETP TF1 (c-statistic = 0.51)")) +
  ggplot2::geom_line(etpTF5, mapping=aes(x=specificity, y=sensitivity, colour="ETP TF5 (c-statistic = 0.56)")) +
  ggplot2::geom_line(peakTF5, mapping=aes(x=specificity, y=sensitivity, colour="Peak TF5 (c-statistic = 0.51)")) +
  ggplot2::geom_line(F12, mapping=aes(x=specificity, y=sensitivity, colour="F1+2 (c-statistic = 0.64)")) +
  ggplot2::geom_line(Fibrinogen, mapping=aes(x=specificity, y=sensitivity, colour="Fibrinogen (c-statistic = 0.59)")) +
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


print(rocplot)

