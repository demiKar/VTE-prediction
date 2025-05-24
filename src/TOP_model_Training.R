# R version 4.3.2 was used to generate the plots.


# ===== Loading libraries =====

library("tidyverse")
library("caret")
library("pROC")
library("cowplot")
library("psych")
library("magrittr")
library("data.table")
library("yardstick")
library("openxlsx")
library(ggplot2)
library(dplyr)
  

# ===== Functions =====

khoranaROCCurves <- function(df=raw, 
                             patients=train_patientid_list, 
                             orig_train=train, label="train") {
  
  
  
  if(label == "train"){
    
    output <- df$khorana    
    
  }else{
    
    df <- as.data.table(df)
    output <- df$Khorana.Score
    
  }
  
  output <- ifelse(output >= 2, 1, 0)
  
  labels <- ifelse(df$VTE=="CAT", 1, 0)
  
  return(list('roc'=pROC::roc(labels, output, ci=TRUE), 
              'auc'=round(pROC::auc(labels, output),2)))
}


khoranaCalculation <- function(df=khoranaDF) {
  
  
    df <- as.data.table(df)
    df$khorana <- 0
    df[Cancer == 'Lung', khorana:=khorana+1]
    df[Cancer == 'Gastric', khorana:=khorana+2]
    df[BMI>=35,khorana:=khorana+1]
    df[`HGB`<10,khorana:=khorana+1]
    df[`PLT`>=350,khorana:=khorana+1]
    df[`WBC`>=11,khorana:=khorana+1]
    
  return(df)
}

#Calculate ROC per fold
createROCObject <- function(dfSet){
  predictions <- dfSet$CAT
  labels <- ifelse(dfSet$obs == "CAT", 1, 0)
  
  roc_obj <- pROC::roc(labels, predictions, ci = TRUE)
  auc_value <- round(pROC::auc(labels, predictions), 2) # round to 2 decimals clearly
  
  return(list('roc' = roc_obj, 'auc' = auc_value))
}

# ===== Initialization =====

source("/src/ggroc.R")

Training_name <- "Lung_Gastric"

set.seed(1)

raw <- readxl::read_xlsx("//data/archive/dimitra/Thrombosis/VTE/Publication/Reviews/src/code_online/data/HYPERCAN.xlsx")
raw <- as.data.table(raw)


protein_features <- c("PatientID","VTE","gal-8","LAT","DAB2","CD200R1","CNDP1","CD84","SPON2","CDH6","FR-gamma","SERPINB8","MARCO","BMI","HGB","Age","Sex","History.of.venous.thrombosis")


train_patientid_list <- raw$PatientID

train <- raw[, c(protein_features), with=F]
  
train$VTE <- as.factor(train$VTE)

# ===== Prepocessing =====

  
normalization <- preProcess(train[,3:ncol(train)], method = c("medianImpute", "center", "scale"))
  
train.normalized <- predict(normalization, train[,3:ncol(train)])
  
train.normalized$VTE <-as.data.frame(train)$VTE
  
x <- train.normalized
  
train_patients <- train$PatientID
 
# ===== Training =====


  # Model Training -----------------------------------------------------------
  folds <- 5
  cvIndex <- caret::createFolds(x$VTE, folds, returnTrain=T)
  
  ctrl <- caret::trainControl(
    index=cvIndex,
    method='cv',
    number=folds,
    classProbs=T,
    savePredictions=T,
    summaryFunction = twoClassSummary)
  
  model <- caret::train(
    VTE ~ .,
    data = x,
    method="naive_bayes",
    trControl = ctrl,
    metric="ROC"
  )
  
# RocList
  
roclist = list()
  
for(i in unique(model$pred$Resample)) {

    foldSet <- filter(model$pred, Resample==i)
    tempROCAUC <- createROCObject(foldSet)
    rocFold <- tempROCAUC[['roc']]
    aucFold <- tempROCAUC[['auc']]
    roclist[[sprintf('%s (c-statistic = %g)', i, aucFold)]] <- rocFold

}
  
#Calculate Mean ROC
  
meanROCAUC <- createROCObject(model$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (c-statistic = %g)', aucAvg)]] <- rocAvg
  

# ===== Compare with khorana =====

khoranaROCAUC_train <- khoranaROCCurves(raw,train_patientid_list,train, label = "train")
khoranaROC_train <- khoranaROCAUC_train[['roc']]
khoranaAUC_train <- khoranaROCAUC_train[['auc']]
roclist[[sprintf('Khorana ROC (c-statistic = %g)', khoranaAUC_train)]] <- khoranaROC_train
  
# Graph ROC curves
optimal.threshold.train <- pROC::coords(rocAvg, x="best", best.method="youden", transpose = T)

rocplot_title <- sprintf(paste0('\nSpec:%g, Sens:%g'), 
                           round(optimal.threshold.train[['specificity']],2),
                         round(optimal.threshold.train[['sensitivity']],2))
  
  
rocplot <- ggroc(roclist, aes='colour', legacy.axes = F, aes(size=Folds)) +
    scale_size_manual(values = c(0.25, 0.25, 0.25, 0.25, 0.25, 1, 1)) +
    geom_segment(aes(x=1, xend=0, y=0, yend=1), color='red', linetype='dashed') +
    geom_point(data=data.frame(as.list(optimal.threshold.train)), 
               mapping=aes(x=specificity,y=sensitivity), colour='red') +
    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          plot.title = element_text(size=11,color="black",face="bold"), 
          axis.title=element_text(size=11,color="black",face="bold"),
          legend.title = element_blank(),
          legend.text = element_text(colour="black", size =9),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(0.5, "cm"), 
          legend.position = c(0.8,0.2), 
          axis.line = element_line(size = 0.5, colour = "black"),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),axis.text.y = element_text(size=10, face="bold"),
          axis.text.x = element_text(size=10, face="bold")) #+ 
    
pred.dt <- as.data.table(model$pred)
pred.dt <- pred.dt[,meanCAT := mean(CAT), by="rowIndex"][,meannon.CAT := mean(non.CAT), by="rowIndex"]
  
pred.dt <- unique(pred.dt[,c("pred","obs","meanCAT","meannon.CAT","rowIndex")])
pred.dt <- pred.dt[order(rowIndex)]
  
pred.dt$pred.new <- ifelse(pred.dt$meanCAT>as.numeric(optimal.threshold.train[['threshold']]),"CAT","non.CAT")
pred.dt$pred <- NULL
pred.dt <- unique(pred.dt)
pred.dt <- pred.dt[order(rowIndex)]
  
df_train <- data.frame(
  obs = train$VTE,
  pred = pred.dt$pred.new)
  
df_train$obs <- as.factor(df_train$obs)
levels(df_train$obs) <- c( 'CAT', 'No CAT')
df_train$pred <- as.factor(df_train$pred)
levels(df_train$pred) <- c( 'CAT', 'No CAT')
  
  
df_train$PatientID <- train_patients
df_train$meanCAT <- pred.dt$meanCAT
df_train$optimal.threshold <- as.numeric(optimal.threshold.train[['threshold']])



print(rocplot)  
  


# For mean ROC AUC
ci_auc_avg <- pROC::ci.auc(rocAvg)
print(ci_auc_avg)


# ===== Test =====


test.in <- read.xlsx("/data/AVERT.xlsx")
test.in <- as.data.table(test.in)
test.in$VTE <- as.character(test.in$VTE)
test.in[VTE == 1]$VTE <- "CAT"
test.in[VTE == 0]$VTE <- "non.CAT"

test_patientid_list <- test.in$PatientID
test <- test.in[, c(protein_features), with=F]
test$VTE <- as.factor(test$VTE)

x_test <- test
test_patients <- x_test$PatientID
x_test$PatientID <- NULL

# Preprocess Testing and Predict
test_prediction.norm = predict(normalization, x_test %>% select(-VTE))
  
test_prediction = predict(model,test_prediction.norm, type='prob', )[,"CAT"]
test_obs = ifelse(x_test$VTE=="CAT", "CAT", "non.CAT")
  
testROC <- pROC::roc(test_obs, test_prediction, ci=TRUE)
testAUC <- round(pROC::auc(test_obs, test_prediction),2)


  
  
