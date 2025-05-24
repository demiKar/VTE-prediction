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
              'auc'=trunc(pROC::auc(labels, output) * 10^2)/10^2))
}


#Calculate ROC per fold
createROCObject <- function(dfSet){
  
  predictions <- dfSet$CAT
  labels <- ifelse(dfSet$obs=="CAT", 1, 0)
  return(list('roc'=pROC::roc(labels, predictions, ci=TRUE), 
              'auc'=trunc(pROC::auc(labels, predictions) * 10^2)/10^2))
}

# ===== Initialization =====


source("/src/ggroc.R")

Training_name <- "Lung_Gastric"

set.seed(1)

raw <- readxl::read_xlsx("/data/HYPERCAN.xlsx")
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
 
# ===== Set up 5-fold cross-validation =====


folds <- 5
cvIndex <- caret::createFolds(x$VTE, folds, returnTrain = TRUE)

ctrl <- caret::trainControl(
  index = cvIndex,
  method = "cv",
  number = folds,
  classProbs = TRUE,
  savePredictions = TRUE,
  summaryFunction = twoClassSummary
)

# ===== GLM =====

model_glm <- caret::train(
    VTE ~ .,
    data = x,
    method = "glm",
    family = "binomial", 
    trControl = ctrl,
    metric="ROC"
  )
  


roclist = list()
  
for(i in unique(model_glm$pred$Resample)) {

    foldSet <- filter(model_glm$pred, Resample==i)
    tempROCAUC <- createROCObject(foldSet)
    rocFold <- tempROCAUC[['roc']]
    aucFold <- tempROCAUC[['auc']]
    roclist[[sprintf('%s (AUC = %g)', i, aucFold)]] <- rocFold

}
  
#Calculate Mean ROC
  
meanROCAUC <- createROCObject(model_glm$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (AUC = %g)', aucAvg)]] <- rocAvg
  

# ===== GBM =====


model_gbm <- caret::train(
  VTE ~ .,
  data = x,
  method = "gbm",
  trControl = ctrl,
  metric = "ROC",
  verbose = FALSE
)

roclist = list()

for(i in unique(model_gbm$pred$Resample)) {
  
  foldSet <- filter(model_gbm$pred, Resample==i)
  tempROCAUC <- createROCObject(foldSet)
  rocFold <- tempROCAUC[['roc']]
  aucFold <- tempROCAUC[['auc']]
  roclist[[sprintf('%s (AUC = %g)', i, aucFold)]] <- rocFold
  
}

#Calculate Mean ROC

meanROCAUC <- createROCObject(model_gbm$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (AUC = %g)', aucAvg)]] <- rocAvg


# ===== RF =====


model_rf <- caret::train(
  VTE ~ .,
  data = x,
  method = "rf",
  trControl = ctrl,
  metric = "ROC"
)



roclist = list()

for(i in unique(model_rf$pred$Resample)) {
  
  foldSet <- filter(model_rf$pred, Resample==i)
  tempROCAUC <- createROCObject(foldSet)
  rocFold <- tempROCAUC[['roc']]
  aucFold <- tempROCAUC[['auc']]
  roclist[[sprintf('%s (AUC = %g)', i, aucFold)]] <- rocFold
  
}

#Calculate Mean ROC

meanROCAUC <- createROCObject(model_rf$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (AUC = %g)', aucAvg)]] <- rocAvg



# ===== SVM =====


model_svm <- caret::train(
  VTE ~ .,
  data = x,
  method = "svmRadial",
  trControl = ctrl,
  metric = "ROC"
)


roclist = list()

for(i in unique(model_svm$pred$Resample)) {
  
  foldSet <- filter(model_svm$pred, Resample==i)
  tempROCAUC <- createROCObject(foldSet)
  rocFold <- tempROCAUC[['roc']]
  aucFold <- tempROCAUC[['auc']]
  roclist[[sprintf('%s (AUC = %g)', i, aucFold)]] <- rocFold
  
}

#Calculate Mean ROC

meanROCAUC <- createROCObject(model_svm$pred)
rocAvg <- meanROCAUC[['roc']]
aucAvg <- meanROCAUC[['auc']]
roclist[[sprintf('Mean ROC (AUC = %g)', aucAvg)]] <- rocAvg



# ===== CHECK IN TEST SET =====


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

x_test <- as.data.table(x_test)
# ----------------------------
# Predict on test set
# ----------------------------

# Preprocess test set
x_test_processed <- predict(normalization, x_test %>% select(-VTE))
test_labels <- ifelse(x_test$VTE == "CAT", "CAT", "non.CAT")

# GLM
test_pred_glm <- predict(model_glm, x_test_processed, type = "prob")[, "CAT"]
auc_glm_test <- as.numeric(pROC::auc(test_labels, test_pred_glm))
