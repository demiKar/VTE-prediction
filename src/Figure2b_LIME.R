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
library(lime)
library(modeltools)

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
              'auc'=trunc(pROC::auc(labels, output) * 10^3)/10^3))
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
  labels <- ifelse(dfSet$obs=="CAT", 1, 0)
  return(list('roc'=pROC::roc(labels, predictions, ci=TRUE), 
              'auc'=trunc(pROC::auc(labels, predictions) * 10^3)/10^3))
}

# ===== Initialization =====

source("/src/ggroc.R")

Training_name <- "Lung_Gastric"

set.seed(1)

raw <- readxl::read_xlsx("/data/HYPERCAN.xlsx")
raw <- as.data.table(raw)

protein_features <- c("PatientID","VTE","gal-8","LAT","DAB2","CD200R1","CNDP1","CD84","SPON2","CDH6","FR-gamma","SERPINB8","MARCO","Age","Sex","HGB","BMI","History.of.venous.thrombosis")


train_patientid_list <- raw$PatientID
train_patientid_list <- as.character(train_patientid_list)

train <- raw[, c(protein_features), with=F]
  
train$VTE <- as.factor(train$VTE)

# ===== Prepocessing =====

  
normalization <- preProcess(train[,3:ncol(train)], method = c("medianImpute", "center", "scale"))

train.normalized <- predict(normalization, train[,3:ncol(train)])
  
train.normalized$VTE <-as.data.frame(train)$VTE
  
x <- train.normalized
  
train_patients <- train$PatientID

  
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
  
  
  
  explainer <- lime(x[,1:16], model, bin_continuous = FALSE, quantile_bins = FALSE)
  explanation <- lime::explain(x[,1:16], explainer, n_labels = 1, n_features = 16)
  
  explanation[explanation$feature_desc == "History.of.venous.thrombosis",]$feature_desc <- "History of VTE"
  explanation[explanation$feature_desc == "gal-8",]$feature_desc <- "LGALS8"
  explanation[explanation$feature_desc == "HGB",]$feature_desc <- "Hb"
  
  explanation[explanation$label == "CAT",]$label <- "CAT"
  explanation[explanation$label == "non.CAT",]$label <- "non-CAT"
  
  explanation <- as.data.table(explanation)
  
  num_cases <- unique(suppressWarnings(as.numeric(explanation$case)))
  if (!anyNA(num_cases)) {
    explanation$case <- factor(explanation$case, levels = as.character(sort(num_cases)))
  }
  explanation$feature_desc <- forcats::fct_reorder(explanation$feature_desc, abs(explanation$feature_weight), .desc = TRUE) %>%
    forcats::fct_rev()  # Reverse the order
  
  # ===== LimePlot =====
  
    ggplot(explanation, aes_(~case, ~feature_desc)) + 
    geom_tile(aes_(fill = ~feature_weight)) + 
    facet_wrap(~label) +
    scale_x_discrete("Samples", expand = c(0, 0)) + 
    scale_y_discrete("", expand = c(0, 0)) + 
    scale_fill_gradient2("Feature\nWeight", 
                         low = "steelblue", mid = "grey", high = "firebrick") + 
    theme(panel.border = element_rect(fill = NA, colour = "white", size = 1), 
          panel.grid = element_blank(), axis.text.y = element_text(colour="black"),
          legend.position = "right", 
          axis.text.x = element_blank())
  
