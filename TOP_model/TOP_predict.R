#!/usr/bin/env Rscript

# =========================================================
# TOP MODEL PREDICTION SCRIPT
# =========================================================

# Usage:
#
# Rscript TOP_predict_apply_args.R \
#   --input INPUT.xlsx \
#   --model TOP.rds \
#   --normalization normalization.rds \
#   --output TOP_predictions.csv
#
# =========================================================

# ---------------------------
# Libraries
# ---------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(caret)
  library(data.table)
  library(dplyr)
  library(openxlsx)
})

# ---------------------------
# Arguments
# ---------------------------

option_list <- list(
  
  make_option(
    c("--input"),
    type = "character",
    help = "Input XLSX file"
  ),
  
  make_option(
    c("--model"),
    type = "character",
    help = "TOP.rds model file"
  ),
  
  make_option(
    c("--normalization"),
    type = "character",
    help = "normalization.rds file"
  ),
  
  make_option(
    c("--output"),
    type = "character",
    default = "TOP_predictions.csv",
    help = "Output CSV file [default= %default]"
  )
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

# ---------------------------
# Check arguments
# ---------------------------

if (is.null(opt$input)) {
  stop("Please provide --input")
}

if (is.null(opt$model)) {
  stop("Please provide --model")
}

if (is.null(opt$normalization)) {
  stop("Please provide --normalization")
}

# ---------------------------
# Load objects
# ---------------------------

cat("Loading normalization object...\n")

normalization <- readRDS(opt$normalization)

cat("Loading TOP model...\n")

TOP <- readRDS(opt$model)

model <- TOP

# ---------------------------
# Fixed Youden threshold
# ---------------------------

TOP_threshold <- 0.09547292

# ---------------------------
# Read input
# ---------------------------

cat("Reading input data...\n")

input_data <- read.xlsx(opt$input)
input_data <- as.data.table(input_data)

# ---------------------------
# Required columns
# ---------------------------

required_columns <- c(
  "PatientID",
  "gal-8",
  "LAT",
  "DAB2",
  "CD200R1",
  "CNDP1",
  "CD84",
  "SPON2",
  "CDH6",
  "FR-gamma",
  "SERPINB8",
  "MARCO",
  "BMI",
  "HGB",
  "Age",
  "Sex",
  "History.of.venous.thrombosis"
)

missing_columns <- setdiff(
  required_columns,
  colnames(input_data)
)

if(length(missing_columns) > 0){
  
  stop(
    paste0(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ---------------------------
# Prepare input
# ---------------------------

x_input <- input_data[, required_columns, with = FALSE]

patient_ids <- x_input$PatientID

x_input$PatientID <- NULL

# ---------------------------
# Normalize
# ---------------------------

cat("Applying normalization...\n")

x_input_norm <- predict(
  normalization,
  x_input
)

# ---------------------------
# Predict probabilities
# ---------------------------

cat("Running TOP predictions...\n")

TOP_probability <- predict(
  model,
  x_input_norm,
  type = "prob"
)[, "CAT"]

# ---------------------------
# High / Low classification
# ---------------------------

TOP_risk_group <- ifelse(
  TOP_probability >= TOP_threshold,
  "High",
  "Low"
)

# ---------------------------
# Final output
# ---------------------------

results <- data.frame(
  PatientID = patient_ids,
  TOP_probability = TOP_probability,
  Youden_threshold = TOP_threshold,
  TOP_risk_group = TOP_risk_group
)

# ---------------------------
# Save output
# ---------------------------

write.csv(
  results,
  opt$output,
  row.names = FALSE
)

cat("\nPrediction completed successfully.\n")
cat(paste0("Results saved to: ", opt$output, "\n"))

# =========================================================
# END
# =========================================================
