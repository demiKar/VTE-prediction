# R version 4.3.2 was used to generate the plots.


# ===== Loading libraries =====
library(data.table)


# ===== Initialization =====

# ===== HYPERCAN =====

training_cohort <- readxl::read_xlsx("/data/HYPERCAN.xlsx")
training_cohort <- as.data.table(training_cohort)
training_cohort[VTE=="non.CAT"]$VTE <- 0
training_cohort[VTE=="CAT"]$VTE <- 1
training_cohort$VTE <- as.numeric(as.character(training_cohort$VTE))  

# Define your biomarkers and traits
biomarkers <- c("CD84","SERPINB8","SPON2","CDH6","MARCO","DAB2","FRgamma","gal8","LAT","CNDP1","CD200R1")
traits <-c( "CRP", "BMI","HGB","PLT","WBC", "History.of.venous.thrombosis", "Sex", "Age")

setnames(training_cohort, c("FR-gamma","gal-8"),c("FRgamma","gal8"))

# imputation with median
for (col in traits) {
  if (!is.numeric(training_cohort[[col]])) {
    training_cohort[[col]] <- as.numeric(as.character(training_cohort[[col]]))
  }
  
  if (any(is.na(training_cohort[[col]]))) {
    med <- median(training_cohort[[col]], na.rm = TRUE)
    training_cohort[[col]][is.na(training_cohort[[col]])] <- med
  }
}

# Collect results here
all_results <- list()

# Loop over biomarkers
for (bm in biomarkers) {
  # --- Unadjusted model ---
  model_unadj <- glm(as.formula(paste("VTE ~", bm)), data = training_cohort, family = "binomial")
  conf_unadj <- suppressMessages(confint(model_unadj))
  
  OR_unadj <- exp(coef(model_unadj)[bm])
  CI_unadj <- exp(conf_unadj[bm, ])
  p_unadj <- summary(model_unadj)$coefficients[bm, 4]
  
  unadj_combined <- sprintf("%.2f [%.2f–%.2f]; p = %.10f", OR_unadj, CI_unadj[1], CI_unadj[2], p_unadj)
  
  biomarker_result <- data.frame(Biomarker = bm, Unadjusted = unadj_combined, stringsAsFactors = FALSE)
  
  OR_ref <- OR_unadj
  log_OR_unadj <- log(OR_ref)
  
  # --- Adjusted (1 trait at a time) ---
  for (tr in traits) {
    form <- as.formula(paste("VTE ~", bm, "+", tr))
    model2 <- glm(form, data = training_cohort, family = "binomial")
    conf <- suppressMessages(confint(model2))
    
    OR <- exp(coef(model2)[bm])
    CI <- exp(conf[bm, ])
    pval <- summary(model2)$coefficients[bm, 4]
    
    combined <- sprintf("%.2f [%.2f–%.2f]; p = %.10f", OR, CI[1], CI[2], pval)
    
    log_OR_adj <- log(OR)
    delta_log_OR <- log_OR_adj - log_OR_unadj
    
    if (abs(delta_log_OR) < 1e-4) {
      percent_change <- 0
    } else {
      percent_change <- round(100 * delta_log_OR / abs(log_OR_unadj), 1)
    }
    
    abs_diff <- round(OR - OR_ref, 2)
    
    if ((OR_ref < 1 && OR < 1) || (OR_ref > 1 && OR > 1)) {
      direction <- "Same direction"
    } else {
      direction <- "Flipped"
    }
    
    biomarker_result[[paste0(tr, "_adj")]] <- combined
    biomarker_result[[paste0(tr, "_%change")]] <- percent_change
    biomarker_result[[paste0(tr, "_abs_diff")]] <- abs_diff
    biomarker_result[[paste0(tr, "_direction")]] <- direction
  }
  
  all_results[[bm]] <- biomarker_result
}

final_results <- do.call(rbind, all_results)

for (tr in traits) {
  colnames(final_results)[colnames(final_results) == paste0(tr, "_adj")] <- paste0(tr, " (adj)")
  colnames(final_results)[colnames(final_results) == paste0(tr, "_%change")] <- paste0(tr, " (% logOR change)")
  colnames(final_results)[colnames(final_results) == paste0(tr, "_abs_diff")] <- paste0(tr, " (abs OR change)")
  colnames(final_results)[colnames(final_results) == paste0(tr, "_direction")] <- paste0(tr, " (direction)")
}


# ===== AVERT =====

validation_cohort <- readxl::read_xlsx("/data/AVERT.xlsx")
validation_cohort <- as.data.table(validation_cohort)
validation_cohort[VTE=="non.CAT"]$VTE <- 0
validation_cohort[VTE=="CAT"]$VTE <- 1
validation_cohort$VTE <- as.numeric(as.character(validation_cohort$VTE))  

setnames(validation_cohort, c("FR-gamma","gal-8"),c("FRgamma","gal8"))

### imputation with median
for (col in traits) {
  if (!is.numeric(validation_cohort[[col]])) {
    validation_cohort[[col]] <- as.numeric(as.character(validation_cohort[[col]]))
  }
  
  if (any(is.na(validation_cohort[[col]]))) {
    med <- median(validation_cohort[[col]], na.rm = TRUE)
    validation_cohort[[col]][is.na(validation_cohort[[col]])] <- med
  }
}

validation_cohort$VTE <- as.numeric(as.character(validation_cohort$VTE))  


all_unadj_results <- list()

# Loop over biomarkers
for (bm in biomarkers) {
  # Fit unadjusted logistic regression model
  model_unadj <- glm(as.formula(paste("VTE ~", bm)), data = validation_cohort, family = "binomial")
  conf_unadj <- suppressMessages(confint(model_unadj))
  
  # Extract OR, CI, and p-value
  OR_unadj <- exp(coef(model_unadj)[bm])
  CI_unadj <- exp(conf_unadj[bm, ])
  p_unadj <- summary(model_unadj)$coefficients[bm, 4]
  
  result_combined <- sprintf("%.2f [%.2f–%.2f]; p = %.10f", OR_unadj, CI_unadj[1], CI_unadj[2], p_unadj)
  
  row <- data.frame(
    Biomarker = bm,
    OR_CI_p = result_combined,
    OR = round(OR_unadj, 2),
    CI_low = round(CI_unadj[1], 2),
    CI_high = round(CI_unadj[2], 2),
    p_value = round(p_unadj, 4),
    stringsAsFactors = FALSE
  )
  
  all_unadj_results[[bm]] <- row
}

unadjusted_results <- do.call(rbind, all_unadj_results)


