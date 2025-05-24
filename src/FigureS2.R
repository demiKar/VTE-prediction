# R version 4.3.2 was used to generate the plots.


# ===== Loading libraries =====

library(data.table)
library(caret)
library(predtools)
library(magrittr)
library(dplyr)
library(tidymodels)
library(probably)
library(ggplot2)
library(Hmisc)

# ===== Initialization =====

predictions <- readRDS("/data/TOP.rds")
predictions <- predictions$pred
setDT(predictions)

data_unique <- predictions %>%
  group_by(rowIndex) %>%  
  summarise(obs = first(obs),  
            predicted_probs = mean(CAT))  

predicted_probs <- data_unique$predicted_probs
actual_outcomes <- factor(as.numeric(data_unique$obs == "CAT"), levels = c(1, 0))

data2 <- data.frame(obs = actual_outcomes, predicted_probs = predicted_probs)
data2$group <- "TOP model"

data <- data2

# Convert obs to numeric
setDT(data)
data$obs <- as.character(data$obs)
data$obs <- as.numeric(data$obs)

data$predicted_probs <- as.character(data$predicted_probs)
data$predicted_probs <- as.numeric(data$predicted_probs)

compute_ci <- function(data, num_breaks = 5) {
  data$bin <- cut(data$predicted_probs, breaks = seq(0, 1, length.out = num_breaks + 1), include.lowest = TRUE)
  summary_data <- data %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(predicted_probs),
      obs_rate = mean(as.numeric(obs)), 
      n = n()
    ) %>%
    mutate(
      lower_ci = Hmisc::binconf(n * obs_rate, n, alpha = 0.05)[,2],
      upper_ci = Hmisc::binconf(n * obs_rate, n, alpha = 0.05)[,3]
    )
  return(summary_data)
}


calibration_data <- compute_ci(data, num_breaks = 4)  # Try 7 instead of 5
calibration_data <- calibration_data %>%
  arrange(mean_pred)  # Order bins properly

# Plot Calibration Curve with 95% Confidence Intervals
ggplot(calibration_data, aes(x = mean_pred, y = obs_rate)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.02, color = "black") +
  geom_line(color = "black", size = 0.5) +
  geom_point(size = 2, color = "black") +
  geom_text(aes(label = n), vjust = -1, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(x = "Mean Predicted VTE Probability", y = "Fraction of VTE Events")

