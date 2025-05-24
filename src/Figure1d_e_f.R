# R version 4.3.2 was used to generate the plots.


# ===== Loading libraries =====

library(casebase)
library(survival)
library(ggplot2)
library("MASS")
library(cmprsk)
library(dplyr)
library(purrr)
library(survminer)
library(data.table)
library(ggsurvfit)
library(dplyr)

# ===== Initialization =====

set.seed(1)

raw <- readxl::read_xlsx("/data/HYPERCAN.xlsx")
raw <- as.data.table(raw)

predictions <- readxl::read_xlsx("//data/TOPmodel_predictions.xlsx")


raw$PatientID <- as.character(raw$PatientID)
predictions$PatientID <- as.character(predictions$PatientID)
raw <- merge(raw, predictions, by = "PatientID")

raw$prediction <- "Predicted CAT"
raw[meanCAT < optimal.threshold]$prediction <- "Predicted non-CAT"

# ===== Lung =====

lung <- raw[Cancer == "Lung"]

lung$event <- 1
lung[prediction == "Predicted non-CAT"]$event <- 0
lung$`Thrombotic event (1=Yes)` <- as.factor(lung$`Thrombotic event (1=Yes)`)

lung$Status <- with(lung, 
                         ifelse(`Thrombotic event (1=Yes)` == 1, 1, 
                                ifelse(`DEATH (0=NO, 1=YES)` == 1, 2, 0)))

lung$Status <- as.factor(lung$Status)

cumulative_incidence <- tidycmprsk::cuminc(Surv(Days, Status) ~ prediction,data = lung)


dt <- as.data.table(cumulative_incidence$tidy)


lung_vte_cat <- dt[time == 180 & strata == "Predicted CAT" & outcome == "1"]

sprintf(
  "Predicted CAT group: %.1f%% (95%% CI: %.1f%%–%.1f%%)",
  lung_vte_cat$estimate * 100,
  lung_vte_cat$conf.low * 100,
  lung_vte_cat$conf.high * 100
)

lung_nonvte_cat <- dt[time == 180 & strata == "Predicted non-CAT" & outcome == "1"]

sprintf(
  "Predicted non-CAT group: %.1f%% (95%% CI: %.1f%%–%.1f%%)",
  lung_nonvte_cat$estimate * 100,
  lung_nonvte_cat$conf.low * 100,
  lung_nonvte_cat$conf.high * 100
)


p1 <- cumulative_incidence %>% 
  ggcuminc() + 
  labs(
    x = "Days",
    y = "Probability of CAT"
  ) + ggtitle(paste0("Lung cancer (n=",nrow(lung),")")) +   scale_color_manual(values = c("red","blue")) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=14,color="black",face="bold"), 
        axis.title=element_text(size=14,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =14),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = c(0.75,0.9), 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.y = element_text(hjust = 1,size=14, face="bold"),
        axis.text.x = element_text(hjust = 1,size=14, face="bold"),
        strip.text = element_text(size=12)) + 
  annotate("text", x = 0, y = 1, hjust = 0,label = paste0("Gray Test P= ",round(cumulative_incidence$cmprsk$Tests[1, 2], 5)))


print(p1)



# ===== Gastric =====

gastric <- raw[Cancer == "Gastric"]

gastric$event <- 1
gastric[prediction == "Predicted non-CAT"]$event <- 0
gastric$`Thrombotic event (1=Yes)` <- as.factor(gastric$`Thrombotic event (1=Yes)`)

gastric$Status <- with(gastric, 
                    ifelse(`Thrombotic event (1=Yes)` == 1, 1, 
                           ifelse(`DEATH (0=NO, 1=YES)` == 1, 2, 0)))

gastric$Status <- as.factor(gastric$Status)

cumulative_incidence <- tidycmprsk::cuminc(Surv(Days, Status) ~ prediction,data = gastric)

dt <- as.data.table(cumulative_incidence$tidy)


dt[time == 180 & strata == "Predicted non-CAT"]
dt[time == 180 & strata == "Predicted CAT"]

gastric_vte_cat <- dt[time == 180 & strata == "Predicted CAT" & outcome == "1"]

sprintf(
  "Predicted CAT group: %.1f%% (95%% CI: %.1f%%–%.1f%%)",
  gastric_vte_cat$estimate * 100,
  gastric_vte_cat$conf.low * 100,
  gastric_vte_cat$conf.high * 100
)

gastric_nonvte_cat <- dt[time == 180 & strata == "Predicted non-CAT" & outcome == "1"]

sprintf(
  "Predicted non-CAT group: %.1f%% (95%% CI: %.1f%%–%.1f%%)",
  gastric_nonvte_cat$estimate * 100,
  gastric_nonvte_cat$conf.low * 100,
  gastric_nonvte_cat$conf.high * 100
)


p1 <- cumulative_incidence %>% 
  ggcuminc() + 
  labs(
    x = "Days",
    y = "Probability of CAT"
  ) + ggtitle(paste0("Gastric cancer (n=",nrow(gastric),")")) +   scale_color_manual(values = c("red","blue")) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=14,color="black",face="bold"), 
        axis.title=element_text(size=14,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =14),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = c(0.75,0.9), 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.y = element_text(hjust = 1,size=14, face="bold"),
        axis.text.x = element_text(hjust = 1,size=14, face="bold"),
        strip.text = element_text(size=12)) + 
  annotate("text", x = 0, y = 1, hjust = 0,label = paste0("Gray Test P= ",round(cumulative_incidence$cmprsk$Tests[1, 2], 2)))


print(p1)


# ===== AVERT =====

test <- readxl::read_xlsx("/data/AVERT.xlsx")

test <- as.data.table(test)
test$VTE <- NULL


predictions <- readxl::read_xlsx("/data/TOPmodel_predictions_AVERT.xlsx")

test <- merge(test, predictions, by = "PatientID")


test$prediction <- "Predicted CAT"
test[pred == "No CAT"]$prediction <- "Predicted non-CAT"

setnames(test, c("obs","meanCAT"), c("VTE","Predicted_Score"))  

test$`Death` <- as.factor(test$`Death` )
test$`Thrombotic_event` <- as.factor(test$`Thrombotic_event`)


test$Status <- with(test, 
                       ifelse(`Thrombotic_event` == 1, 1, 
                              ifelse(`Death` == 1, 2, 0)))

test$Status <- as.factor(test$Status)

cumulative_incidence <- tidycmprsk::cuminc(Surv(Days, Status) ~ prediction,data = test)
dt <- as.data.table(cumulative_incidence$tidy)

dt[time == 180 & strata == "Predicted non-CAT"]
dt[time == 180 & strata == "Predicted CAT"]


avert_vte_cat <- dt[time == 180 & strata == "Predicted CAT" & outcome == "1"]

sprintf(
  "Predicted CAT group: %.1f%% (95%% CI: %.1f%%–%.1f%%)",
  avert_vte_cat$estimate * 100,
  avert_vte_cat$conf.low * 100,
  avert_vte_cat$conf.high * 100
)

avert_nonvte_cat <- dt[time == 180 & strata == "Predicted non-CAT" & outcome == "1"]

sprintf(
  "Predicted non-CAT group: %.1f%% (95%% CI: %.1f%%–%.1f%%)",
  avert_nonvte_cat$estimate * 100,
  avert_nonvte_cat$conf.low * 100,
  avert_nonvte_cat$conf.high * 100
)


p3 <- cumulative_incidence %>% 
  ggcuminc() + 
  labs(
    x = "Days",
    y = "Probability of CAT"
  ) + ggtitle(paste0("Avert study (n=",nrow(test),")")) +   scale_color_manual(values = c("red","blue")) +
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title = element_text(size=14,color="black",face="bold"), 
        axis.title=element_text(size=14,color="black",face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(colour="black", size =14),
        legend.key = element_rect(fill = "transparent", size = 1,colour = "transparent"),
        legend.key.size = unit(0.5, "cm"), 
        legend.position = c(0.75,0.9), 
        axis.line = element_line(size = 0.5, colour = "black"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        axis.text.y = element_text(hjust = 1,size=14, face="bold"),
        axis.text.x = element_text(hjust = 1,size=14, face="bold"),
        strip.text = element_text(size=12)) + 
  annotate("text", x = 0, y = 1, hjust = 0,label = paste0("Gray Test P= ",round(cumulative_incidence$cmprsk$Tests[1, 2], 2)))

print(p3)


