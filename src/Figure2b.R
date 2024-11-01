
##### Import Packages

library("readxl")
library("ggplot2")
library("dplyr")
library("lime")

##### Initialization

explanation <- read_xlsx("/data/Figure2B.xlsx")

##### Plot

num_cases <- unique(suppressWarnings(as.numeric(explanation$PatientID)))
if (!anyNA(num_cases)) {
  explanation$case <- factor(explanation$PatientID, levels = as.character(sort(num_cases)))
}
explanation$feature_desc <- forcats::fct_reorder(explanation$feature_desc, abs(explanation$feature_weight), .desc = TRUE) %>%
  forcats::fct_rev()  # Reverse the order


figure2B <- ggplot(explanation, aes_(~case, ~feature_desc)) + 
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

png(filename="Figure2B.png",width=4,height=6,res=300, units = "in")

print(figure2B)

dev.off()
