
##### Import Packages

library("readxl")
library("data.table")
library("ggplot2")
library("dplyr")
library("ggsurvfit")

##### Initialization

cumulative_lung <- readRDS("/data/Figure1D.RDS")
cumulative_gastric <- readRDS("/data/Figure1E.RDS")
cumulative_avert <- readRDS("/data/Figure1F.RDS")


##### Plot cumulative 

figure1d <- cumulative_lung %>% 
  ggcuminc() + 
  labs(x = "Days",y = "Probability of VTE") + ggtitle(paste0("Lung cancer (n=71)")) +   scale_color_manual(values = c("red","blue")) +
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
  annotate("text", x = 0, y = 1, hjust = 0,label = paste0("Gray Test P= ",round(cumulative_lung$cmprsk$Tests[1, 2], 5)))

png("Figure1D.png",units="in", width=7, height=4, res=500)

print(figure1d)

dev.off()


figure1e <- cumulative_gastric %>% 
  ggcuminc() + 
  labs(x = "Days",y = "Probability of VTE") + ggtitle(paste0("Gastric cancer (n=92)")) +   scale_color_manual(values = c("red","blue")) +
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
  annotate("text", x = 0, y = 1, hjust = 0,label = paste0("Gray Test P= ",round(cumulative_gastric$cmprsk$Tests[1, 2], 4)))

png("Figure1E.png",units="in", width=7, height=4, res=500)

print(figure1e)

dev.off()


figure1f <- cumulative_avert %>% 
  ggcuminc() + 
  labs(x = "Days",y = "Probability of VTE") + ggtitle(paste0("Avert study (n=72)")) +   scale_color_manual(values = c("red","blue")) +
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
  annotate("text", x = 0, y = 1, hjust = 0,label = paste0("Gray Test P= ",round(cumulative_avert$cmprsk$Tests[1, 2], 2)))

png("Figure1F.png",units="in", width=7, height=4, res=500)

print(figure1f)

dev.off()
