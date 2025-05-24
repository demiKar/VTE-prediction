data <- read.table("CD200R1_VTE_MRready.txt", head=T)

#run MR
mrobj <- mr_input(bx=data$BETA_EX, bxse=data$SE_EX, by=data$BETA_OUT, byse=data$SE_OUT)
all <- mr_allmethods(mrobj)


#forest plot
univ_plot <- data.frame("Method" = c("Weighted Median", "Egger", "IVW"),
                        "Estimate" = c(0.037, 0.200, 0.160), 
                        "P" = c(".007", ".140", "0.084"), 
                        "SE" = c(0.014, 0.135, 0.084))
univ_plot$mean <- exp(univ_plot$Estimate)
univ_plot$upper <- exp(ci_normal("u", univ_plot$Estimate, univ_plot$SE, .025))
univ_plot$lower <- exp(ci_normal("l", univ_plot$Estimate, univ_plot$SE, .025))

png(filename = "forest_reversedEffects.png", width = 2560, height = 2560, res=400)
forestplot(
  rep(NA, 3),
  univ_plot$mean,
  univ_plot$lower,
  univ_plot$upper,
  size   = .15,
  zero   = 1,
  clip   = c(0.45, 1.30),                 # ← expand axis range
  xticks = c(0.75, 1, 1.25, 1.50, 1.75),
  xlab   = "Odds Ratio",
  mar    = unit(c(1, 1, 1, 1.5), "cm"),   # ← extra 1.5 cm on the right
  txt_gp = fpTxtGp(
    label = gpar(cex = 3),
    ticks = gpar(cex = 3),
    xlab  = gpar(cex = 4)),
  lwd.ci   = 5,   # CI bars
  lwd.zero = 5,   # vertical zero line
  lwd.xaxis = 5,   # x-axis baseline
  col = fpColors(                # ← make everything black
    box   = "black",
    lines = "black",
    zero  = "black",
    axes  = "black")
)
dev.off()