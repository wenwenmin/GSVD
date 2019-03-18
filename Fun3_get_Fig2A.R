path = "Z:/ZHANGLab_mww/project-4-gSVD/my_code_simulation_v3/Simulation-Fig1"
setwd(path)
ptm <- proc.time()

library(ggplot2)
library(gridExtra)
library(cowplot)

logSNR = -2

source('Fun1_OR.R')
source('Fun2_OGR.R')

plot_grid(p0,p4,p3,p2,p1,p11,p22, ncol=1, nrow = 7, align="hv")
ggsave(paste("Fig.1_logSNR=",logSNR,"_v2.pdf",sep = ""),width = 6, height = 8)

time = proc.time() - ptm; print(time)
save.image("Results.RData")

plot_grid(p0,p4,p3,p2,p1,p11,p22, ncol=1, nrow = 7, align="hv")
ggsave("Fig_2A.pdf",width = 7, height = 11)



