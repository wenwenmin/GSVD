path = "Z:/ZHANGLab_mww/Simulation"
setwd(path)

source('Functions-OR/RSVD_GL0.R')
source('Functions-OR/RSVD_GL1.R')
source('Functions-OR/RSVD_L0.R')
source('Functions-OR/RSVD_ALASSO.R')
source('Functions-OR/OG_auxiliary_functions.R')

library(ggplot2)
library(gridExtra)
library(cowplot)

# logSNR = -2
group.size = 20

# Creat simulation data
n = 100
num.groups = 50 # 50 groups
# group.size = 20 # A group contains 10 members
p = num.groups*group.size

groups = create.group(num.groups,group.size) 
nonzero.groups = c(3,4,13,14,15,33,34,43,44,45) # The groups of nonzero

set.seed(10)
v = rnorm(n, mean = 0, sd = 1)
u = create.u(p, nonzero.groups, group.size, seed=10)

u.true = matrix(u,ncol=1)
v.true = matrix(v,ncol=1)

# logSNR = -2.4
gamma = sqrt(norm(u.true%*%t(v.true),"F")^2/(10^(logSNR)*n*p))

set.seed(20)
X <- u.true%*%t(v.true)+gamma*matrix(rnorm(p*n),ncol=n)

out.svd = svd(X, 1, 1)
kg = length(nonzero.groups)
ku = kg*group.size

out1 = GL0_SVD(X, out.svd, groups, kg)
out2 = GL1_SVD(X, out.svd, groups, kg)
out3 = L0_SVD(X, out.svd, ku)
out4 = AL_SVD(X, out.svd, ku)

# get data.frame for ggplot2
sum.data = data.frame(point = 1:length(u.true), u0= u.true/norm(u.true,"F"), u1 = out1$u, u2 = out2$u, u3 = out3$u, u4 = out4$u)

# plot figure
p0 = ggplot(sum.data, aes(x=point, y = u0)) + geom_line(color = "blue") + xlab(NULL) + ylab(NULL)
p0 = p0 +labs(title = "Original Signal") + theme_bw() #+ scale_y_continuous(breaks=c(-0.1,0,0.1))

p1 = ggplot(sum.data, aes(x=point, y = u1)) + geom_line(color = "blue") + xlab(NULL)+ylab(NULL)
p1 = p1 +labs(title = "GL0-SVD") + theme_bw() + scale_y_continuous(breaks=c(-0.1,0,0.1))

p2 = ggplot(sum.data, aes(x=point, y = u2)) + geom_line(color = "blue") + xlab(NULL)+ylab(NULL) 
p2 = p2 + labs(title = "GL1-SVD") + theme_bw() + scale_y_continuous(breaks=c(-0.1,0,0.1))

p3 = ggplot(sum.data, aes(x=point, y = u3)) + geom_line(color = "blue") + xlab(NULL)+ylab(NULL) 
p3 = p3 + labs(title = "L0-SVD") + theme_bw() + scale_y_continuous(breaks=c(-0.1,0,0.1))

p4 = ggplot(sum.data, aes(x=point, y = u4)) + geom_line(color = "blue") + xlab(NULL)+ylab(NULL)
p4 = p4 + labs(title = "L1-SVD") + theme_bw() + scale_y_continuous(breaks=c(-0.1,0,0.1))

Fig1.OR = plot_grid(p0,p4,p3,p2,p1, ncol=1, nrow =6, align="hv")







