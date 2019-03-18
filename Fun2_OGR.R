path = "Z:/ZHANGLab_mww/project-4-gSVD/my_code_simulation_v3/Simulation-Fig1"
setwd(path)

source('Functions-OGR/OGL1.SVD.R')
source('Functions-OGR/OGL0.SVD.R')
source('Functions-OGR/OGR_auxiliary_functions.R')

library(ggplot2)
library(gridExtra)
library(cowplot)

# logSNR = -2
# Creat simulation data
n = 100
num.groups = 49 #20 groups
overlap.size = 20 #i group and (i+1) group contain 10 overlapping members

overlap.groups = get_overlap_groups(num.groups,overlap.size)
p = num.groups*overlap.size + overlap.size; 

nonzero.groups = c(3, 13,14, 33, 43,44) # The groups of nonzero

set.seed(100)
v = rnorm(n, mean = 0, sd = 1)
u = create.overlap.u(p, nonzero.groups, overlap.groups,seed=10)

u.true = matrix(u,ncol=1)
v.true = matrix(v,ncol=1)

gamma = sqrt(norm(u.true%*%t(v.true),"F")^2/(10^(logSNR)*n*p))

set.seed(20)
X <- u.true%*%t(v.true)+gamma*matrix(rnorm(p*n),ncol=n)
out.svd = svd(X, 1, 1)

kg = length(nonzero.groups)
ku = length(which(u.true!=0))

out1 = OGL1.SVD(X, out.svd, overlap.groups, kg)
out2 = OGL0.SVD(X, out.svd, overlap.groups, kg)

# plot(out2$iter.d)
# plot(out1$iter.d)

sum.data = data.frame(point = 1:length(u.true), u0= u.true/norm(u.true,"F"), u1 = out1$u, u2 = -out2$u)

p11 = ggplot(sum.data, aes(x=point, y = u1)) + geom_line(color = "brown") + xlab(NULL)+ylab(NULL)
p11 = p11 +labs(title = "OGL1-SVD") + theme_bw() + scale_y_continuous(breaks=c(-0.1,0,0.1))

p22 = ggplot(sum.data, aes(x=point, y = u2)) + geom_line(color = "brown") + xlab(NULL)+ylab(NULL)
p22 = p22 +labs(title = "OGL0-SVD") + theme_bw() + scale_y_continuous(breaks=c(-0.1,0,0.1))


