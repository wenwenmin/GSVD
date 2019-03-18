create.overlap.u = function(p, nonzero.groups, overlap.groups,seed){
  set.seed(seed)
  u.sign = sample(c(-1,1), p, replace = T)
  u = rep(0,p)
  for(g in nonzero.groups){
    u[overlap.groups[[g]]]=1
  }
  u.Final = u.sign*u
}

create.group= function(num.groups,group.size){
  # groups = create.group(30,20)
  # num.groups  = 30
  # group.size = 20
  num = num.groups*group.size
  group = rep(0,num)
  for(i in 1:num.groups){
    group[((i-1)*group.size+1):(i*group.size)]=i
  }
  return(group)
}

get_overlap_groups = function(num.groups=50,overlap.size=5){
  # groups2 = get_overlap_groups(50,5)
  overlap.group = list()
  for(i in 1:num.groups){
    overlap.group[[i]] = ((i-1)*overlap.size+1):((i-1)*overlap.size+2*overlap.size)
  } 
  return(overlap.group)
}

# -------------------------------------
evaluation = function(u.true, u.est){
  # https://en.wikipedia.org/wiki/Sensitivity_and_specificity
  # u.est = out1$u
  # (number of) positive samples (P)
  # (number of) negative samples (N)
  # (number of) true positive (TP)
  # (number of) true negative (TN)
  # (number of) false positive (FP)
  # (number of) false negative (FN)
  
  u.est1 = sign(abs(u.est))
  u.true1 = sign(abs(u.true))
  
  P = length(which(u.true1==1))
  N = length(which(u.true1==0))
  
  TP = length(which(u.est1[which(u.true1==1)]==1))
  TN = length(which(u.est1[which(u.true1==0)]==0))
  
  FP = length(which(u.est1[which(u.true1==0)]==1))
  FN = length(which(u.est1[which(u.true1==1)]==0)) 
  
  # sensitivity or true positive rate (TPR)
  TPR = TP/P
  
  # specificity (SPC) or true negative rate (TNR)
  TNR = TN/N
  
  # fall-out or false positive rate (FPR)
  FPR = FP/N
  
  # false discovery rate (FDR)
  FDR = FP/(TP+FP)
  
  # accuracy (ACC)
  ACC = (TP+TN)/(TP+FP+FN+TN)
  
  # F1 score is the harmonic mean of precision and sensitivity
  F1 = 2*TP/(2*TP+FP+FN)
  
  return(list(TPR=TPR,TNR=TNR,FPR=FPR,FDR=FDR,ACC=ACC))
}
# -------------------------------------