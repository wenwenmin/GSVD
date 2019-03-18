create.u = function(p, nonzero.groups, group.size, seed){
  set.seed(seed)
  u = sample(c(-1,1), p, replace = T)
  
  id = rep(0,p)
  for(g in nonzero.groups){
    g.members = ((g-1)*group.size+1):(g*group.size)
    id[g.members] = 1
  }
  u[which(id==0)]=0
  return(u)
}

compute.relative.error = function(X,iter.d){
  # iter.d is a value or vector
  obj = norm(X,"F")^2-iter.d^2
  
  Relative.error = obj/norm(X,"F")^2
  
  return(Relative.error) 
}

create.group= function(num.groups,group.size){
  num = num.groups*group.size
  group = rep(0,num)
  for(i in 1:num.groups){
    group[((i-1)*group.size+1):(i*group.size)]=i
  }
  return(group)
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