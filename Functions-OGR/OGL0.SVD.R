OGL0.SVD = function(X, Initial.uv, groups, kg, niter=100, err=0.0001){
  # ---------------------
  # X: input matrix
  # Initial.uv = out.svd
  # ---------------------
  u0 = Initial.uv$u
  v0 = Initial.uv$v
  d0 = -10; iter.d = NULL
  
  for (i in 1:niter){
    u = OGL0.project(X%*%v0, groups, kg)
    v = (t(X)%*%u)/norm(t(X)%*%u,"F") 
    d = t(u)%*%X%*%v
    
    iter.d = c(iter.d,d)
    
    if(abs(d-d0) < -1){break}
    else {
      u0 = u;v0 = v;d0=d}
  }
  return(list(u=u, v=v, iter.d = iter.d))
}

OGL0.project = function(u, overlap.group, kg){
  group.num = length(overlap.group)
  group.norm = rep(0,group.num)
  for(i in 1:group.num){
    g.set = overlap.group[[i]]
    w_i = 1/sqrt(length(g.set))
    group.norm[i] = norm(w_i*u[g.set],"2")
  }
  active.groups = order(group.norm, decreasing = TRUE)[1:kg]
  active.members = NULL
  for(g in active.groups){
    active.members = c(active.members, overlap.group[[g]])
  }
  u[-active.members] = 0
  u = u/norm(u,"F")
}