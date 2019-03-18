GL1_SVD = function(X, Initial.uv, groups, kg, niter=100, err=0.0001){
  # ---------------------
  # X: input matrix
  # Initial.uv = out.svd
  # ---------------------
  u0 = Initial.uv$u
  v0 = Initial.uv$v
  
  d0 = -10
  # Iterative algorithm to solve u and v values
  iter.d = NULL
  for (i in 1:niter){
    u = GL1_SVD.project(X%*%v0, groups, kg)
    v = (t(X)%*%u)/norm(t(X)%*%u,"F") 
    
    d = t(u)%*%X%*%v
  
    #print(norm(X - d%*%tt,"F"))
    iter.d = c(iter.d,d)
    
    if(abs(d-d0) < -1){break}
    else {
      u0 = u;v0 = v;d0=d}
  }
  
  return(list(u=u, v=v, iter.d = iter.d))
}

# -----------------------------------
GL1_SVD.project = function(z, group, k){ 
  group.index = unique(group)
  
  z.group = rep(0,length(group.index))
  
  for(i in 1:length(group.index)){
    # wg = 1/length(which(group==group.index[i]))^2
    z.group[i] = norm(z[group==group.index[i]],"2")
  }
  
  u = z
  
  threshold = sort(z.group,decreasing=T)[k+1]
  
  for(i in 1:length(group.index)){
    if(z.group[i]==0|z.group[i]<=threshold) {u[which(group==group.index[i])] = 0}
    else{
      u[which(group==group.index[i])] = (1 - threshold/z.group[i])*u[which(group==group.index[i])]
    }
  }
  
  if(sum(u^2)==0){return(u)}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}




