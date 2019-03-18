GL0_SVD = function(X, Initial.uv, groups, kg, niter=100, err=0.0001){
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
    u = GL0_SVD.project(X%*%v0, groups, kg)
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


# The project algorithm of group constrained SVD 
GL0_SVD.project = function(z, groups, k){ 
  groups.index = unique(groups)
  z.groups = rep(0,length(groups.index))
  for(i in 1:length(groups.index)){
    # p_g = length(which(groups==groups.index[i]))
    z.groups[i] = norm(z[groups==groups.index[i]],"2")
  }
  extract.groups = groups.index[order(z.groups,decreasing=T)[1:k]]
  extract.id = groups%in%extract.groups
  u = z
  u[!extract.id] = 0
  
  if(sum(u^2)==0){return(u)}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}




