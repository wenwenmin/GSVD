SCAD_SVD = function(X, Initial.uv, ku, niter=100, err=0.0001){
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

    u = SCAD.project(X%*%v0, ku) 
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

# -----------------------------
SCAD.project = function(z, k){ 
  absz = abs(z); 
  sorz = matrix(sort(absz,decreasing = TRUE))
  delta = sorz[k+1]  
  a = 3.7 # a: default choice for SCAD penalty
  
  u = sign(z)*(abs(z)>=delta)*(abs(z)-delta)*(abs(z)<=2*delta)+
    ((a-1)*z-sign(z)*a*delta)/(a-2)*(2*delta<abs(z))*(abs(z)<=a*delta)+z*(abs(z)>a*delta)
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    return(u)} 
}




