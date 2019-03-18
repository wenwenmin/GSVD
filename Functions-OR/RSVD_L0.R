L0_SVD = function(X, Initial.uv, ku, niter=100, err=0.0001){
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

    u = L0_SVD.project(X%*%v0, ku) 
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
L0_SVD.project = function(z, k){  
  u = abs(z)
  
  u[-order(u,decreasing=T)[1:k]] = 0
  
  if(sum(u^2)==0){return(rep(0,length(u)))}
  else{
    u = u/sqrt(sum(u^2))
    u = sign(z)*u
    return(u)} 
}






