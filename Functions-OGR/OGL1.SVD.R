OGL1.SVD = function(X, Initial.uv, groups, kg, niter=100, err=0.0001){
  # ---------------------
  # X: input matrix
  # Initial.uv = out.svd
  # ---------------------
  u0 = Initial.uv$u
  v0 = Initial.uv$v
  d0 = -10; iter.d = NULL
  
  #node.degree = get.node.degree(groups)
  for (i in 1:niter){
    u = ADMM.project.problem1(X%*%v0, groups, kg)
    v = (t(X)%*%u)/norm(t(X)%*%u,"F") 
    d = t(u)%*%X%*%v
    
    iter.d = c(iter.d,d)
    
    if(abs(d-d0) < -1){break}
    else {
      u0 = u;v0 = v;d0=d}
  }
  return(list(u=u, v=v, iter.d = iter.d))
}

ADMM.project.problem1 = function(z, overlap.group, k.groups){
  rho = 1
  iter.num =10
  Num.groups = length(overlap.group)
  # Initialize theta, y and u
  u = z
  y = theta =  list()
  for(i in 1:Num.groups){
    theta[[i]] =  y[[i]] = rep(0, length(overlap.group[[i]]))
  }
  
  for(iteration in 1:iter.num){
    # (a) Update u vector
    for(i in 1:Num.groups){
      u[overlap.group[[i]]] = z[overlap.group[[i]]] + theta[[i]] + rho*y[[i]]
    }
    u = u/norm(u,"F")
    
    # (b) Update y vector
    y.group.norm = NULL
    for(i in 1:Num.groups){
      y[[i]] = u[overlap.group[[i]]] - theta[[i]]/rho
      y.group.norm = c(y.group.norm, norm(matrix(y[[i]],ncol=1),"F"))
    }
    temp = order(y.group.norm,decreasing=T)
    active.group =temp[1:k.groups] 
    lambda.w = y.group.norm[temp[k.groups+1]]
    
    for(i in 1:Num.groups){
      if(i%in%active.group){
        y[[i]] = (1-lambda.w/y.group.norm[i])*y[[i]]
      }else 
        y[[i]] = 0*y[[i]]
    }
    
    # (c) Update theta vector
    for(i in 1:Num.groups){
      theta[[i]] = theta[[i]] + rho*(y[[i]]-u[overlap.group[[i]]])
    }
  }
  #print(active.group)
  active.membrs = NULL
  for(g in active.group){
    active.membrs = c(active.membrs,overlap.group[[g]])
  }
  #z = get.z(z,overlap.group,k.groups,node.degree)
  z[-active.membrs] = 0
  z = z/norm(z,"F")
}

# get.z = function(z,overlap.group,k.groups,node.degree){
#   z.group.norm = NULL
#   for(g in overlap.group){
#     z.group.norm = c(z.group.norm, norm(z[g],"2"))
#   }
#   lambda.w = order(z.group.norm,decreasing=T)[k.groups+1]
# 
#   z2 = 0*z
#   k = 1
#   for(group in overlap.group){
#     z2[group] = z2[group] + (1-lambda.w/max(z.group.norm[k],0.001))*z[group]
#     print(1-lambda.w/z.group.norm[k])
#     k = k+1
#   }
#   z2 = z2/node.degree
# }
# 
# get.node.degree = function(overlap.groups){
#   Temp = NULL
#   for(group in overlap.groups){
#     Temp = c(Temp,group)
#   }
#   members = unique(Temp)
#   node.degree = rep(0,length(members))
# 
#   for(g in overlap.groups){
#     node.degree[g] = node.degree[g]+1
#   }
#   return(node.degree)
# }



