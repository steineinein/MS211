resolve = function(A, b){
  G = cholesky(A)
  Gt = t(G)
  n = length(b)
  x = numeric(length = n)
  y = numeric(length = n)
  y[1] = b[1]/G[1,1]
  if (n > 1){
    for (i in 2:n){
      y[i] = (b[i] - sum(y[1:i-1]*G[i,1:i-1]))/G[i,i]
    }
  }
  x[n] = y[n]/Gt[n,n]
  if (n > 1){
    for (i in (n-1):1){
      x[i] = (y[i] - sum(x[(i+1):n]*Gt[i,(i+1):n]))/Gt[i,i]
    }
  }
  return (x)
}

cholesky = function(A){
  n_row = nrow(A)
  n_col = ncol(A)
  G = matrix(0, nrow = n_row, ncol = n_col)
  
  G[1,1] = sqrt(A[1,1])
  for (i in 2:n_row){
    G[i,1] = A[i,1]/G[1,1]
  }
  
  for (i in 2:n_row){
    for (j in 2:i){
      if (i == j){
        G[i,j] = sqrt(A[i,i] - sum(G[j,1:(i-1)]^2))
      } else {
        G[i,j] = (A[i,j] - sum(G[i,1:(i-1)]*G[j,1:(i-1)]))/G[j,j]
      }
    }
  }
  return (G)
}