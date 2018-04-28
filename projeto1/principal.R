Posdef = function (n, ev = runif(n, 0, 1)) {
  Z = matrix(ncol=n, rnorm(n^2))
  decomp = qr(Z)
  Q = qr.Q(decomp) 
  R = qr.R(decomp)
  d = diag(R)
  ph = d / abs(d)
  O = Q %*% diag(ph)
  Z = t(O) %*% diag(ev) %*% O
  return(Z)
}

hilbert = function(n){
  H = matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      H[i,j] = 1 /(i + j - 1)
    }
  }
  return(H)
}

source("resolve.R")
n = 1000
A = Posdef(n)
b = rep(1,n)
H = hilbert(n) + diag(n) * 10 ^ (- 14)
x_chol_A = resolve(A,b)
x_chol_H = resolve(H,b)

x_comando_A = solve (A,b)
x_comando_H = solve (H,b)

norma_chol_A = norm(b - A*x_chol_A)/norm(x_chol_A, type =  "2")
norma_chol_H = norm(b - H*x_chol_H)/norm(x_chol_H, type =  "2")
norma_comando_A = norm(b - A*x_comando_A)/norm(x_comando_A, type =  "2")
norma_comando_H = norm(b - H*x_comando_H)/norm(x_comando_H, type =  "2")
