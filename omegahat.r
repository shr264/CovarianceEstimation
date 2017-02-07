library(gRbase)
library(Rgraphviz)
library(ggm)
library(Matrix)
library(MASS)
library(huge)
library(matrixcalc)
library(xtable)

(A = Matrix(runif(n=16,0,1), nrow = 4, ncol = 4))
(U = svd(A)$u)
(D = diag(svd(A)$d))
(A = U%*%D%*%t(U))
A[1,3] = 0
A[3,1] = 0
A[2,1] = A[1,2]
A[4,3] = A[3,4]
A[4,1] = A[1,4]
A[3,2] = A[2,3]
A[2,4] = A[4,2]
(omega = A)
chol(S)

(S = solve(omega))
svd(S)

SC1inv = S
SC1inv[] = 0
SC2inv = S
SC2inv[] = 0
SR2inv = S
SR2inv[] = 0
(SC1inv[c(1,2,4),c(1,2,4)] = solve(S[c(1,2,4),c(1,2,4)]))
(SC2inv[c(2,3,4),c(2,3,4)] = solve(S[c(2,3,4),c(2,3,4)]))
(SR2inv[c(2,4),c(2,4)] = solve(S[c(2,4),c(2,4)]))
omegahat = omega
omegahat[]=0
(omegahat = SC1inv + SC2inv - SR2inv)
solve(omegahat)



x=xtable(SR2inv,align=rep("",ncol(omega)+1)) # We repeat empty string 6 times
print(x, floating=FALSE, tabular.environment="pmatrix", 
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

x=xtable(solve(S[c(1,2,4),c(1,2,4)]),align=rep("",ncol(S[c(1,2,4),c(1,2,4)])+1)) # We repeat empty string 6 times
print(x, floating=FALSE, tabular.environment="pmatrix", 
  hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
