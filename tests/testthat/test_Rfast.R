n=1e4
m=1e2
Knm = matrix(rnorm(n*m),nrow=n)
Kmm = cov(matrix(rnorm(m*m),nrow=m))+diag(1e-5,m)
x = seq(0,1,length.out=n)
x_ind = seq(0,1,length.out=m)
library(microbenchmark)
microbenchmark(chol(Kmm),Rfast::cholesky(Kmm))
microbenchmark(outer(x,x_ind,"-"),Rfast::Outer(x,x_ind,"-"))
microbenchmark(outer(x,x_ind,"-"),outer(x,x_ind,"-")^2,times=1000)
microbenchmark(outer(x,x_ind,"-"),abs(outer(x,x_ind,"-")),times=1000)
microbenchmark(crossprod(Kmm),Rfast::Crossprod(Kmm,Kmm),times=1000)
microbenchmark(Knm%*%Kmm,Rfast::mat.mult(Knm,Kmm))
microbenchmark(Kmm%*%Kmm,Rfast::mat.mult(Kmm,Kmm))
(Knm%*%Kmm)[1:5,1:5]
(Rfast::mat.mult(Knm,Kmm))[1:5,1:5]
microbenchmark(Knm%*%x_ind,mat.mult(Knm,cbind(x_ind)))
microbenchmark(Knm/x_ind)
microbenchmark(colsums(Knm),colSums(Knm))
microbenchmark(t(Knm),transpose(Knm))
microbenchmark(crossprod(Knm,x),crossprod(x,Knm),t(Knm)%*%x)
microbenchmark(t(Kmm))
microbenchmark(Knm%*%t(Kmm),Rfast::Tcrossprod(Knm,Kmm))
Rfast::Tcrossprod(Knm,Kmm)[1:5,1:5]
microbenchmark(rowSums(Knm),Rfast::rowsums(Knm))
