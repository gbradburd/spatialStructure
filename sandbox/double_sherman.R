library(microbenchmark)
library(RcppEigen)
library(Rcpp)

EPS <- 1e-8

double.sherman_r <- function(Ap, u, v,i ){   #, faster=TRUE) {

#	if(!faster) Sherman.correction<-(Ap %*% u %*% t(v) %*% Ap)  /drop(1 + t(v) %*% Ap %*% u)
	
	 ##  t(v) %*% Ap  =  Ap[i,]
	
	
	determinant.correction<- (1 +   sum(Ap[i,]*u))	  #    http://en.wikipedia.org/wiki/Matrix_determinant_lemma
	 Sherman.correction<-(Ap %*% u %*% Ap[i,])/determinant.correction   ##     (Ap u vT Ap)/ (1 + t(v) * Ap * u)

	
	
  return(list( determinant.correction = determinant.correction , new.inverse = Ap - Sherman.correction -  t(Sherman.correction)))
}


dummylikelihood<-function(A,B){
	Ai<-solve(A)
	 ##sum(diag( Ai %*% B )) *determinant(A,logarithm=TRUE)$modulus
}


test.double.sherman<-function(){
	m <- matrix(sample(as.double(1:100000000), 1e6), nrow=1000)
	ms <- solve(m)
	
	F.mat.k<-matrix(sample(as.double(1:100000000), 1e6), nrow=1000)
	F.mat.k.prime<-matrix(sample(as.double(1:100000000), 1e6), nrow=1000)
	
	w.k<-runif(1000)
	w.k.prime<-runif(1000)
	
	##change the inverse matrix to reflect  wi^k -> wi^k + delta.w & wi^k' -> wi^k' - delta.w
	##fast.update.inverse.w

	delta.w<-0.1
	i=10
	##note that we can pull the change is wi, delta.w out of the 
	change.vec<-delta.w * ( w.k * F.mat.k[i,]  - w.k.prime * F.mat.k.prime[i,] )
	
	u <- change.vec
	v <- rep(0,1000)	
	v[i] <- 1
	new.inv<-sherman(ms, u, v)	
	new.inv<-sherman(new.inv, v, u)	

    another.way<-doublesherman(ms, u, v)

	new.inv.old.way<- solve( m+ v %*% t(u) + u %*% t(v))	
    all(abs(new.inv.old.way - new.inv)< EPS)
	all(abs(another.way - new.inv)< EPS)



	use_solve <- function(m, u, v) {
	  solve(m + u %*% t(v) + v %*% t(u))
	}
	
	use_doublesherman <- function(ms, u, v) {
	  doublesherman(ms, u, v)
	}
	
	use_doublesherman_r <- function(ms, u, v,faster) {
	  double.sherman_r(ms, u, v,faster)
	}
	
	res <- microbenchmark(cpp_dsherman={ use_doublesherman(ms, u, v) },rfast_dsherman= {use_doublesherman_r(ms, u, v,faster=TRUE)},r_dsherman= {use_doublesherman_r(ms, u, v,faster=FALSE)},
	                      solve={ use_solve(m, u, v) }, times=30)
                      
}