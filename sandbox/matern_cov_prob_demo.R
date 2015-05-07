#Spatial Covariance Functions
matern <- function(r,nu){
	r^nu * gsl::bessel_Knu(nu,r)
}

spatial.covariance <- function(geo.dist,time.dist,cov.par1,cov.par2,cov.par3,d=2){
	# recover()
	((pi^(d/2)* cov.par1^d) / 
		((2^(cov.par2 + cov.par3*time.dist-1)) * 
			gsl::gsl_sf_gamma(cov.par2 + cov.par3*time.dist + d/2)) ) * 
	matern(cov.par1*geo.dist, cov.par2 + cov.par3*time.dist) #<- this is whence the badness comes
}

spatiotemporal.covariance <- function(geo.dist,time.dist,alpha,nu,zeta){
    d <- 2  # two spatial dimensions
    return( exp(
            d/2 * log(pi) +
            d * log(alpha) -
            ( (nu + zeta*time.dist - 1)*log(2) +
                    gsl::lngamma(nu + zeta*time.dist + d/2) 
                ) +
            nu * log(alpha*geo.dist) +
            gsl::bessel_lnKnu(nu=nu,x=alpha*geo.dist)
    ) )
}


#simulate spatial and temporal sampling
set.seed(11113)
k <- 30
n.loci <- 1e6
spatial.coords <- cbind(runif(k),runif(k))
temporal.coords <- rep(0,k)
geo.dist <- fields::rdist(spatial.coords)
time.dist <- fields::rdist(temporal.coords)

#parameterize the spatial covariance with values that
#	induce NA values on the diagonal of the matrix
problem.cov1 <- spatial.covariance(geo.dist = geo.dist,
									time.dist = time.dist,
									cov.par1 = 0.007757018,
									cov.par2 = 23.20698,
									cov.par3 = 1)
problem.cov2 <- spatiotemporal.covariance(geo.dist = geo.dist,
									time.dist = time.dist,
									alpha = 0.007757018,
									nu = 23.20698,
									zeta = 1)
diag(problem.cov1)
diag(problem.cov2)

#look around this parameter space for problem boundaries
cov.par1.seq <- seq(1e-10,0.01,length.out=100)
cov.par2.seq <- seq(20,30,length.out=100)
the.badness <- matrix(0,nrow=100,ncol=100)
for(i in 1:100){
	for(j in 1:100){
		spatcov <- spatial.covariance(geo.dist,
										time.dist,
										cov.par1.seq[i],
										cov.par2.seq[j],
										1)
		if(any(is.na(spatcov))){
			the.badness[i,j] <- 1
		}
	}
}

#visualize parameter space that makes problems
image(t(the.badness),axes=FALSE)
	axis(1,at=seq(0,1,length.out=10),round(cov.par1.seq[seq(1,100,length.out=10)],4),cex.axis=0.7)
	axis(2,at=seq(0,1,length.out=10),round(cov.par2.seq[seq(1,100,length.out=10)],4),cex.axis=0.7)
	legend(x="topright",pch=22,col=1,pt.bg=c("beige","red"),legend=c("no NAs","NAs"))