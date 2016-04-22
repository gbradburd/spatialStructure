################################
#	covariance of sample frequencies
################################

## PART 1: avoiding absorbing frequency bounds
#number of populations
	k <- 20
#number of loci
	n.loci <- 1e4
#simulate ancestral allele frequencies
	epsilon <- rbeta(n.loci,0.1,0.1)

#generate a population covariance structure based on geography
	cov.func <- function(a0,a1,a2,geoDist){
		geoCov <- a0 * exp(-(a1*geoDist)^a2)
		return(geoCov)
	}
	coords <- cbind(runif(20),runif(20))
	geoDist <- fields::rdist(coords)
	geoCov <- cov.func(0.1,1.9,0.7,geoDist)
	
#simulate population frequencies around the ancestral
#	allele frequencies with a covariance structure dictated 
#	by the geographically determined covariance
require(MASS)
true.freqs <- matrix(NA,nrow=k,ncol=n.loci)
for(i in 1:n.loci){
	true.freqs[,i] <- mvrnorm(1,mu=rep(epsilon[i],k),Sigma=epsilon[i]*(1-epsilon[i])*geoCov)
}

#bound these draws at 0 and 1
bounded.freqs <- true.freqs
bounded.freqs[which(bounded.freqs < 0)] <- 0
bounded.freqs[which(bounded.freqs > 1)] <- 1

#simulate population sampling w/ heterogeneous sample size
sample.sizes <- sample(2:10,20,replace=TRUE)
sample.size.mat <- matrix(sample.sizes,nrow=k,ncol=n.loci)
obs.counts <- matrix(rbinom(k*n.loci,sample.size.mat,bounded.freqs),nrow=k,ncol=n.loci)
obs.freqs <- obs.counts/sample.size.mat

#get various covariances
true.freq.cov <- cov(t(true.freqs))
bounded.freq.cov <- cov(t(bounded.freqs))
obs.cov <- cov(t(obs.freqs))

#plot results
#	note that all frequency-based covariances are higher
#	than the geographic covariance, as a result of the 
#	shared ancestral allele frequency
plot(geoDist,geoCov,pch=20,main="covariance against distance",xlab="geo distance",ylab="covariance",
		ylim=c(range(c(geoCov,true.freq.cov,bounded.freq.cov,obs.cov))))
	points(geoDist,true.freq.cov,pch=20,col="red")
	points(geoDist,bounded.freq.cov,pch=20,col="purple")
	points(geoDist,obs.cov,pch=20,col="blue")
legend(x="right",pch=19,col=c("black","red","purple","blue"),
		legend=c("geoCov","true freq cov","bounded freq cov","observed freq cov"))

#determine the relationship between geoCov and the observed covariance
var.eps <- var(epsilon)
binVar <- mean(epsilon*(1-epsilon))

#currently, our formula is, for a single cluster
	# transformed parameters {
			  # matrix[N,N] Sigma;		// this specifies the parametric, admixed covariance matrix
				# for (i in 1:N){
					# for (j in 1:N){
							# Sigma[i, j] <- gamma  +  binVar * ( alpha0 * exp(-(alphaD* geoDist[i, j])^alpha2 ));
						# if(i==j){
							# Sigma[i, i]<- Sigma[i, i] + binVar * (nugget[i] + 1/sampleSize[i]);
						# }
					# }
				# }
			# }
# so that would look like this
par(mfrow=c(1,3))
plot(geoDist,obs.cov,pch=20,main="covariance against distance",
		xlab="geo distance",ylab="covariance",
		ylim=range(c(obs.cov,c(var.eps+binVar*(geoCov+diag(1/sample.sizes))))))
		points(geoDist,var.eps + binVar*(geoCov + diag(1/sample.sizes)),pch=20,col="red")
index.mat <- lower.tri(geoDist)
plot(geoDist[index.mat],obs.cov[index.mat],pch=20,main="covariance against distance, just off-diagonals",
		xlab="geo distance",ylab="covariance",
		ylim=range(c(obs.cov[index.mat],c(var.eps+binVar*geoCov[index.mat]))))
		points(geoDist[index.mat],var.eps + binVar*(geoCov[index.mat]),pch=20,col="red")
plot(diag(obs.cov),var.eps + binVar*(diag(geoCov)+1/sample.sizes),
		main="just the diagonals now",xlab="diag(obs.cov)",ylab="binVar*(diag(geoCov)+diag(1/sample.sizes))")
		abline(0,1,col="red")


# #but I think that we should only be multiplying the diagonal elements by binVar
# #	so that would look like this:
# par(mfrow=c(1,3))
# plot(geoDist,obs.cov,pch=20,main="covariance against distance",
		# xlab="geo distance",ylab="covariance",
		# ylim=range(c(obs.cov,c(var.eps + (geoCov + diag(binVar/sample.sizes))))))
		# points(geoDist,var.eps + geoCov + diag(binVar/sample.sizes),pch=20,col="green4")
# index.mat <- lower.tri(geoDist)
# plot(geoDist[index.mat],obs.cov[index.mat],pch=20,main="covariance against distance, just off-diagonals",
		# xlab="geo distance",ylab="covariance",
		# ylim=range(c(obs.cov[index.mat],c(var.eps + geoCov[index.mat]))))
		# points(geoDist[index.mat],var.eps + geoCov[index.mat],pch=20,col="green4")
# plot(diag(obs.cov),var.eps + diag(geoCov) + binVar/sample.sizes,
		# main="just the diagonals now",xlab="diag(obs.cov)",ylab="binVar*(diag(geoCov)+diag(1/sample.sizes))")
		# abline(0,1,col="green4")
		
# # which seems like a better fit than the first model


# ##### INTERLUDE ######

# ## PART 2: not avoiding absorbing frequency bounds
# #number of populations
	# k <- 20
# #number of loci
	# n.loci <- 1e4
# #simulate ancestral allele frequencies
	# epsilon <- rbeta(n.loci,0.1,0.1)

# #generate a population covariance structure based on geography
	# cov.func <- function(a0,a1,a2,geoDist){
		# geoCov <- a0 * exp(-(a1*geoDist)^a2)
		# return(geoCov)
	# }
	# coords <- cbind(runif(20),runif(20))
	# geoDist <- fields::rdist(coords)
	# geoCov <- cov.func(0.01,1.9,0.7,geoDist)
	
# #simulate population frequencies around the ancestral
# #	allele frequencies with a covariance structure dictated 
# #	by the geographically determined covariance
# require(MASS)
# true.freqs <- matrix(NA,nrow=k,ncol=n.loci)
# for(i in 1:n.loci){
	# true.freqs[,i] <- mvrnorm(1,mu=rep(epsilon[i],k),Sigma=geoCov)
# }

# #bound these draws at 0 and 1
# bounded.freqs <- true.freqs
# bounded.freqs[which(bounded.freqs < 0)] <- 0
# bounded.freqs[which(bounded.freqs > 1)] <- 1

# #simulate population sampling w/ heterogeneous sample size
# sample.sizes <- sample(2:10,20,replace=TRUE)
# sample.size.mat <- matrix(sample.sizes,nrow=k,ncol=n.loci)
# obs.counts <- matrix(rbinom(k*n.loci,sample.size.mat,bounded.freqs),nrow=k,ncol=n.loci)
# obs.freqs <- obs.counts/sample.size.mat

# #get various covariances
# true.freq.cov <- cov(t(true.freqs))
# bounded.freq.cov <- cov(t(bounded.freqs))
# obs.cov <- cov(t(obs.freqs))

# #plot results
# #	note that all frequency-based covariances are higher
# #	than the geographic covariance, as a result of the 
# #	shared ancestral allele frequency
# #
# #	also note that now the non-bounded covariance is higher than the bounded covariance
# plot(geoDist,geoCov,pch=20,main="covariance against distance",xlab="geo distance",ylab="covariance",
		# ylim=c(range(c(geoCov,true.freq.cov,bounded.freq.cov,obs.cov))))
	# points(geoDist,true.freq.cov,pch=20,col="red")
	# points(geoDist,bounded.freq.cov,pch=20,col="purple")
	# points(geoDist,obs.cov,pch=20,col="blue")
# legend(x="right",pch=19,col=c("black","red","purple","blue"),
		# legend=c("geoCov","true freq cov","bounded freq cov","observed freq cov"))

# #determine the relationship between geoCov and the observed covariance
# var.eps <- var(epsilon)
# binVar <- mean(epsilon*(1-epsilon))

# #currently, our formula is, for a single cluster
	# # transformed parameters {
			  # # matrix[N,N] Sigma;		// this specifies the parametric, admixed covariance matrix
				# # for (i in 1:N){
					# # for (j in 1:N){
							# # Sigma[i, j] <- gamma  +  binVar * ( alpha0 * exp(-(alphaD* geoDist[i, j])^alpha2 ));
						# # if(i==j){
							# # Sigma[i, i]<- Sigma[i, i] + binVar * (nugget[i] + 1/sampleSize[i]);
						# # }
					# # }
				# # }
			# # }
# # so that would look like this
# par(mfrow=c(2,3))
# plot(geoDist,obs.cov,pch=20,main="covariance against distance",
		# xlab="geo distance",ylab="covariance",
		# ylim=range(c(obs.cov,c(var.eps+binVar*(geoCov+diag(1/sample.sizes))))))
		# points(geoDist,var.eps + binVar*(geoCov+diag(1/sample.sizes)),pch=20,col="red")
# index.mat <- lower.tri(geoDist)
# plot(geoDist[index.mat],obs.cov[index.mat],pch=20,main="covariance against distance, just off-diagonals",
		# xlab="geo distance",ylab="covariance",
		# ylim=range(c(obs.cov[index.mat],c(var.eps+binVar*geoCov[index.mat]))))
		# points(geoDist[index.mat],var.eps + binVar*(geoCov[index.mat]),pch=20,col="red")
# plot(diag(obs.cov),binVar*(diag(geoCov)+1/sample.sizes),
		# main="just the diagonals now",xlab="diag(obs.cov)",ylab="binVar*(diag(geoCov)+diag(1/sample.sizes))")
		# abline(0,1,col="red")


# #but I think that we should only be multiplying the diagonal elements by binVar
# #	so that would look like this:
# plot(geoDist,obs.cov,pch=20,main="covariance against distance",
		# xlab="geo distance",ylab="covariance",
		# ylim=range(c(obs.cov,c(var.eps + (geoCov + diag(binVar/sample.sizes))))))
		# points(geoDist,var.eps + geoCov + diag(binVar/sample.sizes),pch=20,col="green4")
# index.mat <- lower.tri(geoDist)
# plot(geoDist[index.mat],obs.cov[index.mat],pch=20,main="covariance against distance, just off-diagonals",
		# xlab="geo distance",ylab="covariance",
		# ylim=range(c(obs.cov[index.mat],c(var.eps + geoCov[index.mat]))))
		# points(geoDist[index.mat],var.eps + geoCov[index.mat],pch=20,col="green4")
# plot(diag(obs.cov),diag(geoCov)+binVar/sample.sizes,
		# main="just the diagonals now",xlab="diag(obs.cov)",ylab="binVar*(diag(geoCov)+diag(1/sample.sizes))")
		# abline(0,1,col="green4")

# # this looks admittedly worse, but I think this has to do with how much drift we're
# #	missing in the observed covariance due to the absorbing bounds, and I'm not sure
# #	the best way to fix that.

# # and, the green lines still look like a better description of the data than the 
# #	red ones, I think



