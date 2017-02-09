admixed_covariance <- function(N, K, alpha0, alphaD, alpha2, geoDist, w, nugget, mu, sampleSize=NULL) {
	# recover()
	Sigma <- matrix(0,nrow=N,ncol=N)
		for (i in 1:(N-1)){
			for (j in (i+1):N){
				for(k in 1:K){
					Sigma[i, j] <- Sigma[i, j]  +  w[i,k] * w[j,k] * alpha0[k] * ( exp( -(alphaD[k]* geoDist[i, j])^alpha2[k] ) + mu[k]);
				}
				Sigma[j, i] <- Sigma[i, j];			
			}
		}
		for(i in 1:N){
			for(k in 1:K){
				Sigma[i, i] <- Sigma[i, i] + w[i,k]^2 * (alpha0[k] * (1 + mu[k]));
			}
			if(!is.null(sampleSize)){
				Sigma[i,i] <- Sigma[i,i] + (1 - Sigma[i,i]) * (1/sampleSize[i]);
			}
		}
		diag(Sigma) <- diag(Sigma) + nugget;
	return(Sigma);
}

truncate.freqs <- function(freqs){
	trunc.freqs <- freqs
	n.dropped <- length(which(freqs < 0 | freqs > 1))
	trunc.freqs[which(freqs < 0)] <- 0
	trunc.freqs[which(freqs > 1)] <- 1
	trunc.output <- list("trunc.freqs" = trunc.freqs,
						 "n.dropped" = n.dropped)
	return(trunc.output)
}

drop.invars <- function(freqs,epsilon,invars){
	if(length(invars!=0)){
		freqs <- freqs[,-invars]
		epsilon <- epsilon[-invars]
	} else {
		freqs <- freqs
		epsilon <- epsilon
	}
	return(list("freqs"=freqs,"epsilon"=epsilon))
}

get.pars <- function(N,K,n.loci){
	nuggets <- rep(0,N)
	coords <- cbind(runif(N),runif(N))
	geoDist <- fields::rdist(coords)
	a0 <- abs(rnorm(K,mean=0,sd=0.01))
	aD <- rexp(K)
	a2 <- runif(K,0.0001,2)
	mu <- abs(rnorm(K,mean=0,sd=0.01))
	admix.props <- gtools::rdirichlet(N,alpha=rep(0.1,K))
	par.list <- list("nuggets" = nuggets,
					 "coords" =	coords,
					 "geoDist" = geoDist,
					 "a0" =	a0,
					 "aD" =	aD,
					 "a2" =	a2,
					 "mu" =	mu,
					 "admix.props" = admix.props)
	return(par.list)
}

get.sim.draws <- function(N,n.loci,par.list,Sigma,sample.size){
	epsilon <- rnorm(n.loci,mean=0.5,sd=0.01)
	true.freqs <- matrix(NA,nrow=N,ncol=n.loci)
	for(i in 1:n.loci){
		true.freqs[,i] <- MASS::mvrnorm(1,mu=rep(epsilon[i],N),
											Sigma=Sigma)
	}
	trunc.output <- truncate.freqs(true.freqs)
	trunc.freqs <- trunc.output$trunc.freqs
	sample.size.mat <- matrix(sample.size,nrow=N,ncol=n.loci)
	obs.counts <- matrix(rbinom(N*n.loci,sample.size.mat,trunc.freqs),nrow=N,ncol=n.loci)
	obs.freqs <- obs.counts/sample.size.mat	
	invars <- which(apply(obs.freqs,2,function(x){length(unique(x))})==1)
	var.data <- drop.invars(obs.freqs,par.list$epsilon,invars)
	freqs <- var.data$freqs
#	epsilon <- var.data$epsilon
#	n.loci <- length(epsilon)
#	sample.true.cov <- cov(t(true.freqs))
	sample.trunc.cov <- cov(t(freqs))
	return(sample.trunc.cov)
}

get.wish.draws <- function(n.reps,n.loci,Sigma){
	require(MCMCpack)
	wish.draws <- lapply(1:n.reps,function(i){rwish(n.loci,Sigma/n.loci)})
	return(wish.draws)
}

mean.var.comp <- function(N,draws){
	means <- matrix(0,nrow=N,ncol=N)
	vars <- matrix(0,nrow=N,ncol=N)
	for(i in 1:N){
		for(j in i:N){
			means[i,j] <- mean(unlist(lapply(draws,"[[",i,j)))
				means[j,i] <- means[i,j]
			vars[i,j] <- var(unlist(lapply(draws,"[[",i,j)))
				vars[j,i] <- vars[i,j]
		}
	}
	return(list("means" = means,
				"vars" = vars))
}

sim.mean.var.comp <- function(n.reps,N,K,n.loci,sample.size){
	par.cov.list <- get.pars(N,K,n.loci)
	Sigma.wish <- admixed_covariance(N,K,par.cov.list$a0,par.cov.list$aD,par.cov.list$a2,
									 par.cov.list$geoDist,par.cov.list$admix.props,
									 par.cov.list$nugget,par.cov.list$mu,
									 sampleSize=rep(sample.size,N))
	Sigma.sim <- admixed_covariance(N,K,par.cov.list$a0,par.cov.list$aD,par.cov.list$a2,
									 par.cov.list$geoDist,par.cov.list$admix.props,
									 par.cov.list$nugget,par.cov.list$mu,
									 sampleSize=NULL)
	wish.draws <- get.wish.draws(n.reps,n.loci=n.loci,Sigma=Sigma.wish)
	sim.draws <- lapply(1:n.reps,function(i){get.sim.draws(N,n.loci,par.cov.list,Sigma.sim,sample.size)})
	wish.comp <- mean.var.comp(N,wish.draws)
	sim.comp <- mean.var.comp(N,sim.draws)
	im <- upper.tri(wish.comp[[1]],diag=TRUE)
	mean.var.comp.output <- list("par.cov.list" = par.cov.list,
								 "Sigma.wish" = Sigma.wish,
								 "Sigma.sim" = Sigma.sim,
								 "wish.draws" = wish.draws,
								 "sim.draws" = sim.draws,
								 "wish.comp" = wish.comp,
								 "sim.comp" = sim.comp,
								 "N" = N,
								 "K" = K,
								 "n.loci" = n.loci,
								 "sample.size" = sample.size)
	pdf(file=paste0("mean_var_comparison_nChromo=",sample.size,".pdf"),width=10,height=5)
		par(mfrow=c(1,2))
			plot(wish.comp$var[im],wish.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(W[i,j])",ylab="mean(W[i,j])",
					main="sims under wishart")
			plot(sim.comp$var[im],sim.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(U[i,j])",ylab="mean(U[i,j])",
					main="sims under model")
	dev.off()
	save(mean.var.comp.output,file=paste0("mean.var.comp.output_nChromo=",sample.size,".Robj"))
}

sim.mean.var.comp(n.reps=100,N=25,K=2,n.loci=1e4,sample.size=100)
sim.mean.var.comp(n.reps=100,N=25,K=2,n.loci=1e4,sample.size=10)
sim.mean.var.comp(n.reps=100,N=25,K=2,n.loci=1e4,sample.size=2)

load("/Users/gburd/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/wishart_mean_var/mean.var.comp.output_nChromo=2.Robj")
mean.var.comp.output2 <- mean.var.comp.output
load("/Users/gburd/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/wishart_mean_var/mean.var.comp.output_nChromo=10.Robj")
mean.var.comp.output10 <- mean.var.comp.output
load("/Users/gburd/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/wishart_mean_var/mean.var.comp.output_nChromo=100.Robj")
mean.var.comp.output100 <- mean.var.comp.output
im <- upper.tri(matrix(NA,nrow=25,ncol=25))

pdf(file="all_comparisons_NoDiag.pdf",width=6,height=8)
par(mfrow=c(3,2))
			plot(mean.var.comp.output2$wish.comp$var[im],mean.var.comp.output2$wish.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(W[i,j])",ylab="mean(W[i,j])",
					main="sims under wishart - samp.size = 2")
			plot(mean.var.comp.output2$sim.comp$var[im],mean.var.comp.output2$sim.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(U[i,j])",ylab="mean(U[i,j])",
					main="sims under model - samp.size = 2")
			plot(mean.var.comp.output10$wish.comp$var[im],mean.var.comp.output10$wish.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(W[i,j])",ylab="mean(W[i,j])",
					main="sims under wishart - samp.size = 10")
			plot(mean.var.comp.output10$sim.comp$var[im],mean.var.comp.output10$sim.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(U[i,j])",ylab="mean(U[i,j])",
					main="sims under model - samp.size = 10")
			plot(mean.var.comp.output100$wish.comp$var[im],mean.var.comp.output100$wish.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(W[i,j])",ylab="mean(W[i,j])",
					main="sims under wishart - samp.size = 100")
			plot(mean.var.comp.output100$sim.comp$var[im],mean.var.comp.output100$sim.comp$mean[im],
					pch=20,col=adjustcolor(1,0.5),
					xlab="var(U[i,j])",ylab="mean(U[i,j])",
					main="sims under model - samp.size = 100")
dev.off()