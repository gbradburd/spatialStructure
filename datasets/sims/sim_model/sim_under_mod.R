covariance <- function(N,K,w,geoDist,alpha0,alphaD,alpha2,mu,nugget){
	Sigma <- matrix(NA,nrow=N,ncol=N)
	for (i in 1:(N-1)){
		for (j in (i+1):N){
			Sigma[i, j] <- 0;
			for(k in 1:K){
				Sigma[i, j] <- Sigma[i, j]  +  w[i,k] * w[j,k] * ( alpha0[k] * exp(-(alphaD[k]* geoDist[i, j])^alpha2[k] ) + mu[k]);
				Sigma[j, i] <- Sigma[i, j]
			}
		}
	}
	for(i in 1:N){
		Sigma[i, i] <- nugget[i];
		for(k in 1:K){
			Sigma[i, i] <- Sigma[i, i] + w[i,k]^2 * (alpha0[k] + mu[k]);
		}
	}
	return(Sigma)
}

cluster.cov <- function(geoDist,alpha0,alphaD,alpha2,mu){
	sigma.k <- alpha0 * exp(-(alphaD* geoDist)^alpha2 ) + mu
	return(sigma.k)
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

sim.spStr.dataset <- function(N,K,n.loci,seed){
#	recover()
	if(!is.null(seed)){
		set.seed(seed)
	}
	no.good <- TRUE
	while(no.good){
		epsilon <- rnorm(n.loci,mean=0.5,sd=0.01)	#rbeta(n.loci,0.5,0.5)
		nuggets <- rep(0,N)	#rexp(N,200)
		coords <- cbind(runif(N),runif(N))
		geoDist <- fields::rdist(coords)
		a0 <- rexp(K,100)
		aD <- rexp(K)
		a2 <- runif(K,0.0001,2)
		mu <- rexp(K,100)
		admix.props <- gtools::rdirichlet(N,alpha=rep(0.1,K))
		geoCov <- covariance(N,K, admix.props,geoDist,a0,aD,a2,mu,nuggets)
		cluster.covs <- lapply(1:K,function(k){cluster.cov(geoDist,a0[k],aD[k],a2[k],mu[k])})
		cluster.freqs <- lapply(1:K,function(k){
								apply(matrix(epsilon,nrow=1),2,function(e){MASS::mvrnorm(1,mu=rep(e,N),Sigma=e * (1-e) * cluster.covs[[k]])})
							})
		true.freqs <- matrix(NA,nrow=N,ncol=n.loci)
		for(i in 1:n.loci){
			true.freqs[,i] <- MASS::mvrnorm(1,mu=rep(epsilon[i],N),
												Sigma=epsilon[i]*(1-epsilon[i])*geoCov)
		}
		admix.freqs <- Reduce("+",lapply(1:K,function(k){apply(cluster.freqs[[k]],2,function(f){admix.props[,k] * f})}))
		trunc.output <- truncate.freqs(admix.freqs)
		trunc.freqs <- trunc.output$trunc.freqs		
		sample.sizes <- sample(10:20,N,replace=TRUE)
		sample.size.mat <- matrix(sample.sizes,nrow=N,ncol=n.loci)
		obs.counts <- matrix(rbinom(N*n.loci,sample.size.mat,trunc.freqs),nrow=N,ncol=n.loci)
		obs.freqs <- obs.counts/sample.size.mat	
		invars <- which(apply(obs.freqs,2,function(x){length(unique(x))})==1)
		if(length(invars) == n.loci | trunc.output$n.dropped/n.loci > 0.05){
			no.good <- TRUE
		} else {
			no.good <- FALSE
		}
	}
	var.data <- drop.invars(obs.freqs,epsilon,invars)
	freqs <- var.data$freqs
	epsilon <- var.data$epsilon
	n.loci <- length(epsilon)
	sample.true.cov <- cov(t(true.freqs))
	sample.trunc.cov <- cov(t(freqs))
	cov.list <- list("sample.true.cov" = sample.true.cov,
					 "sample.trunc.cov" = sample.trunc.cov)
	freq.list <- list("freqs" = freqs,
					  "epsilon" = epsilon,
					  "true.freqs" = true.freqs,
					  "trunc.freqs" = trunc.freqs,
					  "admix.freqs" = admix.freqs)
	spStr.dataset <- list("population.coordinates" = coords,
						  "sample.sizes" = sample.sizes,
						  "D" = geoDist,
						  "n.loci" = n.loci,
						  "cov.list" = cov.list,
						  "freq.list" = freq.list,
						  "sim.pars" = list("geoCov"=geoCov,
						  					"a0"=a0,
					  						"aD"=aD,
					  						"a2"=a2,
					  						"mu"=mu,
					  						"admix.props" = admix.props,
					  						"nuggets"=nuggets,
					  						"seed"=seed))
	return(spStr.dataset)
}

plot.admix.prop.post <- function(K,N,admix.props,sim.admix.props,K.cols,cluster.order){
#	recover()
	qs <- lapply(1:K,function(i){
				apply(admix.props[,,i],2,
					function(x){quantile(x,c(0.025,0.975))})
			})
	mean.props <- lapply(1:K,function(i){
						apply(admix.props[,,i],2,
							function(x){mean(x)})
					})
	cluster.cols <- K.cols[cluster.order[,2]]
	plot(0,ylim=c(0,1),type='n',xlim=c(0,N),ylab="admix proportion",xlab="Sample")
		lapply(1:K,function(i){segments(x0=1:N,x1=1:N,
									    y0=apply(qs[[i]],2,min),
									    y1=apply(qs[[i]],2,max),
										col=cluster.cols[i],lwd=2)})
		lapply(1:K,function(i){points(mean.props[[i]],col=cluster.cols[i])})
		lapply(1:K,function(i){points(sim.admix.props[,i],col=K.cols[i],pch=18)})
	return(invisible("plotted"))
}