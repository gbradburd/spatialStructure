 cov.func <- function(a0,a1,a2,geoDist,nuggets){
	geoCov <- a0 * exp(-(a1*geoDist)^a2)
	diag(geoCov) <- diag(geoCov) + nuggets
	return(geoCov)
}

sim.spStr.dataset <- function(k,n.loci){
	beta <- rbeta(1,1,5)
	gamma <- 1/(8*beta+4)
	zeta <- 0.25-1/(8*beta+4)
	epsilon <- rbeta(n.loci,beta,beta)
	nuggets <- rexp(k,1)
	coords <- cbind(runif(k),runif(k))
	geoDist <- fields::rdist(coords)
	a0 <- rexp(1)
	aD <- rexp(1)
	a2 <- runif(1,0.0001,2)
	geoCov <- cov.func(a0,aD,a2,geoDist,nuggets)
	true.freqs <- matrix(NA,nrow=k,ncol=n.loci)
	for(i in 1:n.loci){
		true.freqs[,i] <- MASS::mvrnorm(1,mu=rep(epsilon[i],k),
											Sigma=epsilon[i]*(1-epsilon[i])*geoCov)
	}
	bounded.freqs <- true.freqs
	bounded.freqs[which(bounded.freqs < 0)] <- 0
	bounded.freqs[which(bounded.freqs > 1)] <- 1
	sample.sizes <- sample(10:20,k,replace=TRUE)
	sample.size.mat <- matrix(sample.sizes,nrow=k,ncol=n.loci)
	obs.counts <- matrix(rbinom(k*n.loci,sample.size.mat,bounded.freqs),nrow=k,ncol=n.loci)
	obs.freqs <- obs.counts/sample.size.mat
	invars <- which(apply(obs.freqs,2,function(x){length(unique(x))})==1)
	dropped.freqs <- obs.freqs[,-invars]
	n.obs.loci <- ncol(dropped.freqs)
	dropped.cov <- cov(t(dropped.freqs))
	spStr.dataset <- list("allele.freqs" = dropped.freqs,
						  "population.coordinates" =  coords,
						  "sample.cov" = dropped.cov,
						  "D" = geoDist,
						  "sample.sizes" = sample.sizes,
						  "simCov" = geoCov,
						  "simEps" = epsilon,
						  "n.loci" = n.obs.loci,
						  "sim.pars" = list("beta"=beta,
						  					"gamma"=gamma,
						  					"zeta"=zeta,
						  					"geoCov"=geoCov,
						  					"a0"=a0,
					  						"aD"=aD,
					  						"a2"=a2,
					  						"nuggets"=nuggets,
					  						"epsilon"=epsilon))
	return(spStr.dataset)
}

run.test <- function(k,n.loci,model.block.file,dir,n.iter){
	model.block <- paste0(readLines(model.block.file,warn=FALSE),collapse="\n")
	spStr.dataset <- sim.spStr.dataset(k,n.loci)
	data.block <- list( "K" = 1,
						"N" = nrow(spStr.dataset$population.coordinates),
						"L" = (spStr.dataset$n.loci-1),
						"obsSigma" = spStr.dataset$sample.cov,
						"geoDist" = spStr.dataset$D,
						"sampleSize" = spStr.dataset$sample.sizes,
						"obsNugs" = spStr.dataset$sim.pars$nuggets)
	require(rstan)
	fit <- stan(model_code = model.block, data =data.block, iter=n.iter, chains = 1,thin=n.iter/500)
	save(spStr.dataset,data.block,fit,file=paste0(dir,"/test_output.Robj",collapse=""))

	gamma <- extract(fit,pars="gamma",inc_warmup=TRUE,permute=FALSE)
	zeta <- extract(fit,pars="zeta",inc_warmup=TRUE,permute=FALSE)
	alpha0 <- extract(fit,pars="alpha0",inc_warmup=TRUE,permute=FALSE)
	alphaD <- extract(fit,pars="alphaD",inc_warmup=TRUE,permute=FALSE)
	alpha2 <- extract(fit,pars="alpha2",inc_warmup=TRUE,permute=FALSE)
	nugget <- extract(fit,pars="nugget",inc_warmup=TRUE,permute=FALSE)
	par.cov <- geoStructure::get.par.cov(fit,1,data.block$N)
		
pdf(file=paste0(dir,"/all_plots.pdf",collapse=""),width=11,height=5)
	# quartz(width=11,height=5)
	x <- seq((length(alpha0)/2 + 1),length(alpha0))
	layout(mat=matrix(c(1,1,2,3,4,1,1,5,6,7),nrow=2,ncol=5,byrow=TRUE))
	par(mar=c(2,2,1,1))
		burnin <- 100
		n.iter <- dim(par.cov)[1]
	    index.mat <- upper.tri(data.block$obsSigma, diag = TRUE)
	    cov.range <- range(c(data.block$obsSigma, par.cov[(burnin + 1):n.iter, , ]))
	    plot(data.block$geoDist, data.block$obsSigma, xlab = "geographic distance", 
	        ylab = "covariance", ylim = cov.range, type = "n")
	    lapply(seq((burnin + 1), 500, length.out = min(100, 500 - 
	        (burnin + 1))), function(i) {
	        points(data.block$geoDist[index.mat], par.cov[i, 
	            , ][index.mat], pch = 20, col = adjustcolor("red", 
	            0.1))
	    		})
	    points(data.block$geoDist[index.mat], data.block$obsSigma[index.mat], 
	        xlab = "geographic distance", ylab = "covariance", ylim = cov.range, 
	        pch = 19)
		mtext("Cov/geoDist",side=3,padj=2)
	plot(gamma[x],
			ylim=range(c(gamma[x],spStr.dataset$sim.pars$gamma)))
		abline(h=spStr.dataset$sim.pars$gamma,col="red")
		mtext("gamma",side=3,padj=2)
	plot(spStr.dataset$sim.pars$nugget,apply(nugget,3,mean),pch=19,col="blue")
	    lapply(seq((burnin + 1), 500, length.out = min(100, 500 - (burnin + 1))), 
	   			function(i) {
	       			points(spStr.dataset$sim.pars$nugget, nugget[i, , ], 
	        					pch = 20, col = adjustcolor("black", 0.1))
		    		})
		abline(0,1,col="red")
	plot(zeta[x],
			ylim=range(c(zeta[x],spStr.dataset$sim.pars$zeta)))
		abline(h=spStr.dataset$sim.pars$zeta,col="red")
		mtext("zeta",side=3,padj=2)
	plot(alpha0[x],
			ylim=range(c(alpha0[x],spStr.dataset$sim.pars$a0)))
		abline(h=spStr.dataset$sim.pars$a0,col="red")
		mtext("a0",side=3,padj=2)
	plot(alphaD[x],
			ylim=range(c(alphaD[x],spStr.dataset$sim.pars$aD)))
		abline(h=spStr.dataset$sim.pars$aD,col="red")
		mtext("aD",side=3,padj=2)
	plot(alpha2[x],
			ylim=range(c(alpha2[x],spStr.dataset$sim.pars$a2)))
		abline(h=spStr.dataset$sim.pars$a2,col="red")
		mtext("a2",side=3,padj=2)
dev.off()
	return(invisible("test!"))
}
