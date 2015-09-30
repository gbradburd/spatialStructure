#validate
#put stan in tryCatch, email me
#message in place of cat()
#how to write a package that uses stan
check.list <- function(data.block){
	if(is.null(data.block$n.loci)){
		stop()
	}
#! "n.loci" %in% names(data.block) | 
}

check.n.samples <- function(data.block){
	n.samples <- nrow(data.block$sample.covariance)
	if(!is.null(data.block$geo.coords)){
		n.samples <- c(n.samples,nrow(data.block$geo.coords))
	}
	if(!is.null(data.block$time.coords)){
		n.samples <- c(n.samples,nrow(data.block$time.coords))
	}
	if(!is.null(data.block$mean.sample.sizes)){
		n.samples <- c(n.samples,length(data.block$mean.sample.sizes))
	}
	if(length(unique(n.samples)) > 1){
		good <- FALSE
	}
	return(good)
}

check.data <- function(data.block){
	cat("checking data \n")
	cat("checking data list structure")
	cat(sprintf("reading %s samples\n",nrow(data.block$geo.coords)))
	test.n.samples <- check.n.samples(data.block)
		if(!test.n.samples){
			stop("the number of samples is not consistent across objects in the data block\n")
		}
	cat("reading data block\n")
	cat("reading data block\n")
	cat("reading data block\n")
}

make.stan.data.block <- function(time,n.clusters){
	if(time & n.clusters == 1){
		stan.data.block <- 
		"data {
			int<lower=1> K;		  		// number of clusters
			int<lower=2> N; 	  		// number of samples
			int<lower=N+1> L;	    		// number of loci
		  	matrix[N,N] obsSigma; 		// observed covar
		  	matrix[N, N] geoDist; 		// matrix of pairwise geographic distance 
		  	matrix[N, N] timeDist; 		// matrix of pairwise temporal distance 
		  	vector[N] sampleSize;		// number of chromosomes genotyped for each sample
		  	real<lower=0> binVar;		// (mean.f * (1-mean.f)) averaged over all loci
		}"
	}
	if(time & n.clusters != 1){
		stan.data.block <- 
		"data {
			int<lower=1> K;		  		// number of clusters
			int<lower=2> N; 	  		// number of samples
			int<lower=N+1> L;	    		// number of loci
		  	matrix[N,N] obsSigma; 		// observed covar
		  	matrix[N, N] geoDist; 		// matrix of pairwise geographic distance 
		  	matrix[N, N] timeDist; 		// matrix of pairwise temporal distance 
		  	vector[K] DirichAlpha; 		// dirichlet prior parameters on admixture proportions
		  	vector[N] sampleSize;		// number of chromosomes genotyped for each sample
		  	real<lower=0> binVar;		// (mean.f * (1-mean.f)) averaged over all loci
		}"
	}
	if(!time & n.clusters == 1){
		stan.data.block <- 
		"data {
			int<lower=1> K;		  		// number of clusters
			int<lower=2> N; 	  		// number of samples
			int<lower=N+1> L;	    		// number of loci
		  	matrix[N,N] obsSigma; 		// observed covar
		  	matrix[N, N] geoDist; 		// matrix of pairwise geographic distance 
		  	vector[N] sampleSize;		// number of chromosomes genotyped for each sample
		  	real<lower=0> binVar;		// (mean.f * (1-mean.f)) averaged over all loci
		}"
	}
	if(!time & n.clusters != 1){
		stan.data.block <- 
		"data {
			int<lower=1> K;		  		// number of clusters
			int<lower=2> N; 	  		// number of samples
			int<lower=N+1> L;	    		// number of loci
		  	matrix[N,N] obsSigma; 		// observed covar
		  	matrix[N, N] geoDist; 		// matrix of pairwise geographic distance 
		  	vector[K] DirichAlpha; 		// dirichlet prior parameters on admixture proportions
		  	vector[N] sampleSize;		// number of chromosomes genotyped for each sample
		  	real<lower=0> binVar;		// (mean.f * (1-mean.f)) averaged over all loci
		}"
	}
	return(stan.data.block)
}

make.stan.par.block <- function(time,n.clusters){
	if(time & n.clusters == 1){
		par.block <- 
			"parameters {
				real<lower=0> alpha0;				// sill of the parametric covariance in cluster k
				real<lower=0> alphaD;				// effect of geographic distance in the parametric covariance in cluster k
				real<lower=0> alphaT;				// effect of temporal distance in the parametric covariance in cluster k
				real<lower=0, upper=2>  alpha2;		// exponential slope parameter in the parametric covariance in cluster k
				real<lower=0> gamma;					// global covariance due to shared ancestral mean between all samles
		  		real<lower=0> nugget[N]; 			// sample-specific variance (allele sampling error + sample-specific drift)
			}"
	}
	if(time & n.clusters != 1){
		par.block <- 
			"parameters {
				real<lower=0> alpha0[K];			// sill of the parametric covariance in cluster k
				real<lower=0> alphaD[K];			// effect of geographic distance in the parametric covariance in cluster k
				real<lower=0> alphaT[K];			// effect of temporal distance in the parametric covariance in cluster k
				real<lower=0, upper=2>  alpha2[K];	// exponential slope parameter in the parametric covariance in cluster k
				real<lower=0> mu[K];					// shared drift effect in cluster k
				real<lower=0> gamma;					// global covariance due to shared ancestral mean between all samles
			  	real<lower=0> nugget[N]; 			// sample-specific variance (allele sampling error + sample-specific drift)
				simplex[K]    w[N];    				// every sample (N in total) has a K simplex (i.e. K clusters) 
			}"
	}
	if(!time & n.clusters == 1){
		par.block <- 
			"parameters {
				real<lower=0> alpha0;				// sill of the parametric covariance in cluster k
				real<lower=0> alphaD;				// effect of geographic distance in the parametric covariance in cluster k
				real<lower=0, upper=2>  alpha2;		// exponential slope parameter in the parametric covariance in cluster k
				real<lower=0> gamma;					// global covariance due to shared ancestral mean between all samles
			  	real<lower=0> nugget[N]; 			// sample-specific variance (allele sampling error + sample-specific drift)
			}"
	}
	if(!time & n.clusters != 1){
		par.block <- 
			"parameters {
				real<lower=0> alpha0[K];			// sill of the parametric covariance in cluster k
				real<lower=0> alphaD[K];			// effect of geographic distance in the parametric covariance in cluster k
				real<lower=0, upper=2>  alpha2[K];	// exponential slope parameter in the parametric covariance in cluster k
				real<lower=0> mu[K];					// shared drift effect in cluster k
				real<lower=0> gamma;					// global covariance due to shared ancestral mean between all samles
			  	real<lower=0> nugget[N]; 			// sample-specific variance (allele sampling error + sample-specific drift)
				simplex[K]    w[N];    				// every sample (N in total) has a K simplex (i.e. K clusters) 
			}"
	}
	return(par.block)
}

make.stan.trans.par.block <- function(time,n.clusters){
	if(time & n.clusters == 1){
		trans.par.block <- 
		"transformed parameters {
		  matrix[N,N] Sigma;		// this specifies the parametric, admixed covariance matrix
			for (i in 1:N){
				for (j in 1:N){
					Sigma[i, j] <- gamma + binVar * (  alpha0 * exp(-(alphaD * geoDist[i, j] + alphaT * timeDist[i, j])^alpha2 ));  
					if(i==j){
						Sigma[i, i]<- Sigma[i, i] + binVar * (nugget[i] + 1/sampleSize[i]);
					}
				}
			}
		}"
	}
	if(time & n.clusters != 1){
		trans.par.block <- 
		"transformed parameters {
		  matrix[N,N] Sigma;		// this specifies the parametric, admixed covariance matrix
			for (i in 1:N){
				for (j in 1:N){
					Sigma[i, j] <- gamma;
					for(k in 1:K){
						Sigma[i, j] <- Sigma[i, j]  +  binVar * ( w[i,k] * w[j,k] * ( alpha0[k] * exp(-(alphaD[k]* geoDist[i, j] + alphaT * timeDist[i, j])^alpha2[k] ) + mu[k]));
					}
					if(i==j){
						Sigma[i, i]<- Sigma[i, i] + binVar * (nugget[i] + 1/sampleSize[i]);
					}
				}
			}
		}"
	}
	if(!time & n.clusters == 1){
		trans.par.block <- 
		"transformed parameters {
		  matrix[N,N] Sigma;		// this specifies the parametric, admixed covariance matrix
			for (i in 1:N){
				for (j in 1:N){
					Sigma[i, j] <- gamma + binVar *  ( alpha0 * exp(-(alphaD * geoDist[i, j])^alpha2 ));  
					if(i==j){
						Sigma[i, i] <- Sigma[i, i] + binVar * (nugget[i] + 1/sampleSize[i]);
					}
				}
			}
		}"
	}
	if(!time & n.clusters != 1){
		trans.par.block <- 
		"transformed parameters {
		  matrix[N,N] Sigma;		// this specifies the parametric, admixed covariance matrix
			for (i in 1:N){
				for (j in 1:N){
					Sigma[i, j] <- gamma;
					for(k in 1:K){
						Sigma[i, j] <- Sigma[i, j]  +  binVar * ( w[i,k] * w[j,k] * ( alpha0[k] * exp(-(alphaD[k]* geoDist[i, j])^alpha2[k] ) + mu[k]));
					}
					if(i==j){
						Sigma[i, i]<- Sigma[i, i] + binVar * (nugget[i] + 1/sampleSize[i]);
					}
				}
			}
		}"
	}
	return(trans.par.block)
}

make.stan.model.block <- function(time,n.clusters){
	if(time & n.clusters == 1){
		stan.model.block <- 
		"model {
			alpha0 ~ exponential(1);								// prior on alpha0
			alphaD ~ exponential(1);								// prior on alphaD
			alphaT ~ exponential(1);								// prior on alphaT
			alpha2 ~ exponential(1);								// prior on alpha2
			for(i in 1:N) nugget[i] ~ exponential(1);				// prior on nugget
			gamma ~ exponential(1);									// prior on global shared drift
			(L*obsSigma) ~ wishart(L,Sigma);						// likelihood function
		}"
	}
	if(time & n.clusters != 1){
		stan.model.block <- 
		"model {
			alpha0 ~ exponential(1);								// prior on alpha0
			alphaD ~ exponential(1);								// prior on alphaD
			alphaT ~ exponential(1);								// prior on alphaT
			alpha2 ~ exponential(1);								// prior on alpha2
			for(i in 1:N) nugget[i] ~ exponential(1);				// prior on nugget
			mu ~ exponential(1);									// prior on cluster shared drift
			gamma ~ exponential(1);									// prior on global shared drift
			for(i in 1:N) w[i] ~ dirichlet(DirichAlpha);			// prior on admixture proportions
			(L*obsSigma) ~ wishart(L,Sigma);						// likelihood function	
		}"
	}
	if(!time & n.clusters == 1){
		stan.model.block <- 
		"model {
			alpha0 ~ exponential(1);								// prior on alpha0
			alphaD ~ exponential(1);								// prior on alphaD
			alpha2 ~ exponential(1);								// prior on alpha2
			for(i in 1:N) nugget[i] ~ exponential(1);				// prior on nugget
			gamma ~ exponential(1);									// prior on global shared drift
			(L*obsSigma) ~ wishart(L,Sigma);						// likelihood function
		}"
	}
	if(!time & n.clusters != 1){
		stan.model.block <- 
		"model {
			alpha0 ~ exponential(1);								// prior on alpha0
			alphaD ~ exponential(1);								// prior on alphaD
			alpha2 ~ exponential(1);								// prior on alpha2
			for(i in 1:N) nugget[i] ~ exponential(1);				// prior on nugget
			mu ~ exponential(1);									// prior on cluster shared drift
			gamma ~ exponential(1);									// prior on global shared drift
			for(i in 1:N) w[i] ~ dirichlet(DirichAlpha);			// prior on admixture proportions
			(L*obsSigma) ~ wishart(L,Sigma);						// likelihood function	
		}"
	}
	return(stan.model.block)
}

make.stan.code.block <- function(time=FALSE,n.clusters){
	stan.data.block <- make.stan.data.block(time,n.clusters)
	stan.par.block <- make.stan.par.block(time,n.clusters)
	stan.trans.par.block <- make.stan.trans.par.block(time,n.clusters)
	stan.model.block <- make.stan.model.block(time,n.clusters)
	stan.code.block <- paste(stan.data.block, 
							 stan.par.block, 
							 stan.trans.par.block, 
							 stan.model.block,
							 sep="\n")
	return(stan.code.block)
}

geoStructure <- function(data.block,n.chains,n.iter,file.name){
	#checks on data block
	#checks on model specification
	#make.stan.code.block
	stan.block <- make.stan.code.block(data.block$timeDist,data.block$K)
		#checks for stan block first
	#run model
	require(rstan)
	model.fit <- stan(model_code = stan.block,
						data = data.block,
						iter = n.iter,
						chains = n.chains)
	#save fit obj
	save(model.fit,file=file.name)
	#return fit obj
	return(model.fit)
}


load("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/trichocarpa_dataset.Robj")
data.block <- list("K" = 2,
					"N" = nrow(trichocarpa.data$geo.coords),
					"L" = trichocarpa.data$n.loci,
					"obsSigma" = trichocarpa.data$sample.cov,
					"geoDist" = fields::rdist(trichocarpa.data$geo.coords),
					"timeDist" = TRUE,
					"DirichAlpha" = rep(0.1,2),
					"sampleSize" = trichocarpa.data$sample.sizes,
					"binVar" = trichocarpa.data$binom.var)

test2 <- geoStructure(data.block=data.block,n.chains=4,n.iter=10000,file.name="~/desktop/stanobj_k1.Robj")

ad.cov <- apply(extract(test2,"Sigma")$Sigma,c(2,3),mean)
plot(data.block$geoDist,data.block$obsSigma)
points(data.block$geoDist,ad.cov,col="red",pch=20)
plot(data.block$obsSigma,ad.cov) ; abline(0,1,col="red")
