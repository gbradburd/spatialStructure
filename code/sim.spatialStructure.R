#SIMULATE TOY DATA
source("code/SpaceTimeStructureMix.R")

generate.sim.param.list <- function(n.clusters,n.samples,time.sampling){
	sim.param.list <- list( "admix.props" = gtools::rdirichlet(n = n.samples,alpha = rep(1,n.clusters)),
							"cov.par1"= replicate(n.clusters,rexp(1,rate=100)),
							"cov.par2"= replicate(n.clusters,rexp(1,rate=10)),
							"cov.par3"= replicate(n.clusters,ifelse(time.sampling,runif(1),0)),
							"cluster.mean" = replicate(n.clusters,rexp(1,rate=20)),
							"shared.mean" = rexp(1,rate=50),
							"nuggets" = rexp(n.samples,rate=100))
	return(sim.param.list)
}

populate.sim.cluster.list <- function(sim.cluster.list,sim.param.list,sampling.data){
	for(i in seq_along(sim.cluster.list)){
		sim.cluster.list[[i]]$covariance.params$cov.par1 <- sim.param.list$cov.par1[i]
		sim.cluster.list[[i]]$covariance.params$cov.par2 <- sim.param.list$cov.par2[i]
		sim.cluster.list[[i]]$covariance.params$cov.par3 <- sim.param.list$cov.par3[i]
		sim.cluster.list[[i]]$cluster.mean <- sim.param.list$cluster.mean[i]
		sim.cluster.list[[i]]$admix.prop.matrix <- sim.param.list$admix.props[,i] %*% t(sim.param.list$admix.props[,i])
		sim.cluster.list[[i]]$covariance <- cluster.covariance(sampling.data$geo.dist,
																sampling.data$time.dist,
																sim.cluster.list[[i]]$covariance.params)
	}
	return(sim.cluster.list)
}

simulate.freq.data <- function(n.loci,n.samples,sim.admixed.cov.mat){
	freq.data <- t(mvrnorm(n=n.loci,mu=numeric(n.samples),Sigma=sim.admixed.cov.mat))
	return(freq.data)
}

simulate.spatialStructure.dataset <- function(n.clusters,n.samples,n.loci,file.name,time.sampling=FALSE,sim.seed=NULL,counts=FALSE){
# recover()
# set seed
	if(is.null(sim.seed)){
		sim.seed <- sample(1:1e5,1)
	}
	cat(sim.seed)
	set.seed(sim.seed)
# simulate locations
	sampling.data <- list("geo.coords" = cbind(runif(n.samples),runif(n.samples)),
							"temporal.coords" = ifelse(rep(time.sampling,n.samples),
													c(sample(1:100, n.samples,replace=TRUE)),
													0))
	sampling.data$geo.dist <- fields::rdist(sampling.data$geo.coords)
	sampling.data$time.dist <- fields::rdist(sampling.data$temporal.coords)
	sim.param.list <- generate.sim.param.list(n.clusters,n.samples,time.sampling)
	sim.cluster.list <- lapply(seq_along(1:n.clusters),function(i){declare.cluster.list()})
	sim.cluster.list <- populate.sim.cluster.list(sim.cluster.list,sim.param.list,sampling.data)
	sim.admixed.cov.mat <- admixed.covariance(sim.cluster.list,n.clusters,sim.param.list$shared.mean,sim.param.list$nuggets)
	sim.data <- list("n.samples" = n.samples,"n.loci" = n.loci,"sampling.data" = sampling.data,
							"sim.param.list" = sim.param.list,"sim.cluster.list" = sim.cluster.list,
							"sim.admixed.cov.mat" = sim.admixed.cov.mat)
	if(counts){
		sim.admixed.cov.mat <- sim.admixed.cov.mat - sim.param.list$shared.mean
		ancestral.freqs <- rnorm(n.loci,mean=0.5,sd=0.01)
		sim.freqs <- simulate.freq.data(n.loci,n.samples,sim.admixed.cov.mat)
		sim.freqs <- sim.freqs + matrix(ancestral.freqs,nrow=n.samples,ncol=n.loci,byrow=TRUE)
		sim.freqs[which(sim.freqs < 0 )] <- 0
		sim.freqs[which(sim.freqs > 1 )] <- 1
		sim.counts <- matrix(rbinom(n.samples*n.loci,100,sim.freqs),nrow=n.samples,ncol=n.loci)
		sample.covariance <- cov(t(sim.counts/100))
		sim.data$ancestral.freqs <- ancestral.freqs
		sim.data$sim.freqs <- sim.freqs
		sim.data$sim.counts <- sim.counts
		sim.data$sample.covariance <- sample.covariance
	}
	save(sim.data,file=paste(file.name,".Robj",sep=""))
}

simulate.spatialStructure.dataset(n.clusters = 3,
									n.samples = 100,
									n.loci = 10000,
									file.name = "~/desktop/test",
									counts = TRUE)
load("~/desktop/test.Robj")
plot(sim.data$sim.admixed.cov.mat+var(sim.data$ancestral.freqs),sim.data$sample.covariance) ; 
abline(0,1,col="red")