#SIMULATE TOY DATA
source("~/Dropbox/InspectorSpaceTime/spatialStructure/code/SpaceTimeStructureMix.R")

generate.sim.param.list <- function(n.clusters,n.samples,sampling.data,time.sampling){
	sim.param.list <- list( "admix.props" = gtools::rdirichlet(n = n.samples,alpha = rep(0.5,n.clusters)),
							"cov.par1"= replicate(n.clusters,rexp(1,rate=100)),
							"cov.par2"= replicate(n.clusters,rnorm(1,2*mean(sampling.data$geo.dist),sd=3)),
							"cov.par3"= replicate(n.clusters,runif(1,0.01,1)),
							"cov.par4"= replicate(n.clusters,ifelse(time.sampling,rnorm(1,2*mean(sampling.data$geo.dist),sd=3),0)),
							"cluster.mean" = replicate(n.clusters,rexp(1,rate=100)),
#							"shared.mean" = rexp(1,rate=100),
							"nuggets" = rexp(n.samples,rate=100))
	return(sim.param.list)
}

populate.sim.cluster.list <- function(sim.cluster.list,sim.param.list,sampling.data){
	for(i in seq_along(sim.cluster.list)){
		sim.cluster.list[[i]]$covariance.params$cov.par1 <- sim.param.list$cov.par1[i]
		sim.cluster.list[[i]]$covariance.params$cov.par2 <- sim.param.list$cov.par2[i]
		sim.cluster.list[[i]]$covariance.params$cov.par3 <- sim.param.list$cov.par3[i]
		sim.cluster.list[[i]]$covariance.params$cov.par4 <- sim.param.list$cov.par4[i]
		sim.cluster.list[[i]]$cluster.mean <- sim.param.list$cluster.mean[i]
		sim.cluster.list[[i]]$admix.prop.matrix <- sim.param.list$admix.props[,i] %*% t(sim.param.list$admix.props[,i])
		sim.cluster.list[[i]]$covariance <- cluster.covariance(sampling.data$geo.dist,
																sampling.data$time.dist,
																sim.cluster.list[[i]]$covariance.params)
	}
	return(sim.cluster.list)
}

simulate.freq.data <- function(n.loci,n.samples,sim.admixed.cov.mat){
	freq.data <- t(MASS::mvrnorm(n=n.loci,mu=numeric(n.samples),Sigma=sim.admixed.cov.mat))
	return(freq.data)
}

simulate.spatialStructure.dataset <- function(n.clusters,n.samples,n.loci,sample.sizes=NULL,file.name,time.sampling=FALSE,sim.seed=NULL,counts=FALSE){
	# recover()
# set seed
	if(is.null(sim.seed)){
		sim.seed <- sample(1:1e5,1)
	}
	cat(sim.seed)
	set.seed(sim.seed)
# simulate locations
	sampling.data <- list("geo.coords" = cbind(runif(n.samples,0,50),runif(n.samples,0,50)),
							"temporal.coords" = ifelse(rep(time.sampling,n.samples),
													c(sample(1:500, n.samples,replace=TRUE)),
													0))
	sampling.data$geo.dist <- fields::rdist(sampling.data$geo.coords)
	sampling.data$time.dist <- fields::rdist(sampling.data$temporal.coords)
	sim.param.list <- generate.sim.param.list(n.clusters,n.samples,sampling.data,time.sampling)
	sim.cluster.list <- lapply(seq_along(1:n.clusters),function(i){declare.cluster.list()})
	sim.cluster.list <- populate.sim.cluster.list(sim.cluster.list,sim.param.list,sampling.data)
	sim.admixed.cov.mat <- admixed.covariance(sim.cluster.list,n.clusters,sim.param.list$shared.mean,sim.param.list$nuggets,rep(0,n.samples),0)
	sim.data <- list("n.samples" = n.samples,"n.loci" = n.loci,"sampling.data" = sampling.data,
							"sim.param.list" = sim.param.list,"sim.cluster.list" = sim.cluster.list,
							"sim.admixed.cov.mat" = sim.admixed.cov.mat,"sim.seed"=sim.seed)
	if(counts){
		if(length(sample.sizes) < n.samples){
			sim.data$sample.sizes <- rep(sample.sizes,n.samples)
		}
		sim.admixed.cov.mat <- sim.admixed.cov.mat - sim.param.list$shared.mean
		ancestral.freqs <- rnorm(n.loci,mean=0.5,sd=sqrt(sim.param.list$shared.mean))
		sim.data$sim.param.list$shared.mean <- var(ancestral.freqs)
		sim.freqs <- simulate.freq.data(n.loci,n.samples,sim.admixed.cov.mat)
		sim.freqs <- sim.freqs + matrix(ancestral.freqs,nrow=n.samples,ncol=n.loci,byrow=TRUE)
		sim.freqs[which(sim.freqs < 0 )] <- 0
		sim.freqs[which(sim.freqs > 1 )] <- 1
		sim.data$sim.admixed.cov.mat <- sim.admixed.cov.mat + sim.data$sim.param.list$shared.mean
		sim.data$ancestral.freqs <- ancestral.freqs
		sim.data$sim.freqs <- switcharoo.data(sim.freqs)
		sim.data$sim.counts <- matrix(rbinom(n.samples*n.loci,sample.sizes,sim.freqs),nrow=n.samples,ncol=n.loci)
		sim.data$sim.sample.freqs <- sim.data$sim.counts/matrix(sample.sizes,nrow=n.samples,ncol=n.loci)
		sim.data$sample.covariance <- cov(t(sim.data$sim.sample.freqs))
	}
	save(sim.data,file=paste(file.name,".Robj",sep=""))
}

# par(mfrow=c(5,5),mar=c(1,1,1,1),oma=c(1,1,1,1))

simulate.spatialStructure.dataset(n.clusters = 1,
									n.samples = 30,
									n.loci = 10000,
									sample.sizes = 100,
									file.name = "~/desktop/k1_exp_dataset", #"~/desktop/test.Robj"
									counts = TRUE,
									time.sampling = FALSE,
									sim.seed=NULL)
# load("~/desktop/test.Robj")
load("~/desktop/k1_exp_dataset.Robj")

quartz(width=10,height=5)
par(mfrow=c(2,3))
y.lim <- range(sim.data$sample.covariance,
				range(unlist(lapply(sim.data$sim.cluster.list,"[[","covariance"))) + 
				range(unlist(lapply(sim.data$sim.cluster.list,"[[","cluster.mean"))))
plot(sim.data$sampling$geo.dist,sim.data$sample.covariance,ylim=y.lim)
	cluster.cols <- c("blue","red","green")
	lapply(1:length(sim.data$sim.cluster.list),function(i){
				points(sim.data$sampling$geo.dist,
						sim.data$sim.cluster.list[[i]]$covariance + 
						sim.data$sim.cluster.list[[i]]$cluster.mean,
						col=cluster.cols[i],pch=20,cex=0.5);
				abline(h=sim.data$sim.cluster.list[[i]]$cluster.mean,col=cluster.cols[i]) ; 
				return(invisible(0))})
	abline(h=sim.data$sim.param.list$shared.mean,col="gray")

plot(sim.data$sampling$time.dist,sim.data$sample.covariance,ylim=y.lim)
	cluster.cols <- c("blue","red","green")
	lapply(1:length(sim.data$sim.cluster.list),function(i){
				points(sim.data$sampling$time.dist,
						sim.data$sim.cluster.list[[i]]$covariance + 
						sim.data$sim.cluster.list[[i]]$cluster.mean,
						col=cluster.cols[i],pch=20,cex=0.5);
				abline(h=sim.data$sim.cluster.list[[i]]$cluster.mean,col=cluster.cols[i]) ;
				return(invisible(0))})
	abline(h=sim.data$sim.param.list$shared.mean,col="gray")
	
# quartz()
# par(mfrow=c(1,3))
plot(sim.data$sim.admixed.cov.mat+var(sim.data$ancestral.freqs),sim.data$sample.covariance,xaxt='n',yaxt='n')
mtext(side=1,text=sim.data$sim.seed,padj=-2)
	abline(0,1,col="red")


plot(sim.data$sampling.data$geo.dist,sim.data$sim.admixed.cov.mat,
		ylim=c(0,max(sim.data$sim.admixed.cov.mat)))
	points(sim.data$sampling.data$geo.dist,sim.data$sample.covariance,col="red",pch=20)
plot(sim.data$sampling.data$time.dist,sim.data$sim.admixed.cov.mat,
		ylim=c(0,max(sim.data$sim.admixed.cov.mat)))
	points(sim.data$sampling.data$time.dist,sim.data$sample.covariance,col="red",pch=20)


if(FALSE){
	91944
# tmp.cov <- spatiotemporal.covariance(sim.data$sampling.data$geo.dist,sim.data$sampling.data$time.dist,0.1,50,1,Inf)
# plot(sim.data$sampling.data$geo.dist,tmp.cov)

$Cluster_1
List of 4
 $ covariance.params:List of 3
  ..$ cov.par1: num 0.0249
  ..$ cov.par2: num 0.606
  ..$ cov.par3: num 1831124
 $ cluster.mean     : num 0.00283

$Cluster_2
List of 4
 $ covariance.params:List of 3
  ..$ cov.par1: num 0.0391
  ..$ cov.par2: num 1.51
  ..$ cov.par3: num 904281
 $ cluster.mean     : num 0.0187

$Cluster_3
List of 4
 $ covariance.params:List of 3
  ..$ cov.par1: num 0.0771
  ..$ cov.par2: num 9.38
  ..$ cov.par3: num 1273751
 $ cluster.mean     : num 0.475


$Cluster_1
$Cluster_1$cov.par1
[1] 0.3545136

$Cluster_1$cov.par2
[1] 2.413376

$Cluster_1$cov.par3
[1] 0


$Cluster_2
$Cluster_2$cov.par1
[1] 0.02812039

$Cluster_2$cov.par2
[1] 1.7381

$Cluster_2$cov.par3
[1] 0


$Cluster_3
$Cluster_3$cov.par1
[1] 0.4119261

$Cluster_3$cov.par2
[1] 17.90571

$Cluster_3$cov.par3
[1] 0




}


88097