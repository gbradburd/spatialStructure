#generate initial parameters 
# that are the same as those used to simulate data,
# EXCEPT for the admixture proportions, which are simulated from a dirichlet


init.admix.props <- 	sim.admix.props
init.nuggets <- 		sim.nuggets
init.cov.par1_1 <- 		sim.cov.par1_1
init.cov.par2_1 <- 		sim.cov.par2_1
init.cov.par3_1 <- 		sim.cov.par3_1
#init.cluster.mean.1 <- 	sim.cluster.mean.1
init.cov.par1_2 <- 		sim.cov.par1_2
init.cov.par2_2 <- 		sim.cov.par2_2
init.cov.par3_2 <- 		sim.cov.par3_2
#init.cluster.mean.2 <- 	sim.cluster.mean.2
#init.shared.mean <- 	sim.shared.mean

# init.admix.props <- 	gtools::rdirichlet(n = k,alpha = rep(1,2))
# init.nuggets 	<- 		rexp(k)
# init.cov.par1_1 <- 		runif(1,0,10)
# init.cov.par2_1 <- 		runif(1,0,10)
# init.cov.par3_1 <- 		rexp(1)
 init.cluster.mean.1 <- 	rexp(1)
# init.cov.par1_2 <- 		runif(1,0,10)
# init.cov.par2_2 <- 		runif(1,0,10)
# init.cov.par3_2 <- 		rexp(1)
 init.cluster.mean.2 <- 	rexp(1)
 init.shared.mean <- 	rexp(1)

initial.parameters <- list("shared.mean" = init.shared.mean,
							"admix.proportions" = init.admix.props,
							"nuggets" = init.nuggets,
							"cluster.list" = generate.clusters(2))
#initial parameters cluster 1
	initial.parameters$cluster.list$Cluster_1$covariance.params$cov.par1 <- init.cov.par1_1
	initial.parameters$cluster.list$Cluster_1$covariance.params$cov.par2 <- init.cov.par2_1
	initial.parameters$cluster.list$Cluster_1$covariance.params$cov.par3 <- init.cov.par3_1
	initial.parameters$cluster.list$Cluster_1$cluster.mean <- init.cluster.mean.1
	initial.parameters$cluster.list$Cluster_1$covariance <- matrix(0,ncol=k,nrow=k)	#spatiotemporal.covariance(geo.dist,time.dist,init.cov.par1_1,init.cov.par2_1,init.cov.par3_1)
	initial.parameters$cluster.list$Cluster_1$admix.prop.matrix <- init.admix.props[,1] %*% t(init.admix.props[,1])

#initial parameters cluster 2
	initial.parameters$cluster.list$Cluster_2$covariance.params$cov.par1 <- init.cov.par1_2
	initial.parameters$cluster.list$Cluster_2$covariance.params$cov.par2 <- init.cov.par2_2
	initial.parameters$cluster.list$Cluster_2$covariance.params$cov.par3 <- init.cov.par3_2
	initial.parameters$cluster.list$Cluster_2$cluster.mean <- init.cluster.mean.2
	initial.parameters$cluster.list$Cluster_2$covariance <- matrix(0,ncol=k,nrow=k)#spatiotemporal.covariance(geo.dist,time.dist,init.cov.par1_2,init.cov.par2_2,init.cov.par3_2)
	initial.parameters$cluster.list$Cluster_2$admix.prop.matrix <- init.admix.props[,2] %*% t(init.admix.props[,2])


load("linepops.ms.dataset.Robj")
sim.sample.covariance <- cov(t(spacemix.dataset$allele.counts/spacemix.dataset$sample.sizes))
sim.data <- list("geo.coords" = spacemix.dataset$population.coordinates,
				"time.coords" = matrix(0,nrow=nrow(spacemix.dataset$population.coordinates),ncol=2),
				"sample.covariance" = sim.sample.covariance,
				"n.loci" = ncol(spacemix.dataset$allele.counts))
model.options = list("round.earth" = FALSE,
						"n.clusters" = 1,
						"temporal.sampling"=FALSE,
						no.st=TRUE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name"="K1test_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
