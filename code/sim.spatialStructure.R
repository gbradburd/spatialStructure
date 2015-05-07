#SIMULATE TOY DATA
sim.seed <- sample(1:1e5,1)
cat(sim.seed)
set.seed(sim.seed)
k <- 50
n.loci <- 1e10
spatial.coords <- cbind(runif(k),runif(k))
temporal.coords <- rep(0,k)  # sample(1:100,k,replace=TRUE)
geo.dist <- fields::rdist(spatial.coords)
time.dist <- fields::rdist(temporal.coords)
sim.nuggets <- rexp(k,rate=10)
sim.admix.props <- gtools::rdirichlet(n = k,alpha = rep(1,2))
# sim.admix.props <- sample(0:1,size=k,replace=TRUE)
	# sim.admix.props[which(sim.admix.props==1)] <- 0.99
	# sim.admix.props[which(sim.admix.props==0)] <- 0.01
	# sim.admix.props <- cbind(sim.admix.props,1-sim.admix.props)
sim.cov.par1_1 <- runif(1)
sim.cov.par2_1 <- runif(1,0,10)
sim.cov.par3_1 <- runif(1)
sim.cluster.mean.1 <- rexp(1,rate=10)
sim.cov.par1_2 <- runif(1)
sim.cov.par2_2 <- runif(1,0,10)
sim.cov.par3_2 <- runif(1)
sim.cluster.mean.2 <- rexp(1,rate=10)
sim.shared.mean <- rexp(1,rate=20)
sim.cluster.list <- list("Cluster_1" = list("admix.prop.matrix" = sim.admix.props[,1] %*% t(sim.admix.props[,1]),
											"cluster.mean" = sim.cluster.mean.1,
											"covariance" = matrix(0,ncol=k,nrow=k)),	#spatiotemporal.covariance(geo.dist,time.dist,sim.cov.par1_1,sim.cov.par2_1,sim.cov.par3_1)),
						 "Cluster_2" = list("admix.prop.matrix" = sim.admix.props[,2] %*% t(sim.admix.props[,2]),
											"cluster.mean" = sim.cluster.mean.2,
											"covariance" = matrix(0,ncol=k,nrow=k)))	#spatiotemporal.covariance(geo.dist,time.dist,sim.cov.par1_2,sim.cov.par2_2,sim.cov.par3_2)))
sim.admixed.cov.mat <- admixed.covariance(sim.cluster.list,2,sim.shared.mean,sim.nuggets)
sim.param.list <- list("k" = k,"n.loci" = n.loci,"spatial.coords" = spatial.coords,
						"temporal.coords" = temporal.coords,"geo.dist" = geo.dist,
						"time.dist" = time.dist,"sim.nuggets" = sim.nuggets,
						"sim.admix.props" = sim.admix.props,"sim.cov.par1_1" = sim.cov.par1_1,
						"sim.cov.par2_1" = sim.cov.par2_1,"sim.cov.par3_1" = sim.cov.par3_1,
						"sim.cluster.mean.1" = sim.cluster.mean.1,
						"sim.cov.par1_2" = sim.cov.par1_2,"sim.cov.par2_2" = sim.cov.par2_2,
						"sim.cov.par3_2" = sim.cov.par3_2,"sim.cluster.mean.2" = sim.cluster.mean.2,
						"sim.shared.mean" = sim.shared.mean,"sim.admixed.cov.mat" = sim.admixed.cov.mat,
						"sim.cluster.list" = sim.cluster.list)
save(sim.param.list,file="simdata.Robj")