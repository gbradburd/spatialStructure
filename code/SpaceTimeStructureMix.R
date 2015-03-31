spatial.covariance <- function(geo.dist,time.dist,geo.effect,time.effect,exponent,sill){
	sill*exp(-(geo.dist/geo.effect + time.dist/time.effect)^exponent)
}

cluster.spatial.covariance <- function(geo.dist,time.dist,params){
	spatial.covariance(geo.dist,
						time.dist,
						params$geo.effect,
						params$time.effect,
						params$exponent,
						params$sill)
}

cluster.admixed.covariance <- function(cluster){
	cluster$admix.prop.matrix * cluster$covariance
}

admixed.covariance <- function(cluster.list,n.clusters){
	Reduce("+",lapply(cluster.list,FUN=cluster.admixed.covariance))
}

calculate.likelihood <- function(data,parameters){
	A <- solve(parameters$admixed.covariance)

	lnL <- -0.5 * data$n.loci * sum( parameters$inverse * data$sample.covariance ) - (data$n.loci/2) * parameters$determinant 
#	lnL <- -0.5 * data$n.loci * sum( A * data$sample.covariance ) - (data$n.loci/2)*determinant(parameters$admixed.covariance,logarithm=TRUE)$modulus
	return(lnL)
}

declare.cluster.list <- function(){
	cluster.list <- list(
		"covariance.params" = list(	"geo.effect" = NULL,
									"time.effect" = NULL,
									"exponent" = NULL,
									"sill" = NULL),
		"covariance" = NULL,
		"admix.prop.matrix" = NULL)
	class(cluster.list) <- "cluster"
	return(cluster.list)
}

print.cluster <- function(cluster){
	print(str(cluster))
}

generate.clusters <- function(n.clusters){
	cluster.names <- unlist(lapply(1:n.clusters,function(i){paste("Cluster_",i,sep="")}))
	clusters <- lapply(1:n.clusters,function(i){declare.cluster.list()})
	names(clusters) <- cluster.names
	class(clusters) <- "clusters"
	return(clusters)
}

print.clusters <- function(clusters){
	print(str(clusters))
}

make.data.list <- function(data,model.options){
	if(model.options$round.earth){
		dist.func <- fields::rdist.earth
	} else {
		dist.func <- fields::rdist
	}
	data.list <- list("geo.dist" = dist.func(data$geo.coords),
						"time.dist" = dist.func(data$time.coords),
						"n.ind" = nrow(data$geo.coords),
						"sample.covariance" = data$sample.covariance,
						"n.loci" = data$n.loci)
	class(data.list) <- "data"
	return(data.list)
}

print.data <- function(data){
	print(str(data))
}

populate.cluster <- function(cluster,geo.dist,time.dist,admix.proportions){
	# recover()
	cluster$covariance.params$geo.effect <- runif(1,1e-10,1e2)
	cluster$covariance.params$time.effect <- runif(1,1e-10,1e2)
	cluster$covariance.params$exponent <- runif(1,0.01,2)
	cluster$covariance.params$sill <- rexp(1,rate=1/100)
	cluster$covariance <- cluster.spatial.covariance(geo.dist,time.dist,cluster$covariance.params)
	cluster$admix.prop.matrix <- admix.proportions%*%t(admix.proportions)
	return(cluster)
}

populate.cluster.list <- function(cluster.list,data,admix.proportions){
	# recover()
	setNames(lapply(seq_along(cluster.list),
						FUN=function(i){
								populate.cluster(cluster.list[[i]],
												data$geo.dist,
												data$time.dist,
												admix.proportions[,i,drop=FALSE])}),
						names(cluster.list))
}

initialize.param.list <- function(data,model.options,initial.parameters=NULL){
	# recover()
	n.ind <- nrow(data$geo.dist)
	parameters <- list("admix.proportions" = NULL,
						"nuggets" = NULL,
						"cluster.list" = generate.clusters(model.options$n.clusters),
						"admixed.covariance" = NULL)
	if(!is.null(initial.parameters)){
		stop("this part isn't built yet")
	} else {
		parameters$admix.proportions <- gtools::rdirichlet(n = n.ind,alpha = rep(1,model.options$n.clusters))
		parameters$nuggets <- rexp(n.ind)
		parameters$cluster.list <- populate.cluster.list(parameters$cluster.list,data,parameters$admix.proportions)
		parameters$admixed.covariance <- admixed.covariance(parameters$cluster.list,model.options$n.clusters)
		parameters$admixed.covariance.inverse <- solve(parameters$admixed.covariance)
		parameters$determinant <- determinant(parameters$admixed.covariance,logarithm=TRUE)$modulus
	}
	class(parameters) <- "parameters"
	return(parameters)
}

print.parameters <- function(parameters){
	print(str(parameters))
}

make.output.list <- function(n.ind,model.options,mcmc.options){
	output.list <- list("likelihood" = rep(NA,mcmc.options$samplefreq),
						"posterior.prob" = rep(NA,mcmc.options$samplefreq),
						"nugget" = matrix(NA,nrow=n.ind,ncol=mcmc.options$samplefreq),
						"admix.proportions" = matrix(NA,nrow=n.ind,ncol=mcmc.options$samplefreq),
						"covariance.params" = list("geo.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq),
													"time.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq),
													"sill" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq),
													"exponent" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq)),
						"acceptance.rates" = list("geo.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq),
													"time.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq),
													"sill" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq),
													"exponent" = matrix(NA,nrow=model.options$n.clusters,ncol=mcmc.options$samplefreq),
													"nugget" = matrix(NA,nrow=n.ind,ncol=mcmc.options$samplefreq),
													"admix.proportions" = matrix(NA,nrow=n.ind,ncol=mcmc.options$samplefreq)))
	class(output.list) <- "spacetimemix.output"
	return(output.list)
}

print.spacetimemix.output <- function(spacetimemix.output){
	print(str(spacetimemix.output))
}

make.mcmc.quantities <- function(n.ind,model.options,mcmc.options){
	mcmc.quantities <- list("likelihood" = NULL,
							"posterior.prob" = NULL,
							"prior.probs" = list("geo.effect" = rep(NA,model.options$n.clusters),
													"time.effect" = rep(NA,model.options$n.clusters),
													"sill" = rep(NA,model.options$n.clusters),
													"exponent" = rep(NA,model.options$n.clusters),
													"nugget" = rep(NA,n.ind),
													"admix.proportions" = matrix(NA,nrow=n.ind,ncol=mcmc.options$samplefreq)),
							"adaptive.mcmc" = list("log.stps" = 
														list("geo.effect" = rep(NA,model.options$n.clusters),
															 "time.effect" = rep(NA,model.options$n.clusters),
 															 "sill" = rep(NA,model.options$n.clusters),
  															 "exponent" = rep(NA,model.options$n.clusters),
															 "nugget" = rep(NA,n.ind),
															 "admix.proportions" = matrix(NA,nrow=n.ind,ncol=mcmc.options$samplefreq)),
													"nth.batch.accept.rates" = 
														list("geo.effect" = rep(NA,model.options$n.clusters),
															 "time.effect" = rep(NA,model.options$n.clusters),
 															 "sill" = rep(NA,model.options$n.clusters),
  															 "exponent" = rep(NA,model.options$n.clusters),
															 "nugget" = rep(NA,n.ind),
															 "admix.proportions" = matrix(NA,nrow=n.ind,ncol=mcmc.options$samplefreq))))
	class(mcmc.quantities) <- "spacetimemix.mcmc.quantities"
	return(mcmc.quantities)
}

≈ <- function(mcmc.quantities){
	print(str(mcmc.quantities))
}

prior.prob.covariance.param <- function(cluster,element,function.name){
	function.name(cluster$covariance.params[[element]])
}

initialize.mcmc.quantities <- function(data.list,parameter.list,model.options,mcmc.options,initial.parameters=NULL){
	mcmc.quantities <- make.mcmc.quantities(data.list$n.ind,model.options,mcmc.options)
	if(!is.null(initial.parameters)){
		stop("not built yet")
	} else {
		mcmc.quantities$likelihood <- calculate.likelihood(data.list,parameter.list)
#		mcmc.quantities$prior.probs <- calculate.prior.probabilities(parameter.list)
	}
	return(mcmc.quantities)
}

##update the w for the ith individual
update.w.i<-function(i, data.list, parameter.list, model.options, mcmc.quantities){
	
	num.clusters<-	model.options$n.clusters
	these.two<-sample(num.clusters,2)
	clst.1<-these.two[1]
	clst.2<-these.two[2]
	
	old.w.1<-parameter.list$admix.proportions[i,clst.1]
	old.w.2<-parameter.list$admix.proportions[i,clst.2]

	delta.w<-rnorm(1,sd=BLAH)
	new.w.1<- old.w.1 + delta.w
	new.w.2<- old.w.2 - delta.w
	
	if( !(new.w.1<0 | new.w.1>1)  &  !(new.w.2<0 | new.w.2>1) )  {
			
		covar.1 <- parameter.list$cluster.list[[clst.1]]$covariance[i,]  #will eventually need the cluster mean
		covar.2 <- parameter.list$cluster.list[[clst.2]]$covariance[i,]  #will eventually need the cluster mean
		
		u <- delta.w * ( parameter.list$admix.proportions[i,clst.1] * covar.1  - parameter.list$admix.proportions[i,clst.2] * covar.2  )
		u[i] <- u[i]  + delta.w^2 * (covar.1[i,i] + covar.2[i,i])/2 	### (wi^k + Dw) (wi^k + Dw) 
		
		v <- rep(0,data.list$n.ind)	
		v[i] <- 1
		
		parameter.list$determinant  #log of determinant
		old.inverse <-pamameter.list$inverse
		old.determinant<-parameter.list$determinant
		matrix.updated<-double.sherman_r(Ap, u, v,i )   #this function is current in Sherman in sandbox
		parameter.list$determinant<-  old.determinant + matrix.updated$determinant.correction   ##update on log-scale
		pamameter.list$inverse  <- matrix.updated$inverse
		new.likelihood <- calculate.likelihood (data.list,parameters.list)  
		new.likelihood - mcmc.quantities$likelihood
	}
	return(parameter.list)
}


#data <- geo.coords,time.coords,sample.covariance
#model.options <- round.earth,n.clusters
#mcmc.options <- ngen,samplefreq,printfreq
#initial.parameters <- initial.params, initial.params.lstps

MCMC <- function(	data,
					model.options,
					mcmc.options,
					initial.parameters=NULL){
	recover()
	data.list <- make.data.list(data,model.options)
	parameter.list <- initialize.param.list(data.list,model.options,initial.parameters)
	output.list <- make.output.list(data.list$n.ind,model.options,mcmc.options)
	mcmc.quantities <- initialize.mcmc.quantities(data.list,parameter.list,model.options,mcmc.options,initial.parameters)
	
	save(parameter.list,mcmc.quantities,file="~/desktop/testobj.Robj")
}

tmpfunc <- function(list1,list2){priorprob <- dexp(list1$geo.effect,log=TRUE) ; }

k <- 100
n.loci <- 1e4
spatial.coords <- cbind(runif(k),runif(k))
temporal.coords <- sample(1:100,k,replace=TRUE)
geo.dist <- fields::rdist(spatial.coords)
time.dist <- fields::rdist(temporal.coords)

par.cov <- spatial.covariance(geo.dist,time.dist,2,3,1,2)

sim.data <- list("geo.coords" = spatial.coords,
			"time.coords" = temporal.coords,
			"sample.covariance" = par.cov,
			"n.loci" = n.loci)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 2)
mcmc.options = list("ngen" = 100,
					"samplefreq" = 10,
					"printfreq" = 5)
MCMC(sim.data,model.options,mcmc.options)












if(FALSE){
geo.dist <- fields::rdist(spatial.coords)
time.dist <- fields::rdist(temporal.coords)

par.cov <- spatial.covariance(geo.dist,time.dist,2,3,1,2)

require(rgl)
plot3d(x = spatial.coords[,1],y = spatial.coords[,2],z = temporal.coords)
plot3d(geo.dist,time.dist,par.cov)

plot.dirichlet <- function(n.ind,n.clust,alphas=NULL){
	# recover()
	if(is.null(alphas)){
		alphas <- runif(n.clust)
	}
	sim.data <- rdirichlet(n.ind,alphas)
	plot(0,xlim=c(0,n.ind),ylim=c(0,1),type='n')	#,main=unlist(lapply(round(alphas,2),paste,","))
		for(i in 1:n.clust){
			points(sim.data[,i],pch=19,col=i)
			text(x=i*(n.ind/(n.clust*2)),y=mean(sim.data[,i]),col=i,labels=round(alphas[i],2))
		}
}
plot.dirichlet(100,3,alphas=c(0.100,0.100,0.100))

}