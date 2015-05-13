spatiotemporal.covariance <- function(geo.dist,time.dist,alpha,nu,zeta){
    d <- 2  # two spatial dimensions
    return( exp(
            d/2 * log(pi) +
            d * log(alpha) -
            ( (nu + zeta*time.dist - 1)*log(2) +
                    gsl::lngamma(nu + zeta*time.dist + d/2) 
                ) +
            nu * log(alpha*geo.dist) +
            gsl::bessel_lnKnu(nu=nu,x=alpha*geo.dist)
    ) )
}

cluster.covariance <- function(geo.dist,time.dist,params){
	spatiotemporal.covariance(geo.dist,
								time.dist,
								params$cov.par1,
								params$cov.par2,
								params$cov.par3)
}

cluster.admixed.covariance <- function(cluster){
	cluster$admix.prop.matrix * (cluster$covariance + cluster$cluster.mean)
}

admixed.covariance <- function(cluster.list,n.clusters,shared.mean,nuggets){
	Reduce("+",lapply(cluster.list,FUN=cluster.admixed.covariance)) + shared.mean + diag(nuggets)
}

calculate.likelihood.2 <- function(data, covar.inverse, covar.det){
	lnL <- -0.5 * data$n.loci * (sum( covar.inverse * data$sample.covariance ) +  covar.det)
	return(lnL)
}

declare.super.list <- function(){
	super.list <- list(
		"parameter.list",
		"output.list",
		"mcmc.quantities")
	class(super.list) <- "super.list"
	return(super.list)
}

print.super.list <- function(super.list){
	print(str(super.list,max.level=1))
}

initialize.super.list <- function(data.list,model.options,mcmc.options,initial.parameters=NULL){
	super.list <- declare.super.list()
	super.list$parameter.list <- initialize.param.list(data.list,model.options,initial.parameters)
	super.list$output.list <- make.output.list(data.list$n.ind,model.options,mcmc.options)
	super.list$mcmc.quantities <- initialize.mcmc.quantities(data.list,super.list$parameter.list,model.options,mcmc.options,initial.parameters)
	super.list$model.options <- model.options
	super.list <- bookkeep(super.list,1,mcmc.options)
	return(super.list)
}

declare.cluster.list <- function(){
	cluster.list <- list(
		"covariance.params" = list(	"cov.par1" = NULL,
									"cov.par2" = NULL,
									"cov.par3" = NULL),
		"cluster.mean" = NULL,
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

make.admix.prop.matrix <- function(admix.proportions){
	admix.proportions%*%t(admix.proportions)
}

populate.cluster <- function(cluster,geo.dist,time.dist,admix.proportions,model.options){
	# recover()
	if(!model.options$no.st){
		cluster$covariance <- NaN
		while(any(is.na(cluster$covariance))){
			cluster$covariance.params$cov.par1 <- runif(1,1e-10,10)
			cluster$covariance.params$cov.par2 <- runif(1,1e-10,10)
			cluster$covariance.params$cov.par3 <- runif(1,1e-10,10)
			cluster$covariance <- cluster.covariance(geo.dist,time.dist,cluster$covariance.params)
		}
	} else {
			cluster$covariance.params$cov.par1 <- 0
			cluster$covariance.params$cov.par2 <- 0
			cluster$covariance.params$cov.par3 <- 0
			cluster$covariance <- matrix(0,nrow=nrow(geo.dist),ncol=ncol(geo.dist))
	}
	cluster$cluster.mean <- rexp(1)
	cluster$admix.prop.matrix <- make.admix.prop.matrix(admix.proportions)
	return(cluster)
}

populate.cluster.list <- function(cluster.list,data,admix.proportions,model.options){
	# recover()
	setNames(lapply(seq_along(cluster.list),
						FUN=function(i){
								populate.cluster(cluster.list[[i]],
												data$geo.dist,
												data$time.dist,
												admix.proportions[,i,drop=FALSE],
												model.options)}),
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
		parameters$shared.mean <- initial.parameters$shared.mean
		parameters$admix.proportions <- initial.parameters$admix.proportions
		parameters$nuggets <- initial.parameters$nuggets
		parameters$cluster.list <- initial.parameters$cluster.list
		parameters$admixed.covariance <- admixed.covariance(parameters$cluster.list,model.options$n.clusters,parameters$shared.mean,parameters$nuggets)
		parameters$inverse <- solve(parameters$admixed.covariance)
		parameters$determinant <- determinant(parameters$admixed.covariance,logarithm=TRUE)$modulus
	} else {
		parameters$shared.mean <- min(data$sample.covariance)
		parameters$admix.proportions <- gtools::rdirichlet(n = n.ind,alpha = rep(1,model.options$n.clusters))
		parameters$nuggets <- rexp(n.ind)
		parameters$cluster.list <- populate.cluster.list(parameters$cluster.list,data,parameters$admix.proportions,model.options)
		parameters$admixed.covariance <- admixed.covariance(parameters$cluster.list,model.options$n.clusters,parameters$shared.mean,parameters$nuggets)
		parameters$inverse <- solve(parameters$admixed.covariance)
		parameters$determinant <- determinant(parameters$admixed.covariance,logarithm=TRUE)$modulus
	}
	class(parameters) <- "parameters"
	return(parameters)
}

print.parameters <- function(parameters){
	print(str(parameters))
}

make.output.list <- function(n.ind,model.options,mcmc.options){
	output.length <- 1 + mcmc.options$ngen/mcmc.options$samplefreq
	output.list <- list("step" = 1,
						"likelihood" = rep(NA,output.length),
						"posterior.prob" = rep(NA,output.length),
						"nuggets" = matrix(NA,nrow=n.ind,ncol=output.length),
						"admix.proportions" = lapply(1:output.length,"==",NA),
						"shared.mean" = rep(NA,output.length),
						"cluster.params" = list("cov.par1" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
												"cov.par2" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
												"cov.par3" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
												"cluster.mean" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length)),
						"acceptance.rates" = list("cov.par1" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"cov.par2" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"cov.par3" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"cluster.mean" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"nuggets" = matrix(NA,nrow=n.ind,ncol=output.length),
													"admix.proportions" = matrix(NA,nrow=n.ind,ncol=output.length),
													"shared.mean" = rep(NA,output.length)))
	class(output.list) <- "spacetimemix.output"
	return(output.list)
}

print.spacetimemix.output <- function(spacetimemix.output){
	print(str(spacetimemix.output,max.level=1))
}

make.mcmc.quantities <- function(n.ind,model.options,mcmc.options){
	mcmc.quantities <- list("gen" = 1,
							"likelihood" = NULL,
							"posterior.prob" = NULL,
							"dirich.conc.params" = NULL,
							"prior.probs" = list("cov.par1" = rep(NA,model.options$n.clusters),
													"cov.par2" = rep(NA,model.options$n.clusters),
													"cov.par3" = rep(NA,model.options$n.clusters),
													"cluster.mean" = rep(NA,model.options$n.clusters),
													"nuggets" = rep(NA,n.ind),
													"admix.proportions" = matrix(NA,nrow=n.ind,ncol=model.options$n.clusters),
													"shared.mean" = NA),
							"acceptance.rates" = list("cov.par1" = rep(0,model.options$n.clusters),
															 "cov.par2" = rep(0,model.options$n.clusters),
 															 "cov.par3" = rep(0,model.options$n.clusters),
  															 "cluster.mean" = rep(0,model.options$n.clusters),
															 "nuggets" = rep(0,n.ind),
															 "admix.proportions" = rep(0,n.ind),
															 "shared.mean" = 0),
							"adaptive.mcmc" = list( "n.batch" = 0,
													"n.moves" = 
														list("cov.par1" = rep(0,model.options$n.clusters),
															 "cov.par2" = rep(0,model.options$n.clusters),
 															 "cov.par3" = rep(0,model.options$n.clusters),
  															 "cluster.mean" = rep(0,model.options$n.clusters),
															 "nuggets" = rep(0,n.ind),
															 "admix.proportions" = rep(0,n.ind),
															 "shared.mean" = 0),
													"log.stps" = 
														list("cov.par1" = rep(0,model.options$n.clusters),
															 "cov.par2" = rep(0,model.options$n.clusters),
 															 "cov.par3" = rep(0,model.options$n.clusters),
  															 "cluster.mean" = rep(0,model.options$n.clusters),
															 "nuggets" = rep(0,n.ind),
															 "admix.proportions" = rep(0,n.ind),
															 "shared.mean" = 0),
													"nth.batch.accept.rates" = 
														list("cov.par1" = rep(0,model.options$n.clusters),
															 "cov.par2" = rep(0,model.options$n.clusters),
 															 "cov.par3" = rep(0,model.options$n.clusters),
  															 "cluster.mean" = rep(0,model.options$n.clusters),
															 "nuggets" = rep(0,n.ind),
															 "admix.proportions" = rep(0,n.ind),
															 "shared.mean" = 0)),
							"covariance.params.list" = declare.covariance.params.list(model.options),
							"smw.numinst.ticker" = 0)
	class(mcmc.quantities) <- "mcmc.quantities"
	return(mcmc.quantities)
}

print.mcmc.quantities <- function(mcmc.quantities){
	print(str(mcmc.quantities))
}

declare.covariance.params.list <- function(model.options){
	if(model.options$temporal.sampling==FALSE){
		covariance.params.list <- list(
					"covariance.params" = list(
						"cov.par1",
						"cov.par2"),
					"covariance.param.prior.functions" = list(
						prior.prob.cov.par1,
						prior.prob.cov.par2))
	} else {
		covariance.params.list <- list(
					"covariance.params" = list(
						"cov.par1",
						"cov.par2",
						"cov.par3"),
					"covariance.param.prior.functions" = list(
						prior.prob.cov.par1,
						prior.prob.cov.par2,
						prior.prob.cov.par3))
	}
	return(covariance.params.list)
}

prior.prob.param.list <- function(parameter.list,element,function.name,other.args=NULL){
	if(is.null(other.args)){
		function.name(parameter.list[[element]])
	} else {
		function.name(parameter.list[[element]],other.args)	
	}
}

prior.prob.nuggets <- function(nuggets){
	dexp(nuggets,log=TRUE)
}

prior.prob.admix.proportions <- function(x,alpha){
	rowSums((matrix(alpha,nrow=nrow(x),ncol=ncol(x),byrow=TRUE) - 1) * log(x))
}

prior.prob.shared.mean <- function(shared.mean){
	dexp(shared.mean,log=TRUE)
}

prior.prob.covariance.param <- function(cluster,element,function.name){
	function.name(cluster$covariance.params[[element]])
}

prior.prob.covariance.param.clusters <- function(cluster.list,element,function.name){
	unlist(lapply(cluster.list,prior.prob.covariance.param,element,function.name))
}

prior.prob.cov.par1 <- function(cov.par1){
	dunif(cov.par1,0,1e20,log=TRUE)
}

prior.prob.cov.par2 <- function(cov.par2){
	dunif(cov.par2,0,1e20,log=TRUE)
}

prior.prob.cov.par3 <- function(cov.par3){
	dunif(cov.par3,0,1e20,log=TRUE)
}

prior.prob.cluster.mean <- function(cluster.mean){
	dexp(cluster.mean,1,log=TRUE)
}

initialize.mcmc.quantities <- function(data.list,parameter.list,model.options,mcmc.options,initial.parameters=NULL){
	# recover()
	mcmc.quantities <- make.mcmc.quantities(data.list$n.ind,model.options,mcmc.options)
	# if(!is.null(initial.parameters)){	#WILL HAVE TO FIX EVENTUALLY
		# stop("not built yet")
	# } else {
		mcmc.quantities$likelihood <- calculate.likelihood.2(data.list,parameter.list$inverse,parameter.list$determinant)
		mcmc.quantities$prior.probs$nuggets <- prior.prob.param.list(parameter.list,"nuggets",prior.prob.nuggets)
		mcmc.quantities$dirich.conc.params <- matrix(rep(0.5,model.options$n.clusters),nrow=data.list$n.ind,ncol=model.options$n.clusters,byrow=TRUE)
		mcmc.quantities$prior.probs$admix.proportions <- prior.prob.param.list(parameter.list,"admix.proportions",prior.prob.admix.proportions,other.args=mcmc.quantities$dirich.conc.params)
		mcmc.quantities$prior.probs$shared.mean <- prior.prob.param.list(parameter.list,"shared.mean",prior.prob.shared.mean)
		if(!model.options$no.st){
			mcmc.quantities$prior.probs$cov.par1 <- prior.prob.covariance.param.clusters(parameter.list$cluster.list,"cov.par1",prior.prob.cov.par1)
			mcmc.quantities$prior.probs$cov.par2 <- prior.prob.covariance.param.clusters(parameter.list$cluster.list,"cov.par2",prior.prob.cov.par2)
			if(model.options$temporal.sampling){
				mcmc.quantities$prior.probs$cov.par3 <- prior.prob.covariance.param.clusters(parameter.list$cluster.list,"cov.par3",prior.prob.cov.par3)
			} else {
				mcmc.quantities$prior.probs$cov.par3 <- 0
			}
		}
		mcmc.quantities$prior.probs$cluster.mean <- prior.prob.cluster.mean(unlist(lapply(parameter.list$cluster.list,"[[","cluster.mean")))
		mcmc.quantities$posterior.prob <- mcmc.quantities$likelihood + sum(unlist(mcmc.quantities$prior.probs),na.rm=TRUE)
	# }
	return(mcmc.quantities)
}

##Graham's improved sherman tank, with added panzer busting power
sherman_r <- function(Ap, u, v, u.i=FALSE, v.i=FALSE) {
	if(u.i & v.i){
		determ.update <- 1 + v[v.i] * Ap[v.i,u.i] * u[u.i] 
		return(list(new.inverse = Ap - ((Ap[,u.i]*u[u.i]) %*% t(v[v.i] * Ap[,v.i]))/determ.update, determ.update=determ.update))
  	}
	if(u.i & !v.i){
  		determ.update <- drop(1 + t(v) %*% (Ap[,u.i]*u[u.i]))  	
		return(list(new.inverse = Ap - ((Ap[,u.i]*u[u.i]) %*% t(v) %*% Ap)/determ.update, determ.update=determ.update))
  	}
	if(!u.i & v.i){
		determ.update <- drop(1 + v[v.i] * Ap[,v.i] %*% u)  	
		return(list(new.inverse = Ap - (Ap %*% u %*% (v[v.i] * Ap[,v.i]))/determ.update, determ.update=determ.update))
  	}
	if(!u.i & !v.i){
  		determ.update <- drop(1 + t(v) %*% Ap %*% u)
  		return(list(new.inverse = Ap - (Ap %*% u %*% t(v) %*% Ap)/determ.update, determ.update=determ.update))
  	}
	stop("none of the options in sherman_r() have been specified")
}

update.cluster.mean.k <- function(data.list,super.list){
#	recover()
	k <- sample(1:length(super.list$parameter.list$cluster.list),1)
	accepted.move <- 0
	super.list$mcmc.quantities$adaptive.mcmc$n.moves$cluster.mean[k] <- super.list$mcmc.quantities$adaptive.mcmc$n.moves$cluster.mean[k] + 1
	delta.cluster.mean <- rnorm(1,sd=exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps$cluster.mean[k]))
	new.cluster.mean <- delta.cluster.mean + super.list$parameter.list$cluster.list[[k]]$cluster.mean 
	new.prior.prob <- prior.prob.cluster.mean(new.cluster.mean)
	if(is.finite(new.prior.prob) ){
		u <- super.list$parameter.list$admix.proportions[,k] 
		v <- super.list$parameter.list$admix.proportions[,k] 	
		u <- u * delta.cluster.mean
		inverse.updated <- sherman_r(super.list$parameter.list$inverse,u,v)
		new.determinant <-  super.list$parameter.list$determinant + log(abs(inverse.updated$determ.update)) 
		new.inverse <- inverse.updated$new.inverse		
		new.likelihood <- calculate.likelihood.2(data.list,new.inverse,new.determinant)[1]
		likelihood.ratio <- new.likelihood - super.list$mcmc.quantities$likelihood
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$cluster.mean[k]
		accept.ratio <- exp(likelihood.ratio + prior.ratio)
		if(accept.ratio >= runif(1)){
			super.list$mcmc.quantities$prior.probs$cluster.mean[k]  <- new.prior.prob 
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$cluster.list[[k]]$cluster.mean  <- new.cluster.mean
			super.list$parameter.list$determinant <- new.determinant	
			super.list$parameter.list$inverse <- new.inverse 
			accepted.move <- 1
		}
	}
	super.list$mcmc.quantities$acceptance.rates$cluster.mean[k] <- (super.list$mcmc.quantities$acceptance.rates$cluster.mean[k]*(super.list$mcmc.quantities$adaptive.mcmc$n.moves$cluster.mean[k]-1) + accepted.move)/super.list$mcmc.quantities$adaptive.mcmc$n.moves$cluster.mean[k]
	return(super.list)
}

update.nugget.i <- function(data.list,super.list){
	i <- sample(1:data.list$n.ind,1)
	accepted.move <- 0
	super.list$mcmc.quantities$adaptive.mcmc$n.moves$nuggets[i] <- super.list$mcmc.quantities$adaptive.mcmc$n.moves$nuggets[i] + 1
	delta.nugget <- rnorm(1,0,sd=exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps$nuggets[i]))
	new.nugget <- delta.nugget + super.list$parameter.list$nuggets[i] 
	new.prior.prob <- prior.prob.nuggets(new.nugget)
	if(is.finite(new.prior.prob) ){
		u <- rep(0, data.list$n.ind)
		u[i] <- delta.nugget
		v <- rep(0, data.list$n.ind)
		v[i] <-  1	
		inverse.updated <- sherman_r(super.list$parameter.list$inverse, u=u, v=v, u.i=i, v.i=i)
		new.determinant <- super.list$parameter.list$determinant + log(abs(inverse.updated$determ.update)) 
		new.inverse <- inverse.updated$new.inverse
		new.likelihood <- calculate.likelihood.2(data.list, new.inverse ,  new.determinant )[1]
		likelihood.ratio <- new.likelihood - super.list$mcmc.quantities$likelihood
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$nuggets[i]
		accept.ratio <- exp(likelihood.ratio + prior.ratio)
		if(accept.ratio >= runif(1)){
			super.list$mcmc.quantities$prior.probs$nuggets[i]  <- new.prior.prob 
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$nuggets[i]  <- new.nugget
			super.list$parameter.list$determinant <- new.determinant	
			super.list$parameter.list$inverse <- new.inverse
			accepted.move <- 1
		}
	}
	super.list$mcmc.quantities$acceptance.rates$nuggets[i] <- (super.list$mcmc.quantities$acceptance.rates$nuggets[i]*(	super.list$mcmc.quantities$adaptive.mcmc$n.moves$nugget[i]-1) + accepted.move)/	super.list$mcmc.quantities$adaptive.mcmc$n.moves$nugget[i]
	return(super.list)
}

update.shared.mean <- function(data.list,super.list){
	accepted.move <- 0
	super.list$mcmc.quantities$adaptive.mcmc$n.moves$shared.mean <- super.list$mcmc.quantities$adaptive.mcmc$n.moves$shared.mean + 1
	delta.shared.mean <- rnorm(1,sd=exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps$shared.mean)) 
	new.shared.mean <- delta.shared.mean+ super.list$parameter.list$shared.mean
	new.prior.prob <- prior.prob.shared.mean(new.shared.mean)
	if(is.finite(new.prior.prob) ){
		u = rep(delta.shared.mean, data.list$n.ind)
		v = rep(1, data.list$n.ind)
		inverse.updated <- sherman_r(super.list$parameter.list$inverse,u,v)   
		new.determinant <- super.list$parameter.list$determinant + log(abs(inverse.updated$determ.update)) 
		new.inverse <- inverse.updated$new.inverse
		new.likelihood <-   calculate.likelihood.2(data.list, new.inverse ,  new.determinant )[1]
		likelihood.ratio <- new.likelihood - super.list$mcmc.quantities$likelihood
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$shared.mean
		accept.ratio <- exp(likelihood.ratio + prior.ratio)
		if(accept.ratio >= runif(1)){
			super.list$mcmc.quantities$prior.probs$shared.mean  <- new.prior.prob 
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$shared.mean <- new.shared.mean
			super.list$parameter.list$determinant <- new.determinant	
			super.list$parameter.list$inverse <- new.inverse
			accepted.move <- 1
		}
	}
	super.list$mcmc.quantities$acceptance.rates$shared.mean <- (super.list$mcmc.quantities$acceptance.rates$shared.mean*(super.list$mcmc.quantities$adaptive.mcmc$n.moves$shared.mean-1) + accepted.move)/super.list$mcmc.quantities$adaptive.mcmc$n.moves$shared.mean
	return(super.list)
}

update.w.i <- function(data.list,super.list){ 
	# recover()
	i <- sample(1:data.list$n.ind,1)
	accepted.move <- 0
	super.list$mcmc.quantities$adaptive.mcmc$n.moves$admix.proportions[i] <- super.list$mcmc.quantities$adaptive.mcmc$n.moves$admix.proportions[i] + 1
	these.two <- sample(super.list$model.options$n.clusters,2)
	clst.1 <- these.two[1]
	clst.2 <- these.two[2]
	delta.w <- rnorm(1,sd=exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps$admix.proportions[i]))
	new.w.1 <- super.list$parameter.list$admix.proportions[i,clst.1] + delta.w
	new.w.2 <- super.list$parameter.list$admix.proportions[i,clst.2] - delta.w
	new.admixture.vec <- super.list$parameter.list$admix.proportions[i,,drop=FALSE]
	new.admixture.vec[these.two] <- c(new.w.1,new.w.2)
	#write our own dirichlet prior prob that doesn't spit "log(x) NaNs produced" warnings
	new.prior.prob <- prior.prob.admix.proportions(new.admixture.vec,super.list$mcmc.quantities$dirich.conc.params[i,these.two,drop=FALSE])
	if(is.finite(new.prior.prob)){
		covar.1 <- super.list$parameter.list$cluster.list[[clst.1]]$covariance[i,]  + super.list$parameter.list$cluster.list[[clst.1]]$cluster.mean
		covar.2 <- super.list$parameter.list$cluster.list[[clst.2]]$covariance[i,]  + super.list$parameter.list$cluster.list[[clst.2]]$cluster.mean
		u <- delta.w * (super.list$parameter.list$admix.proportions[,clst.1] * covar.1  - super.list$parameter.list$admix.proportions[,clst.2] * covar.2)
		u[i] <- u[i] + delta.w^2 * (covar.1[i] + covar.2[i])/2 	### (wi^k + Dw) (wi^k + Dw) 
		v <- rep(0,data.list$n.ind)
		v[i] <- 1
		matrix.updated.once <- sherman_r(super.list$parameter.list$inverse, u=u, v=v ,v.i=i)  
		matrix.updated.twice <- sherman_r(matrix.updated.once$new.inverse, u=v, v=u, u.i=i) ##note deliberate switcharoo  
		new.determinant <-  super.list$parameter.list$determinant + log(abs(matrix.updated.once$determ.update)) + log(abs(matrix.updated.twice$determ.update))
		new.inverse <- matrix.updated.twice$new.inverse
		new.likelihood <- calculate.likelihood.2(data.list, new.inverse, new.determinant)
		likelihood.ratio <- new.likelihood - super.list$mcmc.quantities$likelihood
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$admix.proportions[i]
		accept.ratio <- exp(likelihood.ratio + prior.ratio)
		if(accept.ratio >= runif(1)){
			super.list$mcmc.quantities$prior.probs$admix.proportions[i] <- new.prior.prob
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$admix.proportions[i,] <- new.admixture.vec
			super.list$parameter.list$determinant <- new.determinant
			super.list$parameter.list$inverse <- new.inverse					
			accepted.move <- 1
		}
	}
	super.list$mcmc.quantities$acceptance.rates$admix.proportions[i] <- (super.list$mcmc.quantities$acceptance.rates$admix.proportions[i]*(super.list$mcmc.quantities$adaptive.mcmc$n.moves$admix.proportion[i]-1) + accepted.move)/super.list$mcmc.quantities$adaptive.mcmc$n.moves$admix.proportion[i]
	return(super.list)
}

propose.cov.param.update <- function(super.list,this.cluster,this.param){
	super.list$parameter.list$cluster.list[[this.cluster]]$covariance.params[[this.param]] + 
		rnorm(1,0,sd=exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps[[this.param]][this.cluster]))
}

recalculate.cluster.list <- function(data.list,cluster.list,this.cluster,this.param,new.param){
	cluster.list[[this.cluster]]$covariance.params[[this.param]] <- new.param
	cluster.list[[this.cluster]]$covariance <- cluster.covariance(data.list$geo.dist,
																			data.list$time.dist,
																			cluster.list[[this.cluster]]$covariance.params)
	return(cluster.list)
}

take.cov.inverse <- function(matrix){
	MASS::ginv(matrix)
	# new.inverse <- tryCatch({
						# MASS::ginv(matrix)
							# },error=function(e){
								# return("skip")
							# })
	# return(new.inverse)
}

update.cluster.covariance.param <- function(data.list,super.list){
	accepted.move <- 0
	this.cluster <- sample(super.list$model.options$n.clusters,1)
	this.param <- sample(length(super.list$mcmc.quantities$covariance.params.list$covariance.params),1)
	super.list$mcmc.quantities$adaptive.mcmc$n.moves[[this.param]][this.cluster] <- super.list$mcmc.quantities$adaptive.mcmc$n.moves[[this.param]][this.cluster] + 1
	this.prior <- super.list$mcmc.quantities$covariance.params.list$covariance.param.prior.functions[[this.param]]
	new.param <- propose.cov.param.update(super.list,this.cluster,this.param)
	new.prior.prob <- this.prior(new.param)
	if(is.finite(new.prior.prob)){
		new.cluster.list <- recalculate.cluster.list(data.list,super.list$parameter.list$cluster.list,this.cluster,this.param,new.param)
		new.cluster.list <- recalculate.cluster.admix.prop.matrix(new.cluster.list,super.list)
		new.admixed.covariance <- admixed.covariance(new.cluster.list,super.list$model.options$n.clusters,super.list$parameter.list$shared.mean,super.list$parameter.list$nuggets)
		new.inverse <- take.cov.inverse(new.admixed.covariance)
		# if(is.numeric(new.inverse)){
			new.determinant <- determinant(new.admixed.covariance,logarithm=TRUE)$modulus
			new.likelihood <- calculate.likelihood.2(data.list, new.inverse, new.determinant)
			likelihood.ratio <- new.likelihood - super.list$mcmc.quantities$likelihood
			prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs[[this.param]][this.cluster]
			accept.ratio <- exp(likelihood.ratio + prior.ratio)
			if(accept.ratio >= runif(1)){
				super.list$parameter.list$cluster.list <- new.cluster.list
				super.list$parameter.list$admixed.covariance <- new.admixed.covariance
				super.list$parameter.list$inverse <- new.inverse
				super.list$parameter.list$determinant <- new.determinant
				super.list$mcmc.quantities$prior.probs[[this.param]][this.cluster] <- new.prior.prob
				super.list$mcmc.quantities$likelihood <- new.likelihood
				super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
				accepted.move <- 1
			}
		# }
	}
	super.list$mcmc.quantities$acceptance.rates[[this.param]][this.cluster] <- (super.list$mcmc.quantities$acceptance.rates[[this.param]][this.cluster]*(super.list$mcmc.quantities$adaptive.mcmc$n.moves[[this.param]][this.cluster]-1) + accepted.move)/super.list$mcmc.quantities$adaptive.mcmc$n.moves[[this.param]][this.cluster]
	return(super.list)
}

index.cluster.param <- function(cluster.list,param){
	unlist(lapply(cluster.list,FUN=function(cluster.list){cluster.list$covariance.params[[param]]}))
}

cluster.make.admix.prop.matrix <- function(cluster,admix.proportions){
	cluster$admix.prop.matrix <- make.admix.prop.matrix(admix.proportions)
	return(cluster)
}

recalculate.cluster.admix.prop.matrix <- function(cluster.list,super.list){
	#Gideon tried to do this with an apply, and it changed the structure of the list
	#	and also, the for-loop was faster.  Gideon wants answers!
	for(i in 1:super.list$model.options$n.clusters){
		cluster.list[[i]]$admix.prop.matrix <- make.admix.prop.matrix(super.list$parameter.list$admix.proportions[,i])
	}														
	return(cluster.list)
}

tidy.super.list <- function(super.list){
	super.list$parameter.list$cluster.list <- recalculate.cluster.admix.prop.matrix(super.list$parameter.list$cluster.list,super.list)
	super.list$parameter.list$admixed.covariance <- admixed.covariance(super.list$parameter.list$cluster.list,super.list$model.options$n.clusters,super.list$parameter.list$shared.mean,super.list$parameter.list$nuggets)	
	return(super.list)
}

recalibrate.super.list <- function(super.list,data.list){
	if(sum(abs(super.list$parameter.list$admixed.covariance %*% super.list$parameter.list$inverse - diag(data.list$n.ind))) > data.list$n.ind * 1e-10){
		super.list$parameter.list$inverse <- MASS::ginv(super.list$parameter.list$admixed.covariance)
		super.list$parameter.list$determinant <- determinant(super.list$parameter.list$admixed.covariance,logarithm=TRUE)$modulus
		super.list$mcmc.quantities$smw.numinst.ticker <- super.list$mcmc.quantities$smw.numinst.ticker + 1
	}
	return(super.list)
}	

bookkeep <- function(super.list,generation,mcmc.options,data.list){
	if(generation %% mcmc.options$samplefreq == 0){
		step <- super.list$output.list$step
		super.list <- bookkeep.params(super.list,step)
		super.list <- bookkeep.acceptance.rates(super.list,step)
		super.list$output.list$step <- step + 1
	}
	super.list <- bookkeep.adaptive.mcmc(super.list,generation)
	if(generation %% mcmc.options$printfreq == 0){
		cat("MCMC generation: ",generation," --- ","Post Prob: ",super.list$mcmc.quantities$posterior.prob,"\n")
	}
	if(generation %% mcmc.options$savefreq == 0){
		super.list <- tidy.super.list(super.list)
		save(super.list,file=mcmc.options$output.file.name)
	}
	if(generation %% 1e5 == 0){
		super.list <- tidy.super.list(super.list)
		super.list <- recalibrate.super.list(super.list,data.list)
	}
	return(super.list)
}

bookkeep.adaptive.mcmc <- function(super.list,generation){
	# recover()
	adaption.step <- ifelse(generation%%50==0,0,1)
	if(!adaption.step){
		super.list <- bookkeep.lstps(super.list)
	}
	super.list$mcmc.quantities <- bookkeep.nth.batch.acceptance.rates(super.list$mcmc.quantities,adaption.step)
	return(super.list)
}

bookkeep.nth.batch.acceptance.rates <- function(mcmc.quantities,adaption.step){
	# recover()
		mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cov.par1 <- adaption.step * (mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cov.par1 + mcmc.quantities$acceptance.rates$cov.par1/50) ;
		mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cov.par2 <- adaption.step * (mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cov.par2 + mcmc.quantities$acceptance.rates$cov.par2/50) ;
		mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cov.par3 <- adaption.step * (mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cov.par3 + mcmc.quantities$acceptance.rates$cov.par3/50) ;
		mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cluster.mean <- adaption.step * (mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$cluster.mean + mcmc.quantities$acceptance.rates$cluster.mean/50) ;
		mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$nuggets <- adaption.step * (mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$nuggets + mcmc.quantities$acceptance.rates$nuggets/50) ;
		mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$admix.proportions <- adaption.step * (mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$admix.proportions + mcmc.quantities$acceptance.rates$admix.proportions/50) ;
		mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$shared.mean <- adaption.step * (mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates$shared.mean + mcmc.quantities$acceptance.rates$shared.mean/50)
	return(mcmc.quantities)
}

update.log.stp <- function(log.stp,n.batch,accept.rate){
	ifelse(accept.rate > 0.44,
			min((log.stp + min(n.batch^(-0.5),0.01)),7),
			min((log.stp - min(n.batch^(-0.5),0.01)),7))
}

bookkeep.lstps <- function(super.list){
	super.list$mcmc.quantities$adaptive.mcmc$log.stps <- mapply(update.log.stp,
																log.stp = super.list$mcmc.quantities$adaptive.mcmc$log.stps,
																n.batch=super.list$mcmc.quantities$adaptive.mcmc$n.batch,
																accept.rate=super.list$mcmc.quantities$adaptive.mcmc$nth.batch.accept.rates)
	super.list$mcmc.quantities$adaptive.mcmc$n.batch <- super.list$mcmc.quantities$adaptive.mcmc$n.batch + 1
	return(super.list)
}

bookkeep.acceptance.rates <- function(super.list,step){
	super.list$output.list$acceptance.rates$shared.mean[step] <- super.list$mcmc.quantities$acceptance.rates$shared.mean
	super.list$output.list$acceptance.rates$admix.proportions[,step] <- super.list$mcmc.quantities$acceptance.rates$admix.proportions
	super.list$output.list$acceptance.rates$nuggets[,step] <- super.list$mcmc.quantities$acceptance.rates$nuggets
	super.list$output.list$acceptance.rates$cov.par1[,step] <- super.list$mcmc.quantities$acceptance.rates$cov.par1
	super.list$output.list$acceptance.rates$cov.par2[,step] <- super.list$mcmc.quantities$acceptance.rates$cov.par2
	super.list$output.list$acceptance.rates$cov.par3[,step] <- super.list$mcmc.quantities$acceptance.rates$cov.par3
	super.list$output.list$acceptance.rates$cluster.mean[,step] <- super.list$mcmc.quantities$acceptance.rates$cluster.mean
	return(super.list)
}

bookkeep.params <- function(super.list,step){
	super.list$output.list$likelihood[step] <- super.list$mcmc.quantities$likelihood
	super.list$output.list$posterior.prob[step] <- super.list$mcmc.quantities$posterior.prob
	super.list$output.list$shared.mean[step] <- super.list$parameter.list$shared.mean
	super.list$output.list$admix.proportions[[step]] <- super.list$parameter.list$admix.proportions
	super.list$output.list$nuggets[,step] <- super.list$parameter.list$nuggets
	super.list$output.list$cluster.params$cov.par1[,step] <- index.cluster.param(super.list$parameter.list$cluster.list,"cov.par1")
	super.list$output.list$cluster.params$cov.par2[,step] <- index.cluster.param(super.list$parameter.list$cluster.list,"cov.par2")
	super.list$output.list$cluster.params$cov.par3[,step] <- index.cluster.param(super.list$parameter.list$cluster.list,"cov.par3")
	super.list$output.list$cluster.params$cluster.mean[,step] <- unlist(lapply(super.list$parameter.list$cluster.list,"[[","cluster.mean"))
	return(super.list)
}

make.prob.intervals <- function(function.sample.probs){
	intervals <- 1000*cumsum(c(0, function.sample.probs))
	intervals[1] <- 1
	return(intervals)
}

make.update.parameters.list <- function(model.options){
	# recover()
	update.function.list <- list("update.w.i " = update.w.i,
							"update.cluster.covariance.param" = update.cluster.covariance.param,
							"update.cluster.mean.k" =  update.cluster.mean.k,
							"update.nugget.i" = update.nugget.i,
							"update.shared.mean" = update.shared.mean)
	function.sample.probs <- c(0.8,0.17,0.01,0.01,0.01)
	if(model.options$no.st){
		drop <- grepl("update.cluster.covariance.param",names(update.function.list))
		update.function.list <- update.function.list[-which(drop)]
		function.sample.probs <- function.sample.probs[-which(drop)]
	}
	if(model.options$n.cluster == 1){
		drop <- grepl("update.w.i",names(update.function.list))
		update.function.list <- update.function.list[-which(drop)]
		function.sample.probs <- function.sample.probs[-which(drop)]
	}
	function.sample.probs <- function.sample.probs/sum(function.sample.probs)
	function.sample.intervals <- make.prob.intervals(function.sample.probs)
	update.parameters.list <- list("update.function.list" = update.function.list,
									"function.sample.intervals" = function.sample.intervals)
	return(update.parameters.list)
}

update.parameters <- function(data.list,super.list,update.parameters.list){
	x <- sample(1:1000,1)
	update <- findInterval(x,update.parameters.list[["function.sample.intervals"]],all.inside=TRUE)
	super.list <- update.parameters.list$update.function.list[[update]](data.list,super.list)
	return(super.list)
}

#data <- geo.coords,time.coords,sample.covariance
#model.options <- round.earth,n.clusters,temporal.sampling,no.st
#mcmc.options <- ngen,samplefreq,printfreq
#initial.parameters <- initial.params, initial.params.lstps

#GIDEON PLAYGROUND
MCMC.gid <- function(	data,
						model.options,
						mcmc.options,
						initial.parameters=NULL,
						seed=NULL){
	# recover()
	on.exit(save(super.list,file=mcmc.options$output.file.name))
	if(is.null(seed)){
		seed <- sample(1:10000,1)
	}
	cat(seed,"\n")
	set.seed(seed)
	data.list <- make.data.list(data,model.options)
	save(data.list,file="data.list.Robj")
	super.list <- initialize.super.list(data.list,model.options,mcmc.options,initial.parameters)
	update.parameters.list <- make.update.parameters.list(model.options)
	for(z in 2:mcmc.options$ngen){
		super.list$mcmc.quantities$gen <- z
		super.list <- update.parameters(data.list,super.list,update.parameters.list)
		super.list <- bookkeep(super.list,z,mcmc.options,data.list)
	}
	save(super.list,file=mcmc.options$output.file.name)
	return("run complete")
}