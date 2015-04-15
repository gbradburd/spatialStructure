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

admixed.covariance <- function(cluster.list,n.clusters,shared.mean){
	Reduce("+",lapply(cluster.list,FUN=cluster.admixed.covariance)) + shared.mean   
	
}

update.admixed.covariance <- function(cluster.list,n.clusters,shared.mean,updated.clusters,admix.prop.matrix.primes){
	cluster.list[[updated.clusters[1]]]$admix.prop.matrix <- admix.prop.matrix.primes[[1]]
	cluster.list[[updated.clusters[2]]]$admix.prop.matrix <- admix.prop.matrix.primes[[2]]
	admixed.covariance.prime <- admixed.covariance(cluster.list,n.cluster,shared.mean)
	return(admixed.covariance.prime)
}

calculate.likelihood <- function(data,parameters){
	#A <- solve(parameters$admixed.covariance)

	lnL <- -0.5 * data$n.loci * sum( parameters$admixed.covariance.inverse * data$sample.covariance ) - (data$n.loci/2) * parameters$determinant 

	lnL <- -0.5 * data$n.loci * sum( parameters$inverse * data$sample.covariance ) - (data$n.loci/2) * parameters$determinant 
	return(lnL)
}

##GRAHAM REWROTE THIS TO TAKE COVAR & DET as functions
calculate.likelihood.2 <- function(data, covar.inverse, covar.det){

	lnL <- -0.5 * data$n.loci * (sum( covar.inverse * data$sample.covariance ) +  covar.det)
#	lnL <- -0.5 * data$n.loci * sum( A * data$sample.covariance ) - (data$n.loci/2)*determinant(parameters$admixed.covariance,logarithm=TRUE)$modulus
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
	print(str(super.list))
}

initialize.super.list <- function(data.list,model.options,mcmc.options,initial.parameters=NULL){
	super.list <- declare.super.list()
	super.list$parameter.list <- initialize.param.list(data.list,model.options,initial.parameters)
	super.list$output.list <- make.output.list(data.list$n.ind,model.options,mcmc.options)
	super.list$mcmc.quantities <- initialize.mcmc.quantities(data.list,super.list$parameter.list,model.options,mcmc.options,initial.parameters)
	super.list$model.options <- model.options
	return(super.list)
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
		parameters$shared.mean <- initial.parameters$shared.mean
		parameters$admix.proportions <- initial.parameters$admix.proportions
		parameters$nuggets <- initial.parameters$nuggets
		parameters$cluster.list <- initial.parameters$cluster.list
		parameters$admixed.covariance <- admixed.covariance(parameters$cluster.list,model.options$n.clusters,parameters$shared.mean)
		parameters$inverse <- solve(parameters$admixed.covariance)
		parameters$determinant <- determinant(parameters$admixed.covariance,logarithm=TRUE)$modulus
	} else {
		parameters$shared.mean <- min(data$sample.covariance)
		parameters$admix.proportions <- gtools::rdirichlet(n = n.ind,alpha = rep(1,model.options$n.clusters))
		parameters$nuggets <- rexp(n.ind)
		parameters$cluster.list <- populate.cluster.list(parameters$cluster.list,data,parameters$admix.proportions)
		parameters$admixed.covariance <- admixed.covariance(parameters$cluster.list,model.options$n.clusters,parameters$shared.mean)
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
	output.length <- mcmc.options$ngen/mcmc.options$samplefreq
	output.list <- list("step" = 1,
						"likelihood" = rep(NA,output.length),
						"posterior.prob" = rep(NA,output.length),
						"nuggets" = matrix(NA,nrow=n.ind,ncol=output.length),
						"admix.proportions" = lapply(1:output.length,"==",NA),
						"shared.mean" = rep(NA,output.length),
						"covariance.params" = list("geo.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"time.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"sill" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"exponent" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length)),
						"acceptance.rates" = list("geo.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"time.effect" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"sill" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"exponent" = matrix(NA,nrow=model.options$n.clusters,ncol=output.length),
													"nuggets" = matrix(NA,nrow=n.ind,ncol=output.length),
													"admix.proportions" = matrix(NA,nrow=n.ind,ncol=output.length)))
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
													"nuggets" = rep(NA,n.ind),
													"admix.proportions" = matrix(NA,nrow=n.ind,ncol=model.options$n.clusters),
													"shared.mean" = NA),
							"adaptive.mcmc" = list("log.stps" = 
														list("geo.effect" = rep(0,model.options$n.clusters),
															 "time.effect" = rep(0,model.options$n.clusters),
 															 "sill" = rep(0,model.options$n.clusters),
  															 "exponent" = rep(0,model.options$n.clusters),
															 "nuggets" = rep(0,n.ind),
															 "admix.proportions" = rep(0,n.ind),
															 "shared.mean" = 0),
													"nth.batch.accept.rates" = 
														list("geo.effect" = rep(0,model.options$n.clusters),
															 "time.effect" = rep(0,model.options$n.clusters),
 															 "sill" = rep(0,model.options$n.clusters),
  															 "exponent" = rep(0,model.options$n.clusters),
															 "nuggets" = rep(0,n.ind),
															 "admix.proportions" = rep(0,n.ind),
															 "shared.mean" = 0)),
							"covariance.params.list" = declare.covariance.params.list())
	class(mcmc.quantities) <- "mcmc.quantities"
	return(mcmc.quantities)
}

print.mcmc.quantities <- function(mcmc.quantities){
	print(str(mcmc.quantities))
}

declare.covariance.params.list <- function(){
	covariance.params.list <- list(
				"covariance.params" = list(
					"geo.effect",
					"time.effect",
					"sill",
					"exponent"),
				"covariance.param.prior.functions" = list(
					prior.prob.geo.effect,
					prior.prob.time.effect,
					prior.prob.sill,
					prior.prob.exponent))
	return(covariance.params.list)
}

prior.prob.param.list <- function(parameter.list,element,function.name){
	function.name(parameter.list[[element]])
}

prior.prob.nuggets <- function(nuggets){
	dexp(nuggets,log=TRUE)
}

prior.prob.admix.proportions <- function(admix.proportions){
	log(ddirichlet(admix.proportions,
					alpha=rep(0.5,2))) #WILL HAVE TO CHANGE THIS TO MAKE FLEXIBLE W/R/T ALPHA & N.IND
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

prior.prob.geo.effect <- function(geo.effect){
	dunif(geo.effect,0,1e20,log=TRUE)
}

prior.prob.time.effect <- function(time.effect){
	dunif(time.effect,0,1e20,log=TRUE)
}

prior.prob.sill <- function(sill){
	dunif(sill,0,1e20,log=TRUE)
}

prior.prob.exponent <- function(exponent){
	dexp(exponent,1,log=TRUE)
}

initialize.mcmc.quantities <- function(data.list,parameter.list,model.options,mcmc.options,initial.parameters=NULL){
	# recover()
	mcmc.quantities <- make.mcmc.quantities(data.list$n.ind,model.options,mcmc.options)
	# if(!is.null(initial.parameters)){
		# stop("not built yet")
	# } else {
#
		mcmc.quantities$likelihood <- calculate.likelihood.2(data.list, parameter.list$inverse,  parameter.list$determinant )     #calculate.likelihood(data.list,parameter.list)    #GRAHAM PUT THIS IN
		mcmc.quantities$prior.probs$nuggets <- prior.prob.param.list(parameter.list,"nuggets",prior.prob.nuggets)
		mcmc.quantities$prior.probs$admix.proportions <- prior.prob.param.list(parameter.list,"admix.proportions",prior.prob.admix.proportions)
		mcmc.quantities$prior.probs$shared.mean <- prior.prob.param.list(parameter.list,"shared.mean",prior.prob.shared.mean)
		mcmc.quantities$prior.probs$geo.effect <- prior.prob.covariance.param.clusters(parameter.list$cluster.list,"geo.effect",prior.prob.geo.effect)
		mcmc.quantities$prior.probs$time.effect <- prior.prob.covariance.param.clusters(parameter.list$cluster.list,"time.effect",prior.prob.time.effect)
		mcmc.quantities$prior.probs$sill <- prior.prob.covariance.param.clusters(parameter.list$cluster.list,"sill",prior.prob.sill)
		mcmc.quantities$prior.probs$exponent <- prior.prob.covariance.param.clusters(parameter.list$cluster.list,"exponent",prior.prob.exponent)
		mcmc.quantities$posterior.prob <- mcmc.quantities$likelihood + sum(unlist(mcmc.quantities$prior.probs))
	# }
	return(mcmc.quantities)
}

sherman_r <- function(Ap, u, v) {
  determ.update <- drop(1 + t(v) %*% Ap %*% u)
  
  list(new.inverse = Ap - (Ap %*% u %*% t(v) %*% Ap)/determ.update, determ.update=determ.update)
  
  }

update.cluster.mean.k<-function(k,data.list,super.list){ 
	delta.mean<-rnorm(1,sd=0.05) 
	 recover()
	new.cluster.mean <- delta.cluster.mean + super.list$parameter.list$cluster.list[[k]] 

	new.prior.prob <- prior.prob.nuggets(new.cluster.mean)   ###wanted to do this here, but perhaps not best way
	if(is.finite(new.prior.prob) ){
		u <- super.list$parameter.list$admix.proportions[,k] 
		v <- super.list$parameter.list$admix.proportions[,k] 
		
		u<- u * 	delta.mean
			
		inverse.updated<-sherman_r (super.list$parameter.list$inverse,u,v)   ##YOU DO THE SHERMAN TANK
		new.determinant <-  super.list$parameter.list$determinant + log(abs(inverse.updated$determ.update)) 
		new.inverse <- inverse.updated$new.inverse

		if(FALSE){ #TEST matrix inversion
			new.admixed.covariance<-super.list$parameter.list$admixed.covar
			new.admixed.covariance[i,i]<-new.admixed.covariance[i,i]+delta.nugget
			test <- solve(new.admixed.covariance)	 #GID CHECK
			summary(c(abs(test-new.inverse)))  #GID CHECK
		}
		
		new.likelihood <-   calculate.likelihood.2(data.list, new.inverse ,  new.determinant )[1]
		old.likelihood <- super.list$mcmc.quantities$likelihood
		likelihood.ratio <- new.likelihood - old.likelihood
		
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$nuggets[i]

		accept.ratio<-exp(likelihood.ratio +prior.ratio)
		if( accept.ratio >= runif(1)  ){
			super.list$mcmc.quantities$prior.probs$nuggets[i]  <- new.prior.prob 
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$nuggets[i]  <- new.nugget
			super.list$parameter.list$determinant <- new.determinant	
			super.list$parameter.list$inverse <- new.inverse 
		}
	}
	return(super.list)
}



update.nugget.i<-function(i,data.list,super.list){ 
	delta.nugget<-rnorm(1,sd=0.05) 
	 recover()
	new.nugget <- delta.nugget + super.list$parameter.list$nuggets[i] 

	new.prior.prob <- prior.prob.nuggets(new.nugget)   ###wanted to do this here, but perhaps not best way
	if(is.finite(new.prior.prob) ){
		u <- rep(0,  nrow(super.list$parameter.list$cluster.list[[1]]$covariance ))  ##WILL NEED NUM OF INDS HERE.
		u[i] <-  delta.nugget
		v <- rep(0,  nrow(super.list$parameter.list$cluster.list[[1]]$covariance )) ##WILL NEED NUM OF INDS HERE.
		v[i] <-  1
		
		inverse.updated<-sherman_r (super.list$parameter.list$inverse,u,v)   ##YOU DO THE SHERMAN TANK
		new.determinant <-  super.list$parameter.list$determinant + log(abs(inverse.updated$determ.update)) 
		new.inverse <- inverse.updated$new.inverse

		if(FALSE){ #TEST matrix inversion
			new.admixed.covariance<-super.list$parameter.list$admixed.covar
			new.admixed.covariance[i,i]<-new.admixed.covariance[i,i]+delta.nugget
			test <- solve(new.admixed.covariance)	 #GID CHECK
			summary(c(abs(test-new.inverse)))  #GID CHECK
		}
		
		new.likelihood <-   calculate.likelihood.2(data.list, new.inverse ,  new.determinant )[1]
		old.likelihood <- super.list$mcmc.quantities$likelihood
		likelihood.ratio <- new.likelihood - old.likelihood
		
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$nuggets[i]

		accept.ratio<-exp(likelihood.ratio +prior.ratio)
		if( accept.ratio >= runif(1)  ){
			super.list$mcmc.quantities$prior.probs$nuggets[i]  <- new.prior.prob 
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$nuggets[i]  <- new.nugget
			super.list$parameter.list$determinant <- new.determinant	
			super.list$parameter.list$inverse <- new.inverse 
		}
	}
	return(super.list)
}



update.shared.mean<-function(data.list,super.list){ 
	delta.shared.mean<-rnorm(1,sd=0.05) 
	# recover()
	new.shared.mean <- delta.shared.mean+ super.list$parameter.list$shared.mean
	new.prior.prob <- prior.prob.shared.mean(new.shared.mean)   ###wanted to do this here, but perhaps not best way
	if(is.finite(new.prior.prob) ){
		u = rep(delta.shared.mean,  nrow(super.list$parameter.list$cluster.list[[1]]$covariance ))  ##WILL NEED NUM OF INDS HERE.
		v = rep(1,  nrow(super.list$parameter.list$cluster.list[[1]]$covariance )) ##WILL NEED NUM OF INDS HERE.
		inverse.updated<-sherman_r (super.list$parameter.list$inverse,u,v)   
		new.determinant <-  super.list$parameter.list$determinant + log(abs(inverse.updated$determ.update)) 
		new.inverse <- inverse.updated$new.inverse

		if(FALSE){ #TEST matrix inversion
			new.admixed.covariance <- admixed.covariance(super.list$parameter.list$cluster.list,length(super.list$parameter.list$cluster.list),super.list$parameter.list$shared.mean+ delta.shared.mean)
			test <- solve(new.admixed.covariance)	 #GID CHECK
			summary(c(abs(test-new.inverse)))  #GID CHECK
		}
		
		new.likelihood <-   calculate.likelihood.2(data.list, new.inverse ,  new.determinant )[1]
		old.likelihood <- super.list$mcmc.quantities$likelihood
		likelihood.ratio <- new.likelihood - old.likelihood
		
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$shared.mean

		accept.ratio<-exp(likelihood.ratio +prior.ratio)
		if( accept.ratio >= runif(1)  ){
			super.list$mcmc.quantities$prior.probs$shared.mean  <- new.prior.prob 
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$shared.mean <- new.shared.mean
			super.list$parameter.list$determinant <- new.determinant	
			super.list$parameter.list$inverse <- new.inverse 
		}
	}
	return(super.list)
}

require("gtools")   ##for dirchlet, we could write our own
##update the w for the ith individual
update.w.i<-function( i, data.list,super.list){ #  i, data.list, parameter.list, model.options, mcmc.quantities){
	# recover()
	rejected.move <- 1
	num.clusters<-	length(super.list$parameter.list$cluster.list)
	these.two<-sample(num.clusters,2)
	clst.1<-these.two[1]
	clst.2<-these.two[2]
	
	old.w.1<-super.list$parameter.list$admix.proportions[i,clst.1]
	old.w.2<-super.list$parameter.list$admix.proportions[i,clst.2]


	delta.w<-rnorm(1,sd=0.05)  ##WILL NEED TO DEFINE
#	delta.w <- rnorm(1,sd= exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps$admix.proportions[i]))
	new.w.1 <- old.w.1 + delta.w
	new.w.2 <- old.w.2 - delta.w
	
	if( !(new.w.1<0 | new.w.1>1)  &  !(new.w.2<0 | new.w.2>1) )  {
		 #GID CHECK to see if graham inverse is same as gid inverse
		new.admix.proportions <- super.list$parameter.list$admix.proportions
		new.admix.proportions[i,these.two] <- c(new.w.1,new.w.2)
		covar.1 <- super.list$parameter.list$cluster.list[[clst.1]]$covariance[i,]  #will eventually need the cluster mean
		covar.2 <- super.list$parameter.list$cluster.list[[clst.2]]$covariance[i,]  #will eventually need the cluster mean
		
		u <- delta.w * ( super.list$parameter.list$admix.proportions[,clst.1] * covar.1  - super.list$parameter.list$admix.proportions[,clst.2] * covar.2  )
		u[i] <- u[i]  + delta.w^2 * (covar.1[i] + covar.2[i])/2 	### (wi^k + Dw) (wi^k + Dw) 
		v <- rep(0,data.list$n.ind)
		v[i] <- 1
		
		matrix.updated.once<-  sherman_r (super.list$parameter.list$inverse,u,v)  
		matrix.updated.twice<-  sherman_r (matrix.updated.once$new.inverse,v,u)  
		
		new.determinant <-  super.list$parameter.list$determinant + log(abs(matrix.updated.once$determ.update)) + log(abs(matrix.updated.twice$determ.update))
		new.inverse <- matrix.updated.twice$new.inverse
		
		if(FALSE){ #TEST matrix inversion
			new.admix.prop.matrices <- list(new.admix.proportions[,these.two[1]] %*% t(new.admix.proportions[,these.two[1]]),
											new.admix.proportions[,these.two[2]] %*% t(new.admix.proportions[,these.two[2]]))
			new.admixed.covariance <- update.admixed.covariance(super.list$parameter.list$cluster.list,num.clusters,
															super.list$parameter.list$shared.mean,
															updated.clusters = these.two,
															admix.prop.matrix.primes = new.admix.prop.matrices)
			test <- solve(new.admixed.covariance)	 #GID CHECK

			summary(c(abs(test-new.inverse)))  #GID CHECK
		}
		
#TEST MATRIX stuff, keep for mo.
#				determinant(parameter.list$admixed.covariance+u %*% t(v) + v %*% t(u) )
	#	all(abs(matrix.updated.twice$new.inverse - solve(parameter.list$admixed.covariance+u %*% t(v) + v %*% t(u) )<EPS))   #TEST
	#	calculate.likelihood.2(data.list, solve(parameter.list$admixed.covariance+u %*% t(v) + v %*% t(u) ) ,  determinant(parameter.list$admixed.covariance+u %*% t(v) + v %*% t(u) )$modulus )   #TEST
	#	calculate.likelihood.2(data.list, new.inverse ,  new.determinant )  #TEST

		new.likelihood <-   calculate.likelihood.2(data.list, new.inverse ,  new.determinant )[1]
		old.likelihood <- super.list$mcmc.quantities$likelihood
		likelihood.ratio <- new.likelihood - old.likelihood
		new.admixture.vec<-super.list$parameter.list$admix.proportions[i,]
		new.admixture.vec[these.two] <- c(new.w.1,new.w.2)
		
		new.prior.prob <- log(ddirichlet(new.admixture.vec ,rep(0.5,num.clusters) ))
		prior.ratio <- new.prior.prob - super.list$mcmc.quantities$prior.probs$admix.proportions[i]
		###prior.prob of admixture should be avector not a matrix
		accept.ratio<-exp(likelihood.ratio +prior.ratio)
		if( accept.ratio >= runif(1)  ){
			super.list$mcmc.quantities$prior.probs$admix.proportions[i]  <- new.prior.prob 
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- super.list$mcmc.quantities$posterior.prob + likelihood.ratio + prior.ratio
			super.list$parameter.list$admix.proportions[i,] <- new.admixture.vec
			super.list$parameter.list$determinant <- new.determinant
		
			super.list$parameter.list$cluster.list[[clst.1]]$admix.prop.matrix <-  super.list$parameter.list$admix.proportions[,clst.1] %*%t(super.list$parameter.list$admix.proportions[,clst.1] )
			super.list$parameter.list$cluster.list[[clst.2]]$admix.prop.matrix <-  super.list$parameter.list$admix.proportions[,clst.2] %*%t(super.list$parameter.list$admix.proportions[,clst.2])	
			##WE NEED TO UPDATE W MATRIX AS WELL BUT I DONT KNOW IF THAT"S WORTH WHILE DOIN HERE
			
			super.list$parameter.list$inverse <- new.inverse 
#			cat("updated ", i)	
		}
	}
	return(super.list)
}

#data <- geo.coords,time.coords,sample.covariance
#model.options <- round.earth,n.clusters
#mcmc.options <- ngen,samplefreq,printfreq
#initial.parameters <- initial.params, initial.params.lstps

MCMC.coop <- function(	data,
						model.options,
						mcmc.options,
						initial.parameters=NULL,
						seed=NULL){
	recover()
	if(is.null(seed)){
		seed <- sample(1:10000,1)
	}
	cat(seed)
	set.seed(seed)
	data.list <- make.data.list(data,model.options)
	super.list <- declare.super.list()
	super.list$parameter.list <- initialize.param.list(data.list,model.options,initial.parameters)
	super.list$output.list <- make.output.list(data.list$n.ind,model.options,mcmc.options)
	super.list$mcmc.quantities <- initialize.mcmc.quantities(data.list,super.list$parameter.list,model.options,mcmc.options,initial.parameters)
	super.list$model.options <- model.options
	tmp1 <- numeric(10000)
	tmp2 <- numeric(10000)
	super.list$mcmc.quantities$likelihood
	super.list$mcmc.quantities$posterior.prob

	quartz(width=8,height=5) ; par(mfrow=c(1,2))
	plot(c(data.list$sample.covariance),c(super.list$parameter.list$admixed.covariance)) ; abline(0,1,col="red")
	for(i in 1:20000){
		j <- sample(1:data.list$n.ind,1)
		super.list <-update.w.i( i=j, data.list,super.list)
		#super.list <- slow.update.w.j(j,data.list,super.list)
		tmp1[i] <- sum((super.list$parameter.list$admix.proportions - data$sim.admix.props)^2)
		tmp2[i] <- sum((super.list$parameter.list$admixed.covariance - data.list$sample.covariance)^2)
	}
	super.list$mcmc.quantities$likelihood
	super.list$mcmc.quantities$posterior.prob
	# plot(data.list$sample.covariance,super.list$parameter.list$admixed.covariance) ; abline(0,1,col="red")
#
	new.covar<-admixed.covariance(super.list$parameter.list$cluster.list,super.list$model.options$n.clusters,super.list$parameter.list$shared.mean)
	plot(data.list$sample.covariance,new.covar) ; abline(0,1,col="red")
	# plot(data.list$geo.dist,new.covar)
	# points(data.list$geo.dist,data.list$sample.covariance,col="red")
	# plot(super.list$parameter.list$admix.proportions,data$sim.admix.props,type="n"); text(super.list$parameter.list$admix.proportions,data$sim.admix.props, 1:10)
	# save(parameter.list,mcmc.quantities,file="~/desktop/testobj.Robj")
}

#Gid sims
slow.update.w.j <- function(j, data.list, super.list){
#	recover()
#	these.two <- sample(super.list$model.options$n.clusters,2)
#	clst.1 <- these.two[1]
#	clst.2 <- these.two[2]
#
	delta.w <- rnorm(1,sd= 0.01)#exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps$admix.proportions[i]))
	parameter.list.prime <- super.list$parameter.list
	mcmc.quantities.prime <- super.list$mcmc.quantities
	new.w.1 <- super.list$parameter.list$admix.proportions[j,1] + delta.w
	new.w.2 <- super.list$parameter.list$admix.proportions[j,2] - delta.w
#	
	if( !(new.w.1<0 | new.w.1>1)  &  !(new.w.2<0 | new.w.2>1) )  {
		parameter.list.prime$admix.proportions[j,1] <- new.w.1
		parameter.list.prime$admix.proportions[j,2] <- new.w.2
		parameter.list.prime$cluster.list$Cluster_1$admix.prop.matrix <- parameter.list.prime$admix.proportions[,1] %*% t(parameter.list.prime$admix.proportions[,1])
		parameter.list.prime$cluster.list$Cluster_2$admix.prop.matrix <- parameter.list.prime$admix.proportions[,2] %*% t(parameter.list.prime$admix.proportions[,2])
		parameter.list.prime$admixed.covariance <- admixed.covariance(parameter.list.prime$cluster.list,super.list$model.options$n.clusters,super.list$parameter.list$shared.mean)
		parameter.list.prime$inverse <- MASS::ginv(parameter.list.prime$admixed.covariance)
		parameter.list.prime$determinant <- determinant(parameter.list.prime$admixed.covariance,logarithm=TRUE)$modulus
		mcmc.quantities.prime$likelihood <- calculate.likelihood.2(data.list,parameter.list.prime$inverse,parameter.list.prime$determinant)
		mcmc.quantities.prime$prior.probs$admix.proportions <- prior.prob.admix.proportions(parameter.list.prime,super.list$model.options)
		mcmc.quantities.prime$posterior.prob <- mcmc.quantities.prime$likelihood + sum(unlist(mcmc.quantities.prime$prior.probs))
#		
		if(exp(mcmc.quantities.prime$posterior.prob - super.list$mcmc.quantities$posterior.prob) >= runif(1)){
			super.list$parameter.list <- parameter.list.prime
			super.list$mcmc.quantities <- mcmc.quantities.prime
		}
	}
	return(super.list)
}

propose.cov.param.update <- function(super.list,this.cluster,this.param){
	super.list$parameter.list$cluster.list[[this.cluster]]$covariance.params[[this.param]] + 
		rnorm(1,0,0.1)	#sd=exp(super.list$mcmc.quantities$adaptive.mcmc$log.stps[[this.param]][this.cluster]))
}

recalculate.cluster.list <- function(data.list,cluster.list,this.cluster,this.param,new.param){
	cluster.list[[this.cluster]]$covariance.params[[this.param]] <- new.param
	cluster.list[[this.cluster]]$covariance <- cluster.spatial.covariance(data.list$geo.dist,
																			data.list$time.dist,
																			cluster.list[[this.cluster]]$covariance.params)
	return(cluster.list)
}

recalculate.posterior.prob <- function(super.list,new.likelihood,new.param.prior,this.cluster,this.param){
	as.numeric(super.list$mcmc.quantities$posterior.prob - 
				super.list$mcmc.quantities$likelihood - 
				super.list$mcmc.quantities$prior.probs[[this.param]][this.cluster] +
				new.likelihood + new.param.prior)
}

update.cluster.covariance.param <- function(data.list,super.list){
	this.cluster <- sample(super.list$model.options$n.clusters,1)
	this.param <- sample(1:length(super.list$mcmc.quantities$covariance.params.list$covariance.params),1)
	this.prior <- super.list$mcmc.quantities$covariance.params.list$covariance.param.prior.functions[[this.param]]
	new.param <- propose.cov.param.update(super.list,this.cluster,this.param)
	new.param.prior <- this.prior(new.param)
	if(is.finite(new.param.prior)){
		new.cluster.list <- recalculate.cluster.list(data.list,super.list$parameter.list$cluster.list,this.cluster,this.param,new.param)
		new.admixed.covariance <- admixed.covariance(new.cluster.list,super.list$model.options$n.clusters,super.list$parameter.list$shared.mean)
		new.inverse <- MASS::ginv(new.admixed.covariance)
		new.determinant <- determinant(new.admixed.covariance,logarithm=TRUE)$modulus
		new.likelihood <- calculate.likelihood.2(data.list, new.inverse, new.determinant)
		new.posterior.prob <- recalculate.posterior.prob(super.list,new.likelihood,new.param.prior,this.cluster,this.param)
		if(exp(new.posterior.prob - super.list$mcmc.quantities$posterior.prob) >= runif(1)){
			super.list$parameter.list$cluster.list <- new.cluster.list
			super.list$parameter.list$admixed.covariance <- new.admixed.covariance
			super.list$parameter.list$inverse <- new.inverse
			super.list$parameter.list$determinant <- new.determinant
			super.list$mcmc.quantities$prior.probs[[this.param]][this.cluster] <- new.param.prior
			super.list$mcmc.quantities$likelihood <- new.likelihood
			super.list$mcmc.quantities$posterior.prob <- new.posterior.prob
		}
	}
	return(super.list)
}

index.cluster.param <- function(cluster.list,param){
	unlist(lapply(cluster.list,FUN=function(cluster.list){cluster.list$covariance.params[[param]]}))
}

bookkeep <- function(super.list){
	step <- super.list$output.list$step
	super.list$output.list$likelihood[step] <- super.list$mcmc.quantities$likelihood
	super.list$output.list$posterior.prob[step] <- super.list$mcmc.quantities$posterior.prob
	super.list$output.list$shared.mean[step] <- super.list$parameter.list$shared.mean
	super.list$output.list$admix.proportions[[step]] <- super.list$parameter.list$admix.proportions
	super.list$output.list$nuggets[,step] <- super.list$parameter.list$nuggets
	super.list$output.list$covariance.params$geo.effect[,step] <- index.cluster.param(super.list$parameter.list$cluster.list,"geo.effect")
	super.list$output.list$covariance.params$time.effect[,step] <- index.cluster.param(super.list$parameter.list$cluster.list,"time.effect")
	super.list$output.list$covariance.params$exponent[,step] <- index.cluster.param(super.list$parameter.list$cluster.list,"exponent")
	super.list$output.list$covariance.params$sill[,step] <- index.cluster.param(super.list$parameter.list$cluster.list,"sill")
	super.list$output.list$step <- step + 1
	return(super.list)
}

#GIDEON PLAYGROUND
MCMC.gid <- function(	data,
						model.options,
						mcmc.options,
						initial.parameters=NULL,
						seed=NULL){
	# recover()
	if(is.null(seed)){
		seed <- sample(1:10000,1)
	}
	cat(seed)
	set.seed(seed)
	data.list <- make.data.list(data,model.options)
	super.list <- initialize.super.list(data.list,model.options,mcmc.options,initial.parameters)
	for(i in 2:mcmc.options$ngen){
		if(i%%mcmc.options$samplefreq==0){
			super.list <- update.shared.mean(data.list,super.list)	#update.cluster.covariance.param(data.list,super.list)
			super.list <- bookkeep(super.list)
		}
	}
	save(super.list,file="~/desktop/testobj.Robj")
}

#SIMULATE TOY DATA
k <- 20
n.loci <- 1e5
spatial.coords <- cbind(runif(k),runif(k))
temporal.coords <- rep(0,k)  # sample(1:100,k,replace=TRUE)  #GRAHAM CHANGES
geo.dist <- fields::rdist(spatial.coords)
time.dist <- fields::rdist(temporal.coords)
sim.admix.props <- sample(0:1,size=k,replace=TRUE)
	sim.admix.props[which(sim.admix.props==1)] <- 0.99
	sim.admix.props[which(sim.admix.props==0)] <- 0.01
	sim.admix.props <- cbind(sim.admix.props,1-sim.admix.props)
sim.cluster.list <- list("Cluster_1" = list("admix.prop.matrix" = sim.admix.props[,1] %*% t(sim.admix.props[,1]),
											"covariance" = spatial.covariance(geo.dist,time.dist,2,0.3,1.1,2)),
						 "Cluster_2" = list("admix.prop.matrix" = sim.admix.props[,2] %*% t(sim.admix.props[,2]),
											"covariance" = spatial.covariance(geo.dist,time.dist,1,1.2,1.3,2.5)))
sim.admixed.cov.mat <- admixed.covariance(sim.cluster.list,2,1.5)
											
#generate initial parameters 
# that are the same as those used to simulate data,
# EXCEPT for the admixture proportions, which are simulated from a dirichlet
# replicate(10,{
fake.admix.props <- sim.admix.props
#fake.admix.props <- gtools::rdirichlet(n = k,alpha = rep(1,2))
#fake.admix.props[3:k,] <- sim.admix.props[3:k,] #sets all but one ind's admix prop to true value, to check label-switchery
initial.parameters <- list("shared.mean" = 0,
							"admix.proportions" = fake.admix.props,
							"nuggets" = rep(0,k),
							"cluster.list" = generate.clusters(2))
#initial parameters cluster 1
	initial.parameters$cluster.list$Cluster_1$covariance.params$geo.effect <- 2
	initial.parameters$cluster.list$Cluster_1$covariance.params$time.effect <- 0.3
	initial.parameters$cluster.list$Cluster_1$covariance.params$exponent <- 1.1
	initial.parameters$cluster.list$Cluster_1$covariance.params$sill <- 2
	initial.parameters$cluster.list$Cluster_1$covariance <- spatial.covariance(geo.dist,time.dist,2,0.3,1.1,2)
	initial.parameters$cluster.list$Cluster_1$admix.prop.matrix <- fake.admix.props[,1] %*% t(fake.admix.props[,1])

#initial parameters cluster 2
	initial.parameters$cluster.list$Cluster_2$covariance.params$geo.effect <- 1
	initial.parameters$cluster.list$Cluster_2$covariance.params$time.effect <- 1.2
	initial.parameters$cluster.list$Cluster_2$covariance.params$exponent <- 1.3
	initial.parameters$cluster.list$Cluster_2$covariance.params$sill <- 2.5
	initial.parameters$cluster.list$Cluster_2$covariance <- spatial.covariance(geo.dist,time.dist,1,1.2,1.3,2.5)
	initial.parameters$cluster.list$Cluster_2$admix.prop.matrix <- fake.admix.props[,2] %*% t(fake.admix.props[,2])




sim.data <- list("geo.coords" = spatial.coords,
				"time.coords" = temporal.coords,
				"sample.covariance" = sim.admixed.cov.mat,
				"n.loci" = n.loci,
				"sim.admix.props" = sim.admix.props)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 2)
mcmc.options = list("ngen" = 1e6,
					"samplefreq" = 1000,
					"printfreq" = 5)

MCMC.coop(sim.data,model.options,mcmc.options,initial.parameters)

MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters)
load("~/Desktop/testobj.Robj")






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
tmpfunc <- function(list1,list2){priorprob <- dexp(list1$geo.effect,log=TRUE) ; }
var(colMeans(sim.freqs))
var(sim.means)
quartz()
plot(geo.dist,cov(t(sim.freqs)))
points(geo.dist,par.cov+var(sim.means),col="red",pch=20)
initiate.adaptive.mcmc <- function(mcmc.quantities){
	mcmc.quantities$adaptive.mcmc$log.stps <- rapply(mcmc.quantities$adaptive.mcmc$log.stps,function(x){0},how="replace")
	return(mcmc.quantities)
}
slow.admixed.covariance <- function(cluster.list,n.clusters,shared.mean){
	admixed.cov <- cluster.list[[1]]$covariance * cluster.list[[1]]$admix.prop.matrix
	for(i in 2:n.clusters){
		admixed.cov <- admixed.cov + cluster.list[[i]]$covariance * cluster.list[[i]]$admix.prop.matrix
	}
	return(admixed.cov)
}
slow.update.w.i <- function(i, data.list, parameter.list, model.options, mcmc.quantities){
	recover()
	these.two <- sample(model.options$n.clusters,2)
	clst.1 <- these.two[1]
	clst.2 <- these.two[2]

	delta.w <- rnorm(1,sd= exp(mcmc.quantities$adaptive.mcmc$log.stps$admix.proportions[i]))
	new.w.1 <- parameter.list$admix.proportions[i,clst.1] + delta.w
	new.w.2 <- parameter.list$admix.proportions[i,clst.2] - delta.w
	
	if( !(new.w.1<0 | new.w.1>1)  &  !(new.w.2<0 | new.w.2>1) )  {
		admix.proportions.prime <- parameter.list$admix.proportions
		admix.proportions.prime[i,clst.1] <- new.w.1
		admix.proportions.prime[i,clst.2] <- new.w.2
		admixed.covariance.prime <- admixed.covariance(parameter.list$cluster.list,model.options$n.clusters,parameter.list$shared.mean)
}}
# Update <- function(super.list){
		
 }
 }
k <- 100
n.loci <- 1e4
spatial.coords <- cbind(runif(k),runif(k))
temporal.coords <- sample(1:100,k,replace=TRUE)
geo.dist <- fields::rdist(spatial.coords)
time.dist <- fields::rdist(temporal.coords)

par.cov <- spatial.covariance(geo.dist,time.dist,2,3,1,2)
# sim.means <- rbeta(n.loci,0.10,0.10)
# sim.freqs <- lapply(1:n.loci, FUN=function(i){MASS::mvrnorm(n = 1, mu = rep(sim.means[i],k), Sigma = par.cov)})
# sim.freqs <- do.call(cbind,sim.freqs)
# sample.covariance <- cov(t(sim.freqs))

admix<-sample(0:1,size=k,replace=TRUE)
par.cov.underlying<-par.cov
par.cov<-admix %*% t(admix) * par.cov + (1-admix) %*% t(1-admix) * par.cov

sim.data <- list("geo.coords" = spatial.coords,
			"time.coords" = temporal.coords,
			"sample.covariance" = par.cov,
			"n.loci" = n.loci,
			"true.spatial.covar"=par.cov.underlying
			)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 2)
mcmc.options = list("ngen" = 100,
					"samplefreq" = 10,
					"printfreq" = 5)
MCMC.coop(sim.data,model.options,mcmc.options)

	tmp1 <- numeric(10000)
	tmp2 <- numeric(10000)
	tmp3 <- numeric(10000)
	tmp4 <- numeric(10000)
	super.list$mcmc.quantities$likelihood
	super.list$mcmc.quantities$posterior.prob
	par(mfrow=c(2,2))
		plot(data.list$sample.covariance,
				super.list$parameter.list$admixed.covariance,
				xlab="sim.cov",
				ylab="est.cov",
				main="Initial: sim vs. est cov") ; abline(0,1,col="red")
	for(i in 1:10000){
		j <- sample(1:data.list$n.ind,1)
		#super.list <-update.w.i( i=j, data.list,super.list)
		# super.list <- slow.update.w.j(j,data.list,super.list)
		super.list <- update.cluster.covariance.param(data.list,super.list)
		tmp1[i] <- super.list$mcmc.quantities$likelihood	#sum((super.list$parameter.list$admix.proportions - data$sim.admix.props)^2)
		tmp2[i] <- sum((super.list$parameter.list$admixed.covariance - data.list$sample.covariance)^2)
		tmp3[i] <- super.list$parameter.list$cluster.list$Cluster_1$covariance.params[[1]]
		tmp4[i] <- super.list$parameter.list$cluster.list$Cluster_2$covariance.params[[1]]
	}
	super.list$mcmc.quantities$likelihood
	super.list$mcmc.quantities$posterior.prob
	plot(data.list$sample.covariance,
		super.list$parameter.list$admixed.covariance,
				xlab="sim.cov",
				ylab="est.cov",
				main="Post-run: sim vs. est cov") ; abline(0,1,col="red")
	plot(tmp3,xlab="mcmc iterations",ylab = "sum sq. diff",main="between est and sim admix props") ; abline(h=2,col="red")
	plot(tmp4,xlab="mcmc iterations",ylab = "sum sq. diff", main="between est and sim admixed cov")
	# for(i in 1:mcmc.options$ngen){
		# super.list <- Update(super.list)
		
	# }
plot.output.obj <- function(output.obj){
	for(i in 1:nrow(output.obj)){
		lines(output.obj[i,],col=i)
	}
}
plot(0,xlim=c(0,1000),ylim=c(0,5),type='n',xlab="",ylab="")
lapply(super.list$output.list$covariance.params,plot.output.obj)
points(super.list$output.list$covariance.params$geo.effect[1,],pch=20,col="blue")
}