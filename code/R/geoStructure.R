validate.data.list <- function(data.block){
	if(!"spatial" %in% names(data.block)){
		stop("\nUser must specify a \"spatial\" option\n\n")
	}
	if(!"N" %in% names(data.block)){
		stop("\nUser must specify a \"N\"\n\n")
	}
	if(!"K" %in% names(data.block)){
		stop("\nUser must specify a \"K\"\n\n")
	}
	if(!"L" %in% names(data.block)){
		stop("\nUser must specify a \"L\"\n\n")
	}
	if(!"obsSigma" %in% names(data.block)){
		stop("\nUser must specify a \"obsSigma\"\n\n")
	}
	if(!"projMat" %in% names(data.block)){
		stop("\nUser must specify a \"projMat\"\n\n")
	}
	if(!"geoDist" %in% names(data.block)){
		stop("\nUser must specify a \"geoDist\"\n\n")
	}
	if(!"sampleSize" %in% names(data.block)){
		stop("\nUser must specify a \"sampleSize\"\n\n")
	}
	return(invisible("list elements validated"))	
}

validate.n.samples <- function(data.block){
	n.samples <- data.block$N
	n.samples <- c(data.block$N,nrow(data.block$obsSigma)+1)
	n.samples <- c(n.samples,length(data.block$sampleSize))
	if(!is.null(data.block$geoDist)){
		n.samples <- c(n.samples,nrow(data.block$geoDist))
	}
	if(length(unique(n.samples)) > 1){
		stop("\nthe number of samples is not consistent 
				across entries in the data.block\n\n")
	}
	return(invisible("n.samples validated"))
}

validate.model <- function(data.block){
	if(data.block$spatial){
		if(is.null(data.block$geoDist)){
			stop("\nyou have specified a spatial model,
				  but you have not specified a matrix 
				  of pairwise geographic distances")
		}
	}
	return(invisible("model validated"))
}

make.data.block.S3 <- function(data.block){
	data.block <- data.block
	class(data.block) <- "data.block"
	return(data.block)
}

print.data.block <- function(data.block){
	print(str(data.block,max.level=1))
}

validate.data.block <- function(data.block){
	message("\nchecking data.block\n")
		validate.data.list(data.block)
		validate.n.samples(data.block)
	message(sprintf("\treading %s samples",data.block$N))
	message(sprintf("\treading %s loci",data.block$L))
	message("\nchecking specified model\n")
		validate.model(data.block)
	message(sprintf("\treading %s cluster(s)",data.block$K))
	if(data.block$spatial){
		message(sprintf("\nuser has specified a spatial model\n"))
	}
	if(!data.block$spatial){
		message(sprintf("\nuser has specified a purely discrete model\n"))
	}
	data.block <- make.data.block.S3(data.block)
	return(data.block)
}

make.stan.code.block <- function(spatial,n.clusters){
	stan.code.block.name <- "stan.block"
	if(n.clusters == 1){
		stan.code.block.name<- paste0("oneK.",stan.code.block.name)
	}
	if(n.clusters > 1){
		stan.code.block.name<- paste0("multiK.",stan.code.block.name)
	}
	if(spatial){
		stan.code.block.name <- paste0("space.",stan.code.block.name)
	}
	return(get(stan.code.block.name))
}

get.projection.matrix <- function(mean.sample.sizes){
	k <- length(mean.sample.sizes)
	transformation.matrix <- get.transformation.matrix(mean.sample.sizes)
	qr.transformation.matrix <- qr(t(transformation.matrix))
	projection.matrix <- qr.Q(qr.transformation.matrix)[,1:qr.transformation.matrix$rank]
	stopifnot(qr.transformation.matrix$rank == sum(abs(eigen(transformation.matrix)$values - 1) < 1e-2) )
	return(projection.matrix)
}

get.transformation.matrix <- function(mean.sample.sizes){
	k <- length(mean.sample.sizes)
	transformation.matrix <- diag(k) - matrix(1/k,nrow=k,ncol=k,byrow=TRUE)
#	transformation.matrix <- diag(k) - matrix(mean.sample.sizes/(sum(mean.sample.sizes)),nrow=k,ncol=k,byrow=TRUE)
	return(transformation.matrix)
}

project.sample.covariance <- function(sample.covariance,mean.sample.sizes){
	projection.matrix <- get.projection.matrix(mean.sample.sizes)
	sample.covariance <- t(projection.matrix) %*% sample.covariance %*% projection.matrix
	return(sample.covariance)
}

get.norm.factor <- function(freqs,n.loci){
	pseudo.means <- colMeans(rbind(freqs,rep(0.5,n.loci)))
	norm.factor <- sqrt(pseudo.means * (1-pseudo.means))
	norm.factor <- matrix(norm.factor,nrow=nrow(freqs),ncol=n.loci,byrow=TRUE)
	return(norm.factor)
}

get.mean.freqs <- function(freqs,n.loci){
	mean.freqs <- matrix(colMeans(freqs),nrow=nrow(freqs),ncol=n.loci,byrow=TRUE)
	return(mean.freqs)
}

make.freq.data.list.S3 <- function(freq.data){
	freq.data <- freq.data
	class(freq.data) <- "freq.data"
	return(freq.data)
}

print.freq.data <- function(freq.data){
	print(str(freq.data,max.level=1))
}

standardize.freqs <- function(freqs,sample.sizes){
	invars <- apply(freqs,2,function(x){length(unique(x))==1})	
	freqs <- freqs[,!invars]
	n.loci <- ncol(freqs)
	proj.mat <- get.projection.matrix(sample.sizes)	
	norm.factor <- get.norm.factor(freqs,n.loci)
	mean.freqs <- get.mean.freqs(freqs,n.loci)
	mc.freqs <- freqs-mean.freqs
	norm.freqs <- freqs/norm.factor
	std.freqs <- (freqs-mean.freqs)/norm.factor
	sample.cov <- cov(t(freqs),use="pairwise.complete.obs")
	mc.cov <- cov(t(mc.freqs),use="pairwise.complete.obs")
	norm.cov <- cov(t(norm.freqs),use="pairwise.complete.obs")
	std.cov <- cov(t(std.freqs),use="pairwise.complete.obs")
	proj.std.cov <- project.sample.covariance(std.cov,sample.sizes)
	std.freq.list <- list("freqs" = freqs,
						  "norm.freqs" = norm.freqs,
						  "std.freqs" = std.freqs,
						  "norm.factor" = norm.factor,
						  "mean.freqs" = mean.freqs)
	std.cov.list <- list("sample.cov" = sample.cov,
					   	 "mc.cov" = mc.cov,
					   	 "norm.cov" = norm.cov,
					   	 "std.cov" = std.cov)
	proj.cov.list <- list("proj.mat" = proj.mat,
					   	 "proj.std.cov" = proj.std.cov)
	freq.data.list <- list("n.loci" = n.loci,
						   "std.freq.list" = std.freq.list,
						   "std.cov.list" = std.cov.list,
						   "proj.cov.list" = proj.cov.list)
	freq.data.list <- make.freq.data.list.S3(freq.data.list)
	return(freq.data.list)
}

make.data.block <- function(K,std.freq.list,D,coords,sample.sizes,prefix,spatial){
	data.block <- list( "K" = K,
						"N" = nrow(D),
						"L" = std.freq.list$n.loci-1,
						"obsSigma" = (std.freq.list$n.loci-1) * std.freq.list$proj.cov.list$proj.std.cov,
						"geoDist" = D,
						"projMat" = std.freq.list$proj.cov.list$proj.mat,
						"sampleSize" = sample.sizes,
						"spatial" = spatial,
						"coords" = coords)
	data.block <- validate.data.block(data.block)
	return(data.block)
}

geoStructure <- function(spatial=TRUE,K,freqs,D,coords=NULL,sample.sizes,prefix,n.chains=1,n.iter=1e4,...){
	#validate data block
	std.freq.list <- standardize.freqs(freqs,sample.sizes)
		save(std.freq.list,file=paste0(prefix,"_std.freq.list.Robj"))
	data.block <- make.data.block(K,std.freq.list,D,coords,sample.sizes,prefix,spatial)
		save(data.block,file=paste0(prefix,"_data.block.Robj"))
	#validate model specification
	#make.stan.code.block
	stan.block <- make.stan.code.block(spatial,K)
		#write stan block to file
	#run model
	#put stan in tryCatch, email me
	require(rstan)
	model.fit <- stan(model_code = stan.block,
						data = data.block,
						iter = n.iter,
						chains = n.chains,
						thin = n.iter/500,
						...)
	#save fit obj
	save(model.fit,file=paste(prefix,"model.fit.Robj",sep="_"))
	geoStr.results <- get.geoStructure.results(data.block,model.fit,n.chains)
	save(geoStr.results,file=paste(prefix,"geoStr.results.Robj",sep="_"))
	make.all.the.plots(geoStr.results,n.chains,data.block,std.freq.list,prefix,burnin=100,cluster.colors=NULL)
	return(geoStr.results)
}

get.geoStructure.results <- function(data.block,model.fit,n.chains){
	geoStr.results <- setNames(
						lapply(1:n.chains,
							function(i){
								get.geoStructure.chain.results(data.block,model.fit,i)
							}),
					  paste0("chain_",1:n.chains))
	return(geoStr.results)
}

get.MAP.iter <- function(model.fit,chain.no){
	logpost <- get_logposterior(model.fit)
	MAP.iter <- lapply(logpost,which.max)[[chain.no]]
	return(MAP.iter)
}

get.admix.props <- function(model.fit,chain.no,N,n.clusters){
	# recover()
	admix.props <- array(1,dim=c(length(get_logposterior(model.fit)[[chain.no]]),N,n.clusters))
	if(any(grepl("w",model.fit@model_pars))){
		for(k in 1:n.clusters){
			admix.props[,,k] <- extract(model.fit,
										pars=unlist(lapply(1:N,function(j){sprintf("w[%s,%s]",j,k)})),
										permute=FALSE,inc_warmup=TRUE)[,chain.no,]
		}
	}
	return(admix.props)
}

get.par.cov <- function(model.fit,chain.no,N){
	par.cov <- array(NA,dim=c(length(get_logposterior(model.fit)[[chain.no]]),N,N))
	for(i in 1:N){
		for(j in 1:N){
			my.par <- sprintf("Sigma[%s,%s]",i,j)
			par.cov[,i,j] <- extract(model.fit,pars=my.par,inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
		}
	}
	return(par.cov)
}

get.nuggets <- function(model.fit,chain.no,N){
	nuggets <- extract(model.fit,pars="nugget",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(nuggets)
}

get.alpha.params <- function(model.fit,chain.no,cluster,n.clusters){
	alpha.pars <- model.fit@model_pars[grepl("alpha",model.fit@model_pars)]
	if(n.clusters > 1){
		alpha.params <- setNames(
							lapply(1:length(alpha.pars),
									function(i){
										extract(model.fit,
												pars=paste0(alpha.pars[i],"[",cluster,"]"),
												inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
									}),alpha.pars)
	} else {
		alpha.params <- setNames(
							lapply(1:length(alpha.pars),
									function(i){
										extract(model.fit,
												pars=alpha.pars[i],
												inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
									}),alpha.pars)		
	}
	return(alpha.params)
}

get.cluster.mu <- function(model.fit,chain.no,cluster){
	mu <- extract(model.fit,
					pars=paste0("mu","[",cluster,"]"),
					inc_warmup=TRUE,permute=FALSE)[,chain.no,]
	return(mu)
}

get.mu <- function(model.fit,chain.no){
	mu <- extract(model.fit,pars="mu",
					inc_warmup=TRUE,permute=FALSE)[,chain.no,]
	return(mu)
}

get.cluster.DirichAlpha <- function(model.fit,chain.no,cluster){
	DirichAlpha <- extract(model.fit,
						   pars=paste0("DirichAlpha","[",cluster,"]"),
						   inc_warmup=TRUE,permute=FALSE)[,chain.no,]
	return(DirichAlpha)
}

get.cluster.dirichHP <- function(model.fit,chain.no,cluster){
	dirichHP <- extract(model.fit,
						 pars=paste0("dirichHP","[",cluster,"]"),
						 inc_warmup=TRUE,permute=FALSE)[,chain.no,]
	return(dirichHP)
}

get.cov.function <- function(data.block){
	if(data.block$K == 1){
		if(data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(cluster.params$alphaD*data.block$geoDist)^cluster.params$alpha2))
			}
		}
		if(!data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(matrix(0,nrow=data.block$N,ncol=data.block$N))
			}
		}
	} else {
		if(data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(cluster.params$alphaD*data.block$geoDist)^cluster.params$alpha2) + 
						cluster.params$mu)
			}
		}
		if(!data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(matrix(cluster.params$mu,nrow=data.block$N,ncol=data.block$N))
			}
		}
	}
	return(cov.func)
}

get.cluster.cov <- function(cluster.params,data.block,n.iter){
	cov.function <- get.cov.function(data.block)
	cluster.cov <- lapply(1:n.iter,
							function(i){
								cov.function(cluster.params=
												lapply(cluster.params,"[[",i),
												data.block)
							})
	return(cluster.cov)
}

get.cluster.params <- function(model.fit,data.block,chain.no,cluster,n.clusters,n.iter){
	# recover()
	cluster.params <- list()
	if(data.block$spatial){
		cluster.params <- get.alpha.params(model.fit,chain.no,cluster,n.clusters)
	}
	if(n.clusters > 1){
		cluster.params <- c(cluster.params,
							list("mu" = get.cluster.mu(model.fit,chain.no,cluster)),
							list("DirichAlpha" = get.cluster.DirichAlpha(model.fit,chain.no,cluster)),
							list("dirichHP" = get.cluster.dirichHP(model.fit,chain.no,cluster)))
	}
	cluster.cov <- get.cluster.cov(cluster.params,data.block,n.iter)
	cluster.params <- c(cluster.params,list("cluster.cov"=cluster.cov))
	return(cluster.params)
}

get.cluster.params.list <- function(model.fit,data.block,chain.no,n.iter){
	cluster.params <- setNames(
						lapply(1:data.block$K,
									function(i){
										get.cluster.params(model.fit,data.block,chain.no,i,data.block$K,n.iter)
									}),
						paste("Cluster",1:data.block$K,sep="_"))
	cluster.params <- make.cluster.params.S3(cluster.params)
	return(cluster.params)
}

make.cluster.params.S3 <- function(cluster.params){
	cluster.params <- cluster.params
	class(cluster.params) <- "cluster.params"
	return(cluster.params)
}

print.cluster.params <- function(cluster.params){
	print(str(cluster.params,max.level=1))
}

index.MAP <- function(param,MAP.iter){
	if(class(param) == "numeric"){
		MAP.param <- param[MAP.iter]
	}
	if(class(param) == "list"){
		MAP.param <- param[[MAP.iter]]
	}
	if(class(param) == "array"){
		MAP.param <- param[MAP.iter,,]
	}
	if(class(param) == "matrix"){
		MAP.param <- param[MAP.iter,]
	}
	if(class(param) == "cluster.params"){
		MAP.param <- index.MAP.cluster.params.list(param,MAP.iter)
	}
	return(MAP.param)
}

index.MAP.cluster.params <- function(cluster.params,MAP.iter){
	MAP.cluster.params <- lapply(cluster.params,index.MAP,MAP.iter)
	return(MAP.cluster.params)
}

index.MAP.cluster.params.list <- function(cluster.params.list,MAP.iter){
	MAP.cluster.params.list <- lapply(cluster.params.list,index.MAP.cluster.params,MAP.iter)
	return(MAP.cluster.params.list)
}

get.n.iter <- function(model.fit,chain.no){
	n.iter <- length(get_logposterior(model.fit)[[chain.no]])
	return(n.iter)
}

make.geoStructure.results.S3 <- function(geoStructure.results){
	geoStructure.results <- geoStructure.results
	class(geoStructure.results) <- "geoStructure.results"
	return(geoStructure.results)
}

print.geoStructure.results <- function(geoStructure.results){
	print(str(geoStructure.results,max.level=1))
}

mc.par.cov.post <- function(sample.sizes,par.cov.post){
	MC.mat <- get.transformation.matrix(sample.sizes)
	mc.par.cov <- array(NA,dim=dim(par.cov.post))
	for(i in 1:dim(par.cov.post)[1]){
		mc.par.cov[i,,] <- MC.mat %*% par.cov.post[i,,] %*% t(MC.mat)
	}
	return(mc.par.cov)
}

get.geoStructure.chain.results <- function(data.block,model.fit,chain.no){
	n.iter <- get.n.iter(model.fit,chain.no)
	post <- list("posterior" = get_logposterior(model.fit)[[chain.no]],
				 "nuggets" = get.nuggets(model.fit,chain.no,data.block$N),
				 "par.cov" = get.par.cov(model.fit,chain.no,data.block$N))
	post[["mc.par.cov"]] <- mc.par.cov.post(data.block$sampleSize,post$par.cov)
	if(data.block$spatial | data.block$K > 1){
		post[["cluster.params"]] <- get.cluster.params.list(model.fit,data.block,chain.no,n.iter)
	} else {
		post[["mu"]] <- get.mu(model.fit,chain.no)
	}
	if(data.block$K > 1){
		post[["admix.proportions"]] <- get.admix.props(model.fit,chain.no,data.block$N,data.block$K)
	}
	MAP.iter <- get.MAP.iter(model.fit,chain.no)
	point <- lapply(post,function(X){index.MAP(X,MAP.iter)})
	geoStructure.results <- list("post" = post,"point" = point)
	geoStructure.results <- make.geoStructure.results.S3(geoStructure.results)
	return(geoStructure.results)
}

post.process.par.cov <- function(geoStr.results,samples){
	pp.cov.list <- lapply(samples,
							function(i){
								list("inv" = chol2inv(chol(geoStr.results$post$par.cov[i,,])),
									 "log.det" = determinant(geoStr.results$post$par.cov[i,,])$modulus[[1]])
							})
	return(pp.cov.list)
}

calculate.likelihood <- function(obsSigma,inv.par.cov,log.det,n.loci){
	#recover()
	lnL <- -0.5 * (sum( inv.par.cov * obsSigma) + n.loci * log.det)
	return(lnL)
}

determine.log.shift <- function(chunk.lnls){
	if(diff(range(unlist(chunk.lnls))) > 700){
		message("the difference between the min and max lnLs may be inducing underflow")
	}
	return(max(unlist(chunk.lnls)))
}

shift.chunk.lnls <- function(chunk.lnls,A,n.iter){
	shift.chunk.lnls <- log(sum(exp(chunk.lnls-A))) + A - log(n.iter)
	return(shift.chunk.lnls)
}

calc.lnl.x.MCMC <- function(cov.chunk,pp.par.cov){
	#recover()
	lnl.x.mcmc <- lapply(pp.par.cov,
						function(x){
							calculate.likelihood(cov.chunk,x$inv,x$log.det,n.loci=1)
					})
	return(unlist(lnl.x.mcmc))
}

chunk.freq.data <- function(freqs){
	n.loci <- ncol(freqs)
	chunks <- lapply(1:n.loci,
					  function(i){
					  	freqs[,i,drop=FALSE] %*% t(freqs[,i,drop=FALSE])
					  })
	return(chunks)
}

calculate.lpd <- function(chunk.lnls,n.iter){
	#recover()
	#subtract max lnL from log likelihood to avoid overflow
		A <- determine.log.shift(chunk.lnls)
		chunk.lnls.shifted <- lapply(chunk.lnls,function(x){shift.chunk.lnls(x,A,n.iter)})
	lpd <- sum(unlist(chunk.lnls.shifted))
	return(lpd)
}

calculate.pwaic <- function(chunk.lnls){
	pwaic <- lapply(chunk.lnls,var)
	return(sum(unlist(pwaic)))
}

calculate.waic <- function(freqs,geoStr.results,samples=NULL){
	cat("breaking data into locus-by-locus covariances...\n\n")
	chunks <- chunk.freq.data(freqs)
	if(is.null(samples)){
		samples <- 1:length(geoStr.results$post$posterior)
	}
	n.iter <- length(samples)
	#invert the posterior distn of parametric cov matrices
	cat("inverting posterior distribution of parametric covariance matrices...\n\n")
		pp.par.cov <- post.process.par.cov(geoStr.results,samples)
	#calc likelioods of nth data chunk across all sampled MCMC iterations
	cat("calculating likelihood of each site across posterior distribution of parameters...\n\n")
		chunk.lnls <- lapply(chunks,function(x){calc.lnl.x.MCMC(x,pp.par.cov=pp.par.cov)})
	#calculate log pointwise predictive density
	cat("calculating wAIC score...\n\n\n")
		lpd <- calculate.lpd(chunk.lnls,n.iter)
	#calculate effective number of parameters
	pwaic <- calculate.pwaic(chunk.lnls)
	elpd <- lpd - pwaic
	waic <- -2 * elpd
	return(waic)
}

plot.prob <- function(geoStr.results,burnin=0){
	n.iter <- length(geoStr.results$post$posterior)
	z <- (burnin+1):n.iter
	plot(geoStr.results$post$posterior[z],
			xlab="MCMC iterations",ylab="posterior probability",
			main="Posterior probability",type='l')
	return(invisible(0))
}

plot.nuggets <- function(geoStr.results,burnin){
	n.iter <- length(geoStr.results$post$posterior)
	z <- (burnin+1):n.iter
	matplot(geoStr.results$post$nuggets[z,],type='l',
				main="sample nuggets",
				ylab="nugget value",
				xlab="MCMC iterations")
	return(invisible("nuggets"))
}

get.ylim <- function(cluster.params,n.clusters,param,z){
	y.lim <- range(unlist(
				lapply(
					lapply(1:n.clusters,
						function(i){
							cluster.params[[i]][[param]]
						}),
					function(x){
						range(x[z])
					})))
	y.lim <- y.lim + c(-0.15*diff(y.lim),0.15*diff(y.lim))
	return(y.lim)
}

plot.cluster.param <- function(cluster.param,clst.col,z){
	points(cluster.param[z],type='l',col=clst.col)
	return(invisible(0))
}

plot.cluster.cov.params <- function(data.block,geoStr.results,burnin,cluster.colors){
	#recover()
	n.clusters <- data.block$K
	n.iter <- length(geoStr.results$post$posterior)
	z <- (burnin+1):n.iter
	params <- names(geoStr.results$post$cluster.params$Cluster_1)[!names(geoStr.results$post$cluster.params$Cluster_1)=="cluster.cov"]
	param.ranges <- lapply(params,function(x){get.ylim(geoStr.results$post$cluster.params,n.clusters,x,z)})
	if(length(params) > 0){
		for(i in 1:length(params)){
			plot(0,type='n',main=params[i],
				xlab="MCMC iterations",ylab="parameter value",
				xlim=c(1,length(z)),
				ylim=param.ranges[[i]])
			lapply(1:n.clusters,function(j){plot.cluster.param(geoStr.results$post$cluster.params[[j]][[params[i]]],cluster.colors[j],z)})
			legend(x="topright",col= cluster.colors[1:n.clusters],lty=1,legend=paste0("Cluster_",1:n.clusters))
		}
	}
	return(invisible(0))
}

plot.admix.props <- function(data.block,geoStr.results,cluster.colors,burnin){
#	recover()
	n.clusters <- data.block$K
	n.iter <- length(geoStr.results$post$posterior)
	par(mfrow=c(n.clusters,1),mar=c(3,3,2,2))
		for(i in 1:n.clusters){
			matplot(geoStr.results$post$admix.proportions[(burnin+1):n.iter,,i],type='l',ylim=c(0,1),
					main=paste0("Cluster ",i),ylab="admixture proportion",col=cluster.colors[i])
		}
	return(invisible(0))
}

plot.model.fit <- function(data.block,std.freq.list,geoStr.results,burnin){
	n.iter <- length(geoStr.results$post$posterior)
	z <- (burnin+1):n.iter
	index.mat <- upper.tri(data.block$geoDist, diag = TRUE)
	cov.range <- range(c(std.freq.list$std.cov.list$std.cov,geoStr.results$post$mc.par.cov[z, , ]))
    plot(data.block$geoDist,std.freq.list$std.cov.list$std.cov,
    	xlab = "geographic distance", 
        ylab = "covariance",
        main="Cov/geoDist",
        ylim = cov.range, type = "n")
    lapply(z, function(i) {
        points(data.block$geoDist[index.mat], geoStr.results$post$mc.par.cov[i,,][index.mat],
        	pch = 20, col = adjustcolor(1, 0.1))
    		})
    points(data.block$geoDist[index.mat], std.freq.list$std.cov.list$std.cov[index.mat], 
        xlab = "geographic distance", ylab = "covariance", ylim = cov.range, 
        col=2,pch = 19)
	legend(x="topright",legend=c("observed","parametric"),pch=19,col=c(2,1))
	return(invisible("plotted"))
}

plot.cluster.covs <- function(data.block,geoStr.results,cluster.colors,burnin){
	order.mat <- order(data.block$geoDist)
	n.iter <- length(geoStr.results$post$posterior)
	z <- (burnin+1):n.iter
	y.range <- range(unlist(lapply(1:data.block$K,function(k){geoStr.results$post$cluster.params[[k]]$cluster.cov[z]})))
	plot(0,xlim=range(data.block$geoDist),ylim=y.range,type='n',
			xlab = "geographic distance",
			ylab = "cluster-specific covariances",
			main = "cluster covariances")
		lapply(1:data.block$K,function(k){
			lapply(geoStr.results$post$cluster.params[[k]]$cluster.cov[z],function(x){
				lines(data.block$geoDist[order.mat],
					  x[order.mat],col=adjustcolor(cluster.colors[k],0.1),pch=20)
			})
		})
	legend(x="topright",col= cluster.colors[1:data.block$K],lty=1,legend=paste0("Cluster_",1:data.block$K))
	return(invisible("cluster covs"))	
}

plot.cluster.covariances <- function(data.block,geoStr.results,cluster.colors){
	ind.mat <- upper.tri(data.block$geoDist,diag=TRUE)
	y.range <- range(unlist(lapply(seq_along(1:data.block$K),
										function(i){geoStr.results$point$cluster.params[[i]]$cluster.cov})))
	plot(data.block$geoDist[ind.mat],data.block$geoDist[ind.mat],
			ylab="sample covariance",xlab="geographic distance",
			ylim=y.range + diff(range(y.range))/10 * c(-1,1),type='n')
	for(i in 1:data.block$K){
		points(data.block$geoDist[ind.mat],
				geoStr.results$point$cluster.params[[i]]$cluster.cov[ind.mat],
				col=cluster.colors[i],pch=20,cex=0.6)
	}
	legend(x="topright",lwd=2,lty=1,col=cluster.colors[1:data.block$K],
			legend=paste0("Cluster_",1:data.block$K),cex=0.7)
	return(invisible("plotted"))
}

structure.polygon <- function(plotting.admix.props,i,j,use.colors){
	polygon(x = c(j-1,j,j,j-1),
			y = c(plotting.admix.props[i,j],
					plotting.admix.props[i,j],
					plotting.admix.props[i+1,j],
					plotting.admix.props[i+1,j]),
			col=use.colors[i])
	return(invisible(j))
}

make.structure.polygon.layer <- function(plotting.admix.props,i,use.colors,sample.order){
	# recover()
		lapply(1:ncol(plotting.admix.props),function(j){
			structure.polygon(plotting.admix.props[,sample.order],i,j,use.colors)
		})
	return(invisible(i))
}

make.structure.plot <- function(data.block,geoStr.results,mar=c(2,4,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=NULL){
	# recover()
	# quartz(width=10,height=5)
	par(mar=mar)
	if(is.null(cluster.order)){
		cluster.order <- seq(1:data.block$K)
	}
	if(is.null(sample.order)){
		sample.order <- seq(1:data.block$N)
	}
	if(!is.null(sort.by)){
		sample.order <- order(geoStr.results$point$admix.proportions[,sort.by])
	}
	if(is.null(cluster.colors)){
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	}
	if(data.block$K==1){
		geoStr.results$point$admix.proportions <- matrix(geoStr.results$point$admix.proportions,nrow=data.block$N,ncol=1)
	}
	use.colors <- cluster.colors[1:data.block$K][cluster.order]
	plot(0,xlim=c(0,data.block$N),ylim=c(0,1),type='n',ylab="admixture",xlab="",xaxt='n')
	plotting.admix.props <- apply(cbind(0,geoStr.results$point$admix.proportions[,cluster.order]),1,cumsum)
	lapply(1:data.block$K,function(i){
		make.structure.polygon.layer(plotting.admix.props,i,use.colors,sample.order)
	})
	if(!is.null(sample.names)){
		axis(side=1,at=seq(1:data.block$N)-0.5,labels=sample.names[sample.order],cex.axis=0.5,las=2)
	}
	return(invisible("plotted"))
}

make.admix.pie.plot <- function(data.block,geoStr.results,cluster.colors,stat,radii=2.7,add=FALSE,title=NULL,x.lim=NULL,y.lim=NULL){
	# recover()
	if(is.null(data.block$coords)){
		message("\nuser has not specified sampling coordinates in the data block\n")
	} else {
		require(caroline)
		cluster.names <- paste0("cluster_",1:data.block$K)
		sample.names <- paste0("sample_",1:data.block$N)
		color.tab <- nv(c(cluster.colors[1:data.block$K]),cluster.names)
		if(stat == "MAP"){
			admix.props <- geoStr.results$point$admix.proportions
		} else if(stat == "mean"){
			admix.props <- apply(geoStr.results$post$admix.proportions,c(2,3),mean)
		} else if(stat == "median"){
			admix.props <- apply(geoStr.results$post$admix.proportions,c(2,3),median)		
		}
		pie.list <- lapply(1:data.block$N,function(i){nv(admix.props[i,],cluster.names)})
		names(pie.list) <- sample.names
		if(add){
			par(new=TRUE)
		} else {
			par(mar=c(2,2,2,2))
		}
		if(is.null(title)){
			title <- "Admixture proportion map"
		}
		if(is.null(x.lim)){
			x.lim <- c(min(data.block$coords[,1]) - 1, max(data.block$coords[,1]) + 1)
		}
		if(is.null(y.lim)){
			y.lim <- c(min(data.block$coords[,2]) - 1, max(data.block$coords[,2]) + 1)
		}
		pies(pie.list,x0=data.block$coords[,1],y0=data.block$coords[,2],
					color.table=color.tab,border="black",radii=radii,
					xlab="",ylab="",main=title,lty=1,density=NULL,
					xlim = x.lim, ylim = y.lim)
		box(lwd=2)
	}
	return(invisible(0))
}

get.cluster.order <- function(K,admix.props,ref.admix.props){
	K.combn <- expand.grid(1:K,1:K)
	mean.props <- lapply(1:K,function(i){
						apply(admix.props[,,i],2,
							function(x){mean(x)})
					})	
	admix.prop.cors <- unlist(
						lapply(1:nrow(K.combn),function(i){
							cor(mean.props[[K.combn[i,1]]],
							ref.admix.props[,K.combn[i,2]])}))
	matchups <- matrix(c(1:K,rep(NA,K)),nrow=K,ncol=2)
		colnames(matchups) <- c("sample","reference")
	while(any(is.na(matchups))){
		winner <- as.numeric(K.combn[which.max(admix.prop.cors),])
		matchups[winner[1],2] <- winner[2]
		admix.prop.cors[which.max(admix.prop.cors)] <- NA
	}
	return(matchups)
}

make.all.chain.plots <- function(geoStr.results,chain.no,data.block,std.freq.list,prefix,burnin,cluster.colors,...){
	pdf(file=paste0(prefix,"_trace.plots.chain_",chain.no,".pdf"),...)
		plot.prob(geoStr.results,burnin)
		plot.nuggets(geoStr.results,burnin)
		plot.cluster.cov.params(data.block,geoStr.results,burnin,cluster.colors)
		if(data.block$K > 1){
			plot.admix.props(data.block,geoStr.results,cluster.colors,burnin)
		}
	dev.off()
	pdf(file=paste0(prefix,"_model.fit.chain_",chain.no,".pdf"),...)
		plot.model.fit(data.block,std.freq.list,geoStr.results,burnin)
	dev.off()
	if(data.block$spatial | data.block$K > 1){
		pdf(file=paste0(prefix,"_cluster.cov.curves.chain_",chain.no,".pdf"),width=5,height=5)
			plot.cluster.covs(data.block,geoStr.results,cluster.colors,burnin)
		dev.off()
	}
	if(data.block$K > 1){
		pdf(file=paste0(prefix,"_pie.map.chain_",chain.no,".pdf"),width=6,height=6)	
			make.admix.pie.plot(data.block,geoStr.results,cluster.colors,stat="median",radii=2.7,add=FALSE,title=NULL,x.lim=NULL,y.lim=NULL)
		dev.off()
		pdf(file=paste0(prefix,"_structure.plot.chain_",chain.no,".pdf"),width=10,height=5)
			make.structure.plot(data.block,geoStr.results,mar=c(2,4,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=cluster.colors)
		dev.off()
	}
	return(invisible("made chain plots!"))
}

make.all.the.plots <- function(geoStr.results,n.chains,data.block,std.freq.list,prefix,burnin,cluster.colors=NULL,...){
	if(is.null(cluster.colors)){
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	}
	lapply(1:n.chains,function(i){
		make.all.chain.plots(geoStr.results[[i]],chain.no=i,data.block,std.freq.list,prefix,burnin,cluster.colors,...)
	})
	return(invisible("made chain plots!"))
}

random.switcharoo <- function(x){
	x <- ifelse(rep(runif(1) < 0.5,length(x)),
					x,
					1-x)
	return(x)
}

switcharoo.data <- function(frequencies){
	frequencies <- apply(frequencies,2,random.switcharoo)
	return(frequencies)
}