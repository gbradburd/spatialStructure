get.best.iter <- function(model.fit,chain.no){
	logpost <- get_logposterior(model.fit)
	best.iter <- lapply(logpost,which.max)[[chain.no]]
	return(best.iter)
}

get.admix.props <- function(model.fit,chain.no,N,n.clusters){
	# recover()
	admix.props <- array(1,dim=c(model.fit@stan_args[[chain.no]]$iter,N,n.clusters))
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
	par.cov <- array(NA,dim=c(model.fit@stan_args[[chain.no]]$iter,N,N))
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

get.gamma <- function(model.fit,chain.no){
	gamma <- extract(model.fit,pars="gamma",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(gamma)
}

get.alpha.params <- function(model.fit,chain.no,cluster,n.clusters){
	alpha.pars <- model.fit@model_pars[grepl("alpha",model.fit@model_pars)]
	if(n.clusters > 1){
		alpha.params <- setNames(
							lapply(1:length(alpha.pars),
									function(i){
										extract(model.fit,
												pars=paste(alpha.pars[i],"[",cluster,"]",sep=""),
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
					pars=paste("mu","[",cluster,"]",sep=""),
					inc_warmup=TRUE,permute=FALSE)[,chain.no,]
	return(mu)
}

get.cov.function <- function(data.block){
	if(data.block$K == 1){
		if(data.block$space & !data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(data.block$geoDist*cluster.params$alphaD)^cluster.params$alpha2))
			}
		}
		if(!data.block$space & data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(cluster.params$alphaT*data.block$timeDist)^cluster.params$alpha2))
			}	
		}
		if(data.block$space & data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(data.block$geoDist*cluster.params$alphaD + 
							cluster.params$alphaT*data.block$timeDist)^cluster.params$alpha2))
			}	
		}
		if(!data.block$space & !data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(0)
			}
		}
	} else {
		if(data.block$space & !data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(data.block$geoDist*cluster.params$alphaD)^cluster.params$alpha2) + 
						cluster.params$mu)
			}
		}
		if(!data.block$space & data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(cluster.params$alphaT*data.block$timeDist)^cluster.params$alpha2) + 
						cluster.params$mu)
			}	
		}
		if(data.block$space & data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(exp(log(cluster.params$alpha0) + 
						-(data.block$geoDist*cluster.params$alphaD + 
							cluster.params$alphaT*data.block$timeDist)^cluster.params$alpha2) + 
						cluster.params$mu)
			}	
		}
		if(!data.block$space & !data.block$time){
			cov.func <- function(cluster.params,data.block){
				return(cluster.params$mu)
			}
		}
	}
	return(cov.func)
}

get.cluster.cov <- function(cluster.params,data.block){
	cov.function <- get.cov.function(data.block)
	cluster.cov <- lapply(1:length(cluster.params[[1]]),
							function(i){
								cov.function(cluster.params=lapply(cluster.params,"[[",i),data.block)
							})
	return(cluster.cov)
}

get.cluster.params <- function(model.fit,data.block,chain.no,cluster,n.clusters){
	cluster.params <- list()
	if(data.block$space | data.block$time){
		cluster.params <- get.alpha.params(model.fit,chain.no,cluster,n.clusters)
	}
	if(n.clusters > 1){
		cluster.params <- c(cluster.params,
							list("mu" = get.cluster.mu(model.fit,chain.no,cluster)))
	}
	cluster.cov <- get.cluster.cov(cluster.params,data.block)
	cluster.params <- c(cluster.params,list("cluster.cov"=cluster.cov))
	return(cluster.params)
}

get.cluster.params.list <- function(model.fit,data.block,chain.no){
	cluster.params <- setNames(
						lapply(1:data.block$K,
									function(i){
										get.cluster.params(model.fit,data.block,chain.no,i,data.block$K)
									}),
						paste("Cluster",1:data.block$K,sep="_"))
	return(cluster.params)
}

index.best <- function(param,best){
	if(class(param) == "numeric"){
		best.param <- param[best]
	}
	if(class(param) == "list"){
		best.param <- param[[best]]
	}
	if(class(param) == "array"){
		best.param <- param[best,,]
	}
	if(class(param) == "matrix"){
		best.param <- param[best,]
	}
	return(best.param)
}

index.best.cluster.params <- function(cluster.params,best){
	best.cluster.params <- lapply(cluster.params,index.best,best)
	return(best.cluster.params)
}

index.best.cluster.params.list <- function(cluster.params.list,best){
	best.cluster.params.list <- lapply(cluster.params.list,index.best.cluster.params,best)
	return(best.cluster.params.list)
}

get.geoStructure.results <- function(model.fit,chain.no,data.block){
	require(rstan)
	load(model.fit)
	load(data.block)
	post <- list("posterior" = get_logposterior(model.fit)[[chain.no]],
				 "admix.proportions" = get.admix.props(model.fit,chain.no,data.block$N,data.block$K),
				 "nuggets" = get.nuggets(model.fit,chain.no,data.block$N),
				 "cluster.params" = get.cluster.params.list(model.fit,data.block,chain.no),
				 "gamma" = get.gamma(model.fit,chain.no)	,
				 "par.cov" = get.par.cov(model.fit,chain.no,data.block$N))
	best <- get.best.iter(model.fit,chain.no)
	point <- list("admix.proportions" = index.best(post$admix.proportions,best),
				 "nuggets" = index.best(post$nuggets,best),
				 "cluster.params" = index.best.cluster.params.list(post$cluster.params,best),
				 "gamma" = index.best(post$gamma,best),
				 "par.cov" = index.best(post$par.cov,best))
	geoStructure.results <- list("post" = post,"point" = point)
	return(geoStructure.results)
}

plot.cluster.covariances <- function(data.block,geoStr.results,time,cluster.colors){
	ind.mat <- upper.tri(data.block$geoDist,diag=TRUE)
	y.range <- range(data.block$obsSigma,
						unlist(lapply(seq_along(1:data.block$K),
										function(i){geoStr.results$point$cluster.params[[i]]$cluster.cov + 
													geoStr.results$point$cluster.params[[i]]$mu})))
	n.col <- ifelse(time,2,1)
	par(mfrow=c(1,n.col))
	plot(data.block$geoDist[ind.mat],data.block$obsSigma[ind.mat],
			ylab="sample covariance",xlab="geographic distance",
			ylim=y.range + diff(range(y.range))/10 * c(-1,1))
	abline(h=geoStr.results$point$mu,col="gray",lty=2)
	for(i in 1:data.block$K){
		points(data.block$geoDist[ind.mat],
				geoStr.results$point$cluster.params[[i]]$cluster.cov[ind.mat] + 
				geoStr.results$point$cluster.params[[i]]$mu,
				col=cluster.colors[i],pch=20,cex=0.6)
	}
	legend(x="topright",lwd=2,lty=1,col=c("black",cluster.colors[1:data.block$K]),
			legend=c("Observed Cov",paste0("Cluster_",1:data.block$K)),cex=0.7)
	if(time){
		plot(data.block$timeDist[ind.mat],data.block$obsSigma[ind.mat],
			ylab="sample covariance",xlab="temporal distance",
			ylim=y.range + diff(range(y.range))/10 * c(-1,1))
		abline(h=geoStr.results$point$mu,col="gray",lty=2)
		for(i in 1:data.block$K){
			points(data.block$timeDist[ind.mat],
				geoStr.results$point$cluster.params[[i]]$cluster.cov[ind.mat] + 
				geoStr.results$point$cluster.params[[i]]$mu,
					col=cluster.colors[i],pch=20,cex=0.7)
		}
	}
	return(invisible("plotted"))
}

plot.model.fit <- function(data.block,geoStr.results,time){
	n.col <- ifelse(time,3,2)
	par(mfrow=c(1,n.col),mar=c(4.5,4.5,1,1))
		plot(data.block$geoDist,data.block$obsSigma,
			xlab="geographic distance",
			ylab="covariance")
		points(data.block$geoDist,geoStr.results$point$par.cov,pch=20,col="red",cex=0.8)
		legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"),cex=0.7)
		if(time){
			plot(data.block$timeDist,data.block$obsSigma,
				xlab="temporal distance",
				ylab="covariance")
			points(data.block$timeDist,geoStr.results$point$par.cov,pch=20,col="red")
		}
		plot(data.block$obsSigma,geoStr.results$point$par.cov,
				xlab="sample covariance",
				ylab ="parametric admixed covariance")
			abline(0,1,col="red")
	return(invisible(0))
}

plot.lnl <- function(geoStr.results,burnin=3000){
	mcmc.length <- length(geoStr.results$post$posterior)
	plot.vec <- (burnin+1):mcmc.length
	plot(geoStr.results$post$posterior[plot.vec],xlab="",ylab="",type='l')
	return(invisible(0))
}

get.ylim <- function(cluster.params,param,n.iter,burnin=0){
	y.lim <- range(unlist(lapply(
					lapply(cluster.params,"[[",param),
				"[",(burnin+1):n.iter)))
	return(y.lim)
}

plot.cluster.param <- function(cluster.param,n.iter,clst.col,burnin=0){
	points(cluster.param[(burnin+1):n.iter],type='l',col=clst.col)
	return(invisible(0))
}

plot.cluster.cov.params <- function(data.block,geoStr.results,time,burnin,cluster.colors){
#	recover()
	n.col <- ifelse(time,5,4)
	n.clusters <- data.block$K
	n.iter <- length(geoStr.results$post$posterior)
	params <- names(geoStr.results$post$cluster.params$Cluster_1)[-which(names(geoStr.results$post$cluster.params$Cluster_1)=="cluster.cov")]
	param.ranges <- lapply(params,function(x){get.ylim(geoStr.results$post$cluster.params,x,burnin)})
	par(mfrow=c(1,n.col),mar=c(1,2,2,1))
		for(i in 1:length(params)){
			plot(0,type='n',main=params[i],
				xlab="MCMC iterations",ylab="parameter value",
				xlim=c(1,n.iter-burnin),
				ylim=param.ranges[[i]])
			lapply(1:n.clusters,function(j){plot.cluster.param(geoStr.results$post$cluster.params[[j]][[params[i]]],n.iter,cluster.colors[j],burnin)})
		}
		legend(x="topright",col= cluster.colors[1:n.clusters],lty=1,legend=paste0("Cluster_",1:n.clusters))
	return(invisible(0))
}

plot.admix.props <- function(data.block,geoStr.results,cluster.colors){
#	recover()
	n.clusters <- data.block$K
	par(mfrow=c(1,n.clusters),mar=c(2,2,2,2))
		for(i in 1:n.clusters){
			matplot(geoStr.results$post$admix.proportions[,,i],type='l',ylim=c(0,1),
					main=paste("Cluster ",i,sep=""),ylab="admixture proportion",col=cluster.colors[i])
		}
	return(invisible(0))
}

plot.nuggets <- function(geoStr.results,burnin){
	par(mar=c(2,2,2,2))
	n.iter <- length(geoStr.results$post$posterior)
	matplot(geoStr.results$post$nuggets,type='l',
				ylim=range(geoStr.results$post$nuggets[(burnin+1):n.iter,]),
				main="sample nuggets",
				ylab="nugget value")
	return(invisible(0))
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

make.structure.plot <- function(data.block,geoStr.results,mar=c(2,2,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=NULL){
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

make.admix.pie.plot <- function(data.block,geoStr.results,cluster.colors,cluster.names,radii=2.7,add=FALSE,title=NULL,xlim=NULL,ylim=NULL){
	# recover()
	require(caroline)
	sample.names <- paste0("sample_",1:data.block$N)
	color.tab <- nv(c(cluster.colors[1:data.block$K]),cluster.names)
	pie.list <- lapply(1:data.block$N,function(i){nv(geoStr.results$point$admix.proportions[i,],cluster.names)})
	names(pie.list) <- sample.names
	if(add){
		par(new=TRUE)
	} else {
		par(mar=c(2,2,2,2))
	}
	if(is.null(title)){
		title <- "admixture proportion map"
	}
	if(is.null(xlim)){
		xlim <- c(min(data.block$geoCoords[,1]) - 1, max(data.block$geoCoords[,1]) + 1)
	}
	if(is.null(xlim)){
		y.lim <- c(min(data.block$geoCoords[,2]) - 1, max(data.block$geoCoords[,2]) + 1)
	}
	pies(pie.list,x0=data.block$geoCoords[,1],y0=data.block$geoCoords[,2],
				color.table=color.tab,border="black",radii=radii,
				xlab="",ylab="",main=title,lty=1,density=NULL,
				xlim = xlim, ylim = ylim)
	box(lwd=2)
	return(invisible(0))
}

make.spatiotemporal.sampling.box <- function(data.block,aspect.vec=NULL){
	require(rgl)
	require(maps)
	world.map.lines <- map(database="world",plot=FALSE,xlim=range(data.block$geoCoords[,1]),ylim=c(range(data.block$geoCoords[,2])))
	plotting.y.range <- range(world.map.lines$y[which(!is.na(world.map.lines$y))])
	plotting.x.range <- range(world.map.lines$x[which(!is.na(world.map.lines$x))])
	view.matrix1 <- matrix(c(0.6739255,-0.7387993,0.0004163443,0,
								0.1659686,0.1519439,0.9743548632,0,
								-0.7199159,-0.6565735,0.2250164002,0,
								0,0,0,1),nrow=4,ncol=4,byrow=TRUE)
	view.matrix2 <- matrix(c(0.999946296,-0.007756224,-0.006859419,0,
								0.006820573,0.991850674,-0.127225429,0,
								0.007790213,0.127171844,0.991850257,0,
								0,0,0,1),nrow=4,ncol=4,byrow=TRUE)
	windowRect <- c(100,100,900,900)
	open3d(windowRect=windowRect,zoom=1.1,userMatrix=view.matrix1)
	plot3d(x=data.block$geoCoords[,1],y=data.block$geoCoords[,2],
			z=data.block$timeCoords,type='n',xlab="",ylab="",zlab="",axes=FALSE)
		planes3d(0,0,-1,d=c(seq(min(data.block$timeCoords),max(data.block$timeCoords),length.out=5)),alpha=0.1,col="gray")
			lines3d(world.map.lines$x,world.map.lines$y,z=min(data.block$timeCoords),col=1,add=TRUE)
			lines3d(world.map.lines$x,world.map.lines$y,z=max(data.block$timeCoords),col=1,add=TRUE)
		box3d(lwd=2)
				mtext3d("x",text="Longitude",
				pos=c(mean(plotting.x.range),
						plotting.y.range[1]-diff(plotting.y.range)/7,
						min(data.block$timeCoords)),cex=1.5)
		mtext3d("y",text="Latitude",
				pos=c(plotting.x.range[1]-diff(plotting.x.range)/7,
						mean(plotting.y.range),
						min(data.block$timeCoords)),cex=1.5)
		mtext3d("z",text="Sampling year",
				pos=c(plotting.x.range[1]-diff(plotting.y.range)/7,
						plotting.y.range[2]+diff(plotting.y.range)/7,
						mean(data.block$timeCoords)),cex=1.5)
		axis3d("x++")
		axis3d("y++")
		axis3d("z+")
		return(invisible("3D plotbox"))
}

plot.sample.ellipses2 <- function(data.block,radii,alpha=0.5,fineness=100,thickness=1,col=1){
	save.shiny <- material3d("shininess")
		material3d("shininess"=30)
	cx <- 0.045 * par3d()$scale[1]
	cy <- 0.019 * par3d()$scale[2]
	angles <- seq(0,2*pi,length.out=fineness)
	n.samples <- nrow(data.block$geoCoords)
	if(length(radii) < n.samples){
		radii <- rep(radii,n.samples)
	}
	if(length(col) < n.samples){
		col <- rep(col,n.samples)
	}
    for (j in seq_along(1:n.samples)) {
		xv <- cx * cos(angles)*radii[j]
		yv <- cy * sin(angles)*radii[j]
		cyl <- translate3d(extrude3d(xv,yv,thickness=100),
							data.block$geoCoords[j,1],
							data.block$geoCoords[j,2],
							data.block$timeCoords[j])
		shade3d(cyl,col=col[j])
	}
	material3d("shininess"=save.shiny)
	return(invisible("done"))
}

pies3d.2 <- function (x, n.clusters, cluster.colors,data.block, radii, edges=100) {
    	# recover()
    if (length(radii) < length(x)){ 
        radii <- rep(radii, length.out = length(x))
	}
	use.colors <- cluster.colors[1:n.clusters]
	cx <- 0.045 * par3d()$scale[1]
	cy <- 0.019 * par3d()$scale[2]
    for (j in seq(along = x)) {
        X <- x[[j]]
		X <- c(0, cumsum(X)/sum(X))
		dx <- diff(X)
		nx <- length(dx)
		twopi <- 2 * pi
		t2xy <- function(t) {
			t2p <- twopi * t
			return(list(x = radii[j] * cx * cos(t2p), y = radii[j] * cy * sin(t2p)))
		}
		for (i in 1:nx) {
			fineness <- max(2, floor(edges * dx[i]))
			P <- t2xy(seq.int(X[i], X[i + 1], length.out = fineness))
			if(fineness != 2){
					xv = c(data.block$geoCoords[j,1] + P$x,data.block$geoCoords[j,1])
					yv = c(data.block$geoCoords[j,2] + P$y,data.block$geoCoords[j,2])
				slice <- translate3d(extrude3d(xv,yv,thickness=100),0,0,data.block$timeCoords[j])
				shade3d(slice,col=use.colors[i])
			}
		}
	}
	return(invisible(0))
}

make.3D.pie.plot <- function(data.block,radii=1,cluster.colors,admix.props=NULL,aspect.vec=NULL){
	require(caroline)
	make.spatiotemporal.sampling.box(data.block,aspect.vec)
	if(is.null(admix.props)){
		plot.sample.ellipses2(data.block,radii,alpha=0.5,fineness=100,thickness=1,col=1)
	} else {
		sample.names <- unlist(lapply(1:data.block$N,function(i){paste("sample_",i,sep="")}))
		cluster.names <- unlist(lapply(1:data.block$K,function(i){paste("Cluster_",i,sep="")}))
		color.tab <- nv(c(cluster.colors[1:data.block$K]),cluster.names)
		pie.list <- lapply(1:data.block$N,function(i){nv(admix.props[i,],cluster.names)})
		pies3d.2(x=pie.list,n.clusters=data.block$K,cluster.colors=cluster.colors,data.block=data.block,radii=radii,edges=100)
	}
}


make.all.the.plots <- function(dir,output.dir,burnin=0){
	setwd(dir)
	model.fit.file <- list.files(pattern="model.fit")
	data.block.file <- list.files(pattern="data.block")
	geoStr.results <- get.geoStructure.results(model.fit.file,chain.no=1,data.block.file)
	load(data.block.file)
	K <- data.block$K
	time <- data.block$time
	cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	cluster.names <- paste0("Cluster_",1:K)
	pdf(file=paste(output.dir,"/","model.fit.",K,".pdf",sep=""),width=(8+time*4),height=4,pointsize=14)
		plot.model.fit(data.block,geoStr.results,time)
	dev.off()
	pdf(file=paste(output.dir,"/","cluster.cov.curves.",K,".pdf",sep=""),width=(5+time*5),height=5,pointsize=18)
		plot.cluster.covariances(data.block,geoStr.results,time,cluster.colors)
	dev.off()
	pdf(file=paste(output.dir,"/","lnl.and.prob.",K,".pdf",sep=""),width=6,height=5,pointsize=18)
		plot.lnl(geoStr.results,burnin)
	dev.off()
	pdf(file=paste(output.dir,"/","cluster.cov.params.",K,".pdf",sep=""),width=12,height=4,pointsize=18)
		plot.cluster.cov.params(data.block,geoStr.results,time,burnin,cluster.colors)
	dev.off()
	pdf(file=paste(output.dir,"/","admix.props.",K,".pdf",sep=""),width=(4*data.block$K),height=4,pointsize=18)
		plot.admix.props(data.block,geoStr.results,cluster.colors)
	dev.off()
	pdf(file=paste(output.dir,"/","nuggets.",K,".pdf",sep=""),width=6,height=4,pointsize=18)
		plot.nuggets(geoStr.results,burnin)
	dev.off()
	pdf(file=paste(output.dir,"/","pie.chart.map.",K,".pdf",sep=""),width=6,height=6,pointsize=18)	
		make.admix.pie.plot(data.block,geoStr.results,cluster.colors,cluster.names,radii=2.6,add=FALSE,title=NULL,xlim=NULL,ylim=NULL)
	dev.off()
	pdf(file=paste(output.dir,"/","structure.plot.",K,".pdf",sep=""),width=10,height=5,pointsize=18)
		make.structure.plot(data.block,geoStr.results,mar=c(2,2,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=cluster.colors)
	dev.off()
}