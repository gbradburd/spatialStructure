geoStructure.ML <- function(data.block,prefix,...){
	#validate on data block
	data.block <- validate.data.block(data.block)
		save(data.block,file=paste(prefix,"data.block.Robj",sep="_"))
	#validate on model specification
	#make.stan.code.block
	stan.block <- make.stan.code.block(data.block$space,data.block$time,data.block$K)
		#write stan block to file
	#run model
	#put stan in tryCatch, email me
	require(rstan)
	model.fit <- optimizing(object = stan_model(model_code= stan.block),
							data = data.block,
							...)
	#save fit obj
	save(model.fit,file=paste(prefix,"model.fit.Robj",sep="_"))
	#return fit obj
	return(model.fit)
}

get.geoStructure.results <- function(model.fit,data.block,chain.no=NULL){
	require(rstan)
	load(model.fit)
	load(data.block)
	if(class(model.fit)=="stanfit"){
		geoStructure.results <- get.geoStructure.Bayes.results(model.fit,chain.no,data.block)
	} else {
		geoStructure.results  <- get.geoStructure.ML.results(model.fit,data.block)
	}
	return(geoStructure.results)
}

get.ML.par.cov <- function(model.fit,data.block){
	N <- data.block$N
	par.cov <- matrix(NA,nrow=N,ncol=N)
	for(i in 1:N){
		for(j in 1:N){
			my.par <- sprintf("Sigma[%s,%s]",i,j)
			par.cov[i,j] <- model.fit$par[my.par]
		}
	}
	return(par.cov)
}

get.ML.nuggets <- function(model.fit,data.block){
	N <- data.block$N
	nuggets <- unlist(lapply(1:N,
						function(i){
							model.fit$par[paste0("nugget[",i,"]")]
						}))
	return(nuggets)
}

get.ML.admix.props <- function(model.fit,data.block){
	N <- data.block$N
	K <- data.block$K
	admix.props <- matrix(1,nrow=N,ncol=K)
	if(any(grepl("w",names(model.fit$par)))){
		for(k in 1:K){
			admix.props[,k] <- unlist(lapply(1:N,function(i){model.fit$par[paste0("w[",i,",",k,"]")]}))
		}
	}
	return(admix.props)
}

get.ML.alpha.params <- function(model.fit,cluster,n.clusters){
	alpha.pars <- unique(gsub("\\[.+\\]","",
							names(model.fit$par)[
								grepl("alpha",names(model.fit$par))
							]))
	if(n.clusters > 1){
		alpha.params <- setNames(
							lapply(1:length(alpha.pars),
									function(i){
										model.fit$par[paste0(alpha.pars[i],"[",cluster,"]")]
									}),alpha.pars)
	} else {
		alpha.params <- setNames(
							lapply(1:length(alpha.pars),
									function(i){
										model.fit$par[alpha.pars[i]]
									}),alpha.pars)		
	}
	return(alpha.params)
}

get.ML.cluster.mu <- function(model.fit,cluster){
	mu <- model.fit$par[paste0("mu","[",cluster,"]")]
	return(mu)
}

get.ML.cluster.params <- function(model.fit,data.block,cluster,n.clusters){
	cluster.params <- list()
	if(data.block$space | data.block$time){
		cluster.params <- get.ML.alpha.params(model.fit,cluster,n.clusters)
	}
	if(n.clusters > 1){
		cluster.params <- c(cluster.params,
							list("mu" = get.ML.cluster.mu(model.fit,cluster)))
	}
	cluster.cov <- get.cluster.cov(cluster.params,data.block)
	cluster.params <- c(cluster.params,"cluster.cov"=cluster.cov)
	return(cluster.params)
}

get.ML.cluster.params.list <- function(model.fit,data.block){
	cluster.params <- setNames(
						lapply(1:data.block$K,
									function(i){
										get.ML.cluster.params(model.fit,data.block,i,data.block$K)
									}),
						paste("Cluster",1:data.block$K,sep="_"))
	return(cluster.params)
}

get.geoStructure.ML.results <- function(model.fit,data.block){
	point <- list("admix.proportions" = get.ML.admix.props(model.fit,data.block),
				  "nuggets" = get.ML.nuggets(model.fit,data.block),
				  "cluster.params" = get.ML.cluster.params.list(model.fit,data.block),
				  "gamma" = model.fit$par["gamma"],
				  "zeta" = model.fit$par["zeta"],
				  "par.cov" = get.ML.par.cov(model.fit,data.block))
	geoStructure.results <- list("post"=NULL,"point"=point)
	return(geoStructure.results)
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
		sample.names <- unlist(lapply(1:data.block$N,function(i){paste0("sample_",i)}))
		cluster.names <- unlist(lapply(1:data.block$K,function(i){paste0("Cluster_",i)}))
		color.tab <- nv(c(cluster.colors[1:data.block$K]),cluster.names)
		pie.list <- lapply(1:data.block$N,function(i){nv(admix.props[i,],cluster.names)})
		pies3d.2(x=pie.list,n.clusters=data.block$K,cluster.colors=cluster.colors,data.block=data.block,radii=radii,edges=100)
	}
}
