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

make.structure.plot <- function(data.list,super.list,sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=NULL){
	# recover()
	quartz(width=10,height=5)
	if(is.null(cluster.order)){
		cluster.order <- seq(1:ncol(super.list$parameter.list$admix.proportions))
	}
	if(is.null(sample.order)){
		sample.order <- seq(1:data.list$n.ind)
	}
	if(!is.null(sort.by)){
		sample.order <- order(super.list$parameter.list$admix.proportions[,sort.by])
	}
	if(is.null(cluster.colors)){
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen")
	}
	all.colors <- cluster.colors
	use.colors <- all.colors[1:length(super.list$parameter.list$cluster.list)][cluster.order]
	plot(0,xlim=c(0,data.list$n.ind),ylim=c(0,1),type='n',ylab="admixture",xlab="",xaxt='n')
	plotting.admix.props <- apply(cbind(0,super.list$parameter.list$admix.proportions[,cluster.order]),1,cumsum)
	lapply(1:length(super.list$parameter.list$cluster.list),function(i){
		make.structure.polygon.layer(plotting.admix.props,i,use.colors,sample.order)
	})
	if(!is.null(sample.names)){
		axis(side=1,at=seq(1:data.list$n.ind)-0.5,labels=sample.names[sample.order],cex.axis=0.5,las=2)
	}
	return(invisible("plotted"))
}

make.spatiotemporal.sampling.plot <- function(metadata,aspect.vec=NULL){
	require(rgl)
	require(maps)
	world.map.lines <- map(database="world",plot=FALSE,xlim=range(metadata$lon),ylim=c(range(metadata$lat)))
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
	plot3d(x=metadata$lon,y=metadata$lat,
			z=metadata$time,type='n',xlab="",ylab="",zlab="",axes=FALSE)
		planes3d(0,0,-1,d=c(seq(min(metadata$time),max(metadata$time),length.out=5)),alpha=0.1,col="gray")
			lines3d(world.map.lines$x,world.map.lines$y,z=min(metadata$time),col=1,add=TRUE)
			lines3d(world.map.lines$x,world.map.lines$y,z=max(metadata$time),col=1,add=TRUE)
		box3d(lwd=2)
				mtext3d("x",text="Longitude",
				pos=c(mean(plotting.x.range),
						plotting.y.range[1]-diff(plotting.y.range)/7,
						min(metadata$time)),cex=1.5)
		mtext3d("y",text="Latitude",
				pos=c(plotting.x.range[1]-diff(plotting.x.range)/7,
						mean(plotting.y.range),
						min(metadata$time)),cex=1.5)
		mtext3d("z",text="Sampling year",
				pos=c(plotting.x.range[1]-diff(plotting.y.range)/7,
						plotting.y.range[2]+diff(plotting.y.range)/7,
						mean(metadata$time)),cex=1.5)
		axis3d("x++")
		axis3d("y++")
		axis3d("z+")
				#plot3d(x=metadata$lon,y=metadata$lat,z=metadata$time,
		#		type='s',size=0.7,shininess=128,add=TRUE)


		# rgl.snapshot(filename=paste(taxon.name,"_3D_1",".png",sep=""))
		# par3d(userMatrix=view.matrix2)
		# rgl.snapshot(filename=paste(taxon.name,"_3D_2",".png",sep=""))
		# par3d(userMatrix=view.matrix1)
		#play3d(spin3d(rpm=3),duration=5)
		# return(0)
}

pies3d <- function (x, n.clusters, locations, radii, edges=100) {
    	# recover()
    if (length(radii) < length(x)){ 
        radii <- rep(radii, length.out = length(x))
	}
	all.colors <- c("blue","red","green","yellow","purple","brown")
	use.colors <- all.colors[1:n.clusters]
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
				polygon3d(x = c(locations[j,]$lon + P$x,locations[j,]$lon),
							y = c(locations[j,]$lat + P$y,locations[j,]$lat),
							z = rep(locations[j,]$time,fineness+1),
							col=use.colors[i])
			}
		}
	}
	return(invisible(0))
}

pies3d.2 <- function (x, n.clusters, locations, radii, edges=100) {
    	# recover()
    if (length(radii) < length(x)){ 
        radii <- rep(radii, length.out = length(x))
	}
	all.colors <- c("blue","red","green","yellow","purple","brown")
	use.colors <- all.colors[1:n.clusters]
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
					xv = c(locations[j,]$lon + P$x,locations[j,]$lon)
					yv = c(locations[j,]$lat + P$y,locations[j,]$lat)
				slice <- translate3d(extrude3d(xv,yv,thickness=100),0,0,locations$time[j])
				shade3d(slice,col=use.colors[i])
			}
		}
	}
	return(invisible(0))
}

plot.sample.ellipses <- function(metadata,scale,alpha=0.5){
	cx <- 0.045 * par3d()$scale[1]
	cy <- 0.019 * par3d()$scale[2]
		ellipse.matrix <- matrix(scale*c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3,byrow=TRUE)
	lapply(1:nrow(metadata),function(i){
			plot3d(ellipse3d(ellipse.matrix, scale = c(cx,cy,1), centre=c(metadata[i,]$lon,metadata[i,]$lat,metadata[i,]$time)), col=1, lty=2, alpha=alpha, add = TRUE) ; 
			return(invisible(i))
		})
	return(invisible("done"))
}

plot.sample.ellipses2 <- function(metadata,radii,alpha=0.5,fineness=100,thickness=1,col=1){
	save.shiny <- material3d("shininess")
		material3d("shininess"=30)
	cx <- 0.045 * par3d()$scale[1]
	cy <- 0.019 * par3d()$scale[2]
	angles <- seq(0,2*pi,length.out=fineness)
	n.samples <- nrow(metadata)
	if(length(radii) < n.samples){
		radii <- rep(radii,n.samples)
	}
	if(length(col) < n.samples){
		col <- rep(col,n.samples)
	}
    for (j in seq_along(1:n.samples)) {
		xv <- cx * cos(angles)*radii[j]
		yv <- cy * sin(angles)*radii[j]
		cyl <- translate3d(extrude3d(xv,yv,thickness=100),metadata$lon[j],metadata$lat[j],metadata$time[j])
		shade3d(cyl,col=col[j])
	}
	material3d("shininess"=save.shiny)
	return(invisible("done"))
}

plot.cluster.covariances <- function(data.list,super.list,time.dist=FALSE){
	dist.mat <- ifelse(matrix(time.dist,nrow=data.list$n.ind,ncol=data.list$n.ind),
						data.list$time.dist,
						data.list$geo.dist)
	cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen")
	ind.mat <- upper.tri(dist.mat,diag=TRUE)
	y.min <- min(data.list$sample.covariance,
						unlist(lapply(seq_along(1:super.list$model.options$n.clusters),
										function(i){super.list$parameter.list$cluster.list[[i]]$covariance + 
													super.list$parameter.list$cluster.list[[i]]$cluster.mean})))
	y.max <- max(data.list$sample.covariance,
						unlist(lapply(seq_along(1:super.list$model.options$n.clusters),
										function(i){super.list$parameter.list$cluster.list[[i]]$covariance + 
													super.list$parameter.list$cluster.list[[i]]$cluster.mean})))
	y.range <- c(y.min,y.max)
	plot(dist.mat[ind.mat],data.list$sample.covariance[ind.mat],
			ylab="sample covariance",xlab="distance",
			ylim=y.range + diff(range(y.range))/10 * c(-1,1))
	abline(h=super.list$parameter.list$shared.mean,col="gray",lty=2)
	for(i in 1:super.list$model.options$n.clusters){
		points(dist.mat[ind.mat],
				super.list$parameter.list$cluster.list[[i]]$covariance[ind.mat] + 
				super.list$parameter.list$cluster.list[[i]]$cluster.mean,
				col=cluster.colors[i],pch=20,cex=0.7)
	}
	return(invisible("plotted"))
}