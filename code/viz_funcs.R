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
	# quartz(width=10,height=5)
	par(mar=c(2,2,2,2))
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
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
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

make.admix.pie.plot <- function(super.list,data.list,all.colors,cluster.names,radii=2.7,add=FALSE,title=NULL,xlim=NULL,ylim=NULL){
	require(caroline)
	sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
	color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
	pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
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
		xlim <- c(min(data.list$geo.coords[,1]) - 1, max(data.list$geo.coords[,1]) + 1)
	}
	if(is.null(xlim)){
		y.lim <- c(min(data.list$geo.coords[,2]) - 1, max(data.list$geo.coords[,2]) + 1)
	}
	pies(pie.list,x0=data.list$geo.coords[,1],y0=data.list$geo.coords[,2],
				color.table=color.tab,border="black",radii=2.7,
				xlab="",ylab="",main=title,lty=1,density=NULL,
				xlim = xlim, ylim = ylim)
	box(lwd=2)
	return(invisible(0))
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
	all.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
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
	all.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
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

plot.cluster.covariances <- function(data.list,super.list,time){
	# dist.mat <- ifelse(matrix(time.dist,nrow=data.list$n.ind,ncol=data.list$n.ind),
						# data.list$time.dist,
						# data.list$geo.dist)
	cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen")
	ind.mat <- upper.tri(data.list$geo.dist,diag=TRUE)
	y.min <- min(data.list$sample.covariance,
						unlist(lapply(seq_along(1:super.list$model.options$n.clusters),
										function(i){super.list$parameter.list$cluster.list[[i]]$covariance + 
													super.list$parameter.list$cluster.list[[i]]$cluster.mean})))
	y.max <- max(data.list$sample.covariance,
						unlist(lapply(seq_along(1:super.list$model.options$n.clusters),
										function(i){super.list$parameter.list$cluster.list[[i]]$covariance + 
													super.list$parameter.list$cluster.list[[i]]$cluster.mean})))
	y.range <- c(y.min,y.max)
	n.col <- ifelse(time,2,1)
	par(mfrow=c(1,n.col))
	plot(data.list$geo.dist[ind.mat],data.list$sample.covariance[ind.mat],
			ylab="sample covariance",xlab="geographic distance",
			ylim=y.range + diff(range(y.range))/10 * c(-1,1))
	abline(h=super.list$parameter.list$shared.mean,col="gray",lty=2)
	for(i in 1:super.list$model.options$n.clusters){
		points(data.list$geo.dist[ind.mat],
				super.list$parameter.list$cluster.list[[i]]$covariance[ind.mat] + 
				super.list$parameter.list$cluster.list[[i]]$cluster.mean,
				col=cluster.colors[i],pch=20,cex=0.7)
	}
	if(time){
		plot(data.list$time.dist[ind.mat],data.list$sample.covariance[ind.mat],
			ylab="sample covariance",xlab="temporal distance",
			ylim=y.range + diff(range(y.range))/10 * c(-1,1))
		abline(h=super.list$parameter.list$shared.mean,col="gray",lty=2)
		for(i in 1:super.list$model.options$n.clusters){
			points(data.list$time.dist[ind.mat],
					super.list$parameter.list$cluster.list[[i]]$covariance[ind.mat] + 
					super.list$parameter.list$cluster.list[[i]]$cluster.mean,
					col=cluster.colors[i],pch=20,cex=0.7)
		}
	}
	return(invisible("plotted"))
}

model.fit.plot <- function(data.list,super.list,time){
	n.col <- ifelse(time,3,2)
	par(mfrow=c(1,n.col),mar=c(4,4,2,2))
		plot(data.list$geo.dist,data.list$sample.covariance,
			xlab="geographic distance",
			ylab="covariance")
		points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
		legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
		if(time){
			plot(data.list$time.dist,data.list$sample.covariance,
				xlab="temporal distance",
				ylab="covariance")
			points(data.list$time.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red")
		}
		plot(data.list$sample.covariance,super.list$parameter.list$admixed.covariance,
				xlab="sample covariance",
				ylab ="parametric admixed covariance")
			abline(0,1,col="red")
	return(invisible(0))
}

plot.lnl <- function(super.list){
	mcmc.length <- length(which(!is.na(super.list$output.list$likelihood)))
	burnin <- 0.3*mcmc.length
	plot.vec <- burnin:mcmc.length
	plot(super.list$output.list$likelihood[plot.vec],xlab="",ylab="",type='l',
		ylim=c(min(super.list$output.list$likelihood[plot.vec],
					super.list$output.list$posterior.prob[plot.vec]),
				max(super.list$output.list$likelihood[plot.vec],
					super.list$output.list$posterior.prob[plot.vec])))
	points(super.list$output.list$posterior.prob[plot.vec],type='l',col="red")
	legend(x="bottomright",lty=1,col=c(1,2),legend=c("log likelihood","log posterior probability"),cex=0.6)
	return(invisible(0))
}

get.ylim <- function(super.list,output){
	mcmc.length <- length(which(!is.na(super.list$output.list$likelihood)))
	burnin <- 0.3*mcmc.length
	ylim <- range(output[,burnin:mcmc.length])
	return(ylim)
}

plot.cluster.cov.params <- function(super.list,all.colors,cluster.names,time){
	n.col <- ifelse(time,5,4)
	n.clusters <- length(super.list$parameter.list$cluster.list)
	par(mfrow=c(1,n.col),mar=c(1,2,2,1))
		matplot(t(super.list$output.list$cluster.params$cov.par1),type='l',lty=1,
				ylab="",col=all.colors[1:n.clusters],
				ylim=c(get.ylim(super.list,super.list$output.list$cluster.params$cov.par1)))
			mtext("cov.par1",side=3,padj=-0.1)
		matplot(t(super.list$output.list$cluster.params$cov.par2),type='l',lty=1,
				ylab="",col=all.colors[1:n.clusters],
				ylim=c(get.ylim(super.list,super.list$output.list$cluster.params$cov.par2)))
			mtext("cov.par2",side=3,padj=-0.1)
		matplot(t(super.list$output.list$cluster.params$cov.par3),type='l',lty=1,
				ylab="",col=all.colors[1:n.clusters],
				ylim=c(get.ylim(super.list,super.list$output.list$cluster.params$cov.par3)))
			mtext("cov.par3",side=3,padj=-0.1)
	if(time){
		matplot(t(super.list$output.list$cluster.params$cov.par4),type='l',lty=1,
				ylab="",col=all.colors[1:n.clusters],
				ylim=c(get.ylim(super.list,super.list$output.list$cluster.params$cov.par4)))
			mtext("cov.par4",side=3,padj=-0.1)
	}
		matplot(t(super.list$output.list$cluster.params$cluster.mean),type='l',lty=1,
				ylab="",col=all.colors[1:n.clusters],
				ylim=c(get.ylim(super.list,super.list$output.list$cluster.params$cluster.mean)))
			mtext("cluster.mean",side=3,padj=-0.1)
		legend(x="topright",col=all.colors[1:n.clusters],lty=1,legend=cluster.names)
	return(invisible(0))
}

make.admix.prop.mat <- function(super.list,k){
	x <- length(which(!is.na(super.list$output.list$posterior.prob)))
	admix.prop.mat <- matrix(NA,nrow=nrow(super.list$output.list$admix.proportions[[1]]),ncol=x)
	for(i in 1:x){
		admix.prop.mat[,i] <- super.list$output.list$admix.proportions[[i]][,k]
	}
	return(admix.prop.mat)
}

plot.admix.props <- function(super.list){
	n.clusters <- length(super.list$parameter.list$cluster.list)
	admix.prop.mat.list <- lapply(1:n.clusters,function(i){make.admix.prop.mat(super.list,i)})
	par(mfrow=c(1,n.clusters),mar=c(2,2,2,2))
		for(i in 1:length(admix.prop.mat.list)){
			matplot(t(admix.prop.mat.list[[i]]),type='l',ylim=c(0,1),
					main=paste("Cluster ",i,sep=""),ylab="admixture proportion")
		}
	return(invisible(0))
}

plot.nuggets <- function(super.list){
	par(mar=c(2,2,2,2))
	matplot(t(super.list$output.list$nuggets),type='l',
				ylim=get.ylim(super.list,super.list$output.list$nuggets),
				main="sample nuggets",
				ylab="nugget value")
	return(invisible(0))
}

plot.shared.mean.effects <- function(super.list,all.colors,cluster.names){
	n.clusters <- length(cluster.names)
	par(mar=c(2,2,2,2))
	matplot(t(super.list$output.list$cluster.params$cluster.mean),type='l',lty=1,
				ylab="",col=all.colors[1:n.clusters],
				ylim=c(get.ylim(super.list,super.list$output.list$cluster.params$cluster.mean)))
	points(super.list$output.list$shared.mean,type='l',col=1,lwd=2)
		legend(x="topright",col=c(1,all.colors[1:n.clusters]),
					lty=1,lwd=c(2,rep(1,n.clusters)),
					legend=c("shared mean",cluster.names),
					cex=0.6)
	return(invisible(0))
}

plot.acceptance.rates <- function(super.list,all.colors){
	n.clusters <- length(super.list$parameter.list$cluster.list)
	layout(mat=matrix(c(1,2,3,4,5,6,rep(7,6),rep(8,6)),nrow=3,ncol=6,byrow=TRUE))
	par(mar=c(0.5,3,2,0.5))
		plot(super.list$output.list$acceptance.rates$shared.mean,type='l',xlab="",ylab="",ylim=c(0,1),main="shared mean")
	par(mar=c(0.5,0.5,2,0.5))
		matplot(t(super.list$output.list$acceptance.rates$cov.par1),
				type='l',xlab="",ylab="",col=all.colors[1:n.clusters],
				ylim=c(0,1),main="cov.par 1",yaxt='n')
		matplot(t(super.list$output.list$acceptance.rates$cov.par2),
				type='l',xlab="",ylab="",col=all.colors[1:n.clusters],
				ylim=c(0,1),main="cov.par 2",yaxt='n')
		matplot(t(super.list$output.list$acceptance.rates$cov.par3),
				type='l',xlab="",ylab="",col=all.colors[1:n.clusters],
				ylim=c(0,1),main="cov.par 3",yaxt='n')
		matplot(t(super.list$output.list$acceptance.rates$cov.par4),
				type='l',xlab="",ylab="",col=all.colors[1:n.clusters],
				ylim=c(0,1),main="cov.par 4",yaxt='n')
		matplot(t(super.list$output.list$acceptance.rates$cluster.mean),
				type='l',xlab="",ylab="",col=all.colors[1:n.clusters],
				ylim=c(0,1),main="cluster.mean",yaxt='n')
	par(mar=c(0.5,3,3,0.5))
		matplot(t(super.list$output.list$acceptance.rates$admix.proportions),
				type='l',xlab="",ylab="",ylim=c(0,1))
			mtext(side=3,padj=2,text="admix proportions")
		matplot(t(super.list$output.list$acceptance.rates$nuggets),
				type='l',xlab="",ylab="",ylim=c(0,1))
			mtext(side=3,padj=2,text="nuggets")
	return(invisible(0))
}

make.all.the.plots <- function(dir,output.dir,K){
	setwd(dir)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
	time <- ifelse(length(unique(c(data.list$time.dist))) > 1, TRUE, FALSE)
	all.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
	#sample covariance fit
	pdf(file=paste(output.dir,"/","model.fit.",K,".pdf",sep=""),width=(8+time*4),height=5,pointsize=18)
		model.fit.plot(data.list,super.list,time)
	dev.off()
	pdf(file=paste(output.dir,"/","cluster.cov.curves.",K,".pdf",sep=""),width=(5+time*5),height=5,pointsize=18)
		plot.cluster.covariances(data.list,super.list,time)
	dev.off()
	pdf(file=paste(output.dir,"/","lnl.and.prob.",K,".pdf",sep=""),width=6,height=5,pointsize=18)
		plot.lnl(super.list)
	dev.off()
	pdf(file=paste(output.dir,"/","cluster.cov.params.",K,".pdf",sep=""),width=12,height=4,pointsize=18)
		plot.cluster.cov.params(super.list,all.colors,cluster.names,time)
	dev.off()
	pdf(file=paste(output.dir,"/","admix.props.",K,".pdf",sep=""),width=(4*super.list$model.options$n.clusters),height=4,pointsize=18)
		plot.admix.props(super.list)
	dev.off()
	pdf(file=paste(output.dir,"/","nuggets.",K,".pdf",sep=""),width=6,height=4,pointsize=18)
		plot.nuggets(super.list)
	dev.off()
	pdf(file=paste(output.dir,"/","shared.mean.effects.",K,".pdf",sep=""),width=6,height=4,pointsize=18)
		plot.shared.mean.effects(super.list,all.colors,cluster.names)
	dev.off()
	pdf(file=paste(output.dir,"/","acceptance.rates.",K,".pdf",sep=""),width=14,height=8.5,pointsize=18)
		plot.acceptance.rates(super.list,all.colors)
	dev.off()
	pdf(file=paste(output.dir,"/","pie.chart.map.",K,".pdf",sep=""),width=6,height=6,pointsize=18)	
		make.admix.pie.plot(super.list,data.list,all.colors,cluster.names,add=FALSE)
	dev.off()
	pdf(file=paste(output.dir,"/","structure.plot.",K,".pdf",sep=""),width=10,height=5,pointsize=18)
		make.structure.plot(data.list,super.list,sort.by=NULL,cluster.colors=all.colors)
	dev.off()
}

if(FALSE){
#spatial
for(i in 1:10){
	dir <- paste("~/Desktop/hgdp_analyses/spatial/k_",i,sep="")	#paste("~/Desktop/globe/globe_analyses/spatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
	pdf(file=paste("sp_admix_map_k=",i,".pdf",sep=""),width=10,height=5)
		map()
		make.admix.pie.plot(super.list,data.list,all.colors=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray"),
							cluster.names=paste(1:i),radii=2,add=TRUE,title=NULL,xlim=NULL,ylim=NULL)
	dev.off()
}
#nonspatial
for(i in 1:10){
	dir <- paste("~/Desktop/hgdp_analyses/nonspatial/k_",i,sep="")	#paste("~/Desktop/globe/globe_analyses/nonspatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
	pdf(file=paste("nsp_admix_map_k=",i,".pdf",sep=""),width=10,height=5)
		map()
		make.admix.pie.plot(super.list,data.list,all.colors=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray"),
							cluster.names=paste(1:i),radii=2,add=TRUE,title=NULL,xlim=NULL,ylim=NULL)
	dev.off()
}
# make.all.the.plots("~/Desktop/globe/spatial/k_1","~/Desktop/globe/spatial/k_1",1)
}




