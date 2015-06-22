source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/SpaceTimeStructureMix.R")
source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/viz_funcs.R")

make.nice.mfrow <- function(n.k){
	dims <- list("n.row" = ceiling(n.k/4),
 				 "n.col" = min(n.k,4))
	return(dims)	
}

rank.cor.matches <- function(pairwise,cors){
	tmp.cors <- cors
	tmp.pairwise <- pairwise
	old.K <- max(pairwise[,1])
	matches <- numeric(old.K)
	for(i in 1:old.K){
		old.K.match <- tmp.pairwise[which.max(tmp.cors),1]
		new.K.match <- tmp.pairwise[which.max(tmp.cors),2]
		matches[old.K.match] <- new.K.match
		drop.index <- which(tmp.pairwise[,1]==old.K.match | tmp.pairwise[,2]==new.K.match)
		tmp.cors <- tmp.cors[-drop.index]
		tmp.pairwise <- tmp.pairwise[-drop.index,]
	}
	return(matches)
}

get.cluster.matches <- function(ref.ad.props,super.list){
	new.ad.props <- super.list$parameter.list$admix.proportions
	old.K <- ncol(ref.ad.props)
	new.K <- ncol(new.ad.props)
	pairwise <- expand.grid(1:old.K,
							1:new.K,stringsAsFactors=FALSE)
	cors <- numeric(nrow(pairwise))
	for(i in 1:length(cors)){
		cors[i] <- cor(ref.ad.props[,pairwise[i,1]],
						new.ad.props[,pairwise[i,2]])
	}
	matches <- rank.cor.matches(pairwise,cors)
	return(matches)
}

get.color.order <- function(ref.ad.props,super.list,use.colors,all.colors){
	# recover()
	matches <- get.cluster.matches(ref.ad.props,super.list)
	new.K <- length(matches) + 1
	color.order <- c(matches,c(1:new.K)[!1:new.K %in% matches])
	use.colors <- use.colors[order(color.order)]
	use.colors[which(is.na(use.colors))] <- all.colors[new.K]
	return(use.colors)
}

plot.pie.grid <- function(dir,all.colors,n.k){
	# recover()
	setwd(dir)
	dims <- make.nice.mfrow(n.k)
	use.colors <- all.colors
	par(mfrow=c(dims$n.row,dims$n.col))
	for(i in 1:n.k){
		setwd(paste(dir,"/k_",i,sep=""))
		load(list.files(pattern="data.list"))
		load(list.files(pattern="output"))
		if(i > 2){
			use.colors <- get.color.order(ref.ad.props,super.list,use.colors,all.colors)
		}
		# cat(i,use.colors,"\n")
		ref.ad.props <- super.list$parameter.list$admix.proportions
		cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
		make.admix.pie.plot(super.list,data.list,use.colors,cluster.names,radii=2.7,add=FALSE,title=paste("K=",i,sep=""))
		setwd("..")
	}
}

plot.structure.plot.grid <- function(dir,n.k,all.colors,sample.order=NULL,sample.names=NULL,sort.by=NULL){
	#recover()
	setwd(dir)
	use.colors <- all.colors
	par(mfrow=c(n.k,1),mar=c(0.1,1.5,0.1,0.1))
	for(i in 1:n.k){
		setwd(paste(dir,"/k_",i,sep=""))
		load(list.files(pattern="data.list"))
		load(list.files(pattern="output"))
		if(i > 2){
			use.colors <- get.color.order(ref.ad.props,super.list,use.colors,all.colors)
		}
		# cat(i,use.colors,"\n")
		if(i > 2){
			cluster.order <- match(all.colors,use.colors)[which(!is.na(match(all.colors,use.colors)))]
		} else {
			cluster.order <- NULL
		}
		make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=cluster.order,sample.names=NULL,sort.by=NULL,cluster.colors=use.colors)
		ref.ad.props <- super.list$parameter.list$admix.proportions
		setwd("..")
	}
}

get.avg.params <- function(n.clusters,super.list,data.list,sampled.gens,no.IBD){
	avg.parameters <- list("admix.proportions" = NULL,
							"nuggets" = NULL,
							"cluster.list" = generate.clusters(n.clusters),
							"admixed.covariance" = NULL)
	avg.parameters$shared.mean <- mean(super.list$output.list$shared.mean[sampled.gens])
	avg.parameters$admix.proportions <-   Reduce("+",super.list$output.list$admix.proportions[sampled.gens])/length(sampled.gens)
	avg.parameters$nuggets <- rowMeans(super.list$output.list$nuggets[,sampled.gens])
		for(i in 1:n.clusters){
			if(!no.IBD){
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par1 <- mean(super.list$output.list$cluster.params$cov.par1[i,sampled.gens])
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par2 <- mean(super.list$output.list$cluster.params$cov.par2[i,sampled.gens])
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par3 <- mean(super.list$output.list$cluster.params$cov.par3[i,sampled.gens])
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par4 <- mean(super.list$output.list$cluster.params$cov.par4[i,sampled.gens])
				avg.parameters$cluster.list[[i]]$covariance <- cluster.covariance(data.list$geo.dist,data.list$time.dist,avg.parameters$cluster[[i]]$covariance.params)
			}else{
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par1 <- 0
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par2 <- 0
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par3 <- 0
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par4 <- 0
				avg.parameters$cluster.list[[i]]$covariance <- matrix(0,nrow=data.list$n.ind,ncol=data.list$n.ind)
			}
			avg.parameters$cluster.list[[i]]$cluster.mean <- mean(super.list$output.list$cluster.params$cluster.mean[i,sampled.gens])
			avg.parameters$cluster.list[[i]]$admix.prop.matrix<-  avg.parameters$admix.proportions[,i] %*% t(avg.parameters$admix.proportions[,i])
		}
	avg.parameters$admixed.covariance <- admixed.covariance(avg.parameters$cluster.list,n.clusters,avg.parameters$shared.mean,avg.parameters$nuggets)
	avg.parameters$inverse <- solve(avg.parameters$admixed.covariance)
	avg.parameters$determinant <- determinant(avg.parameters$admixed.covariance,logarithm=TRUE)$modulus
	return(avg.parameters)
}

n.params <- function(n.ind,k,time,spatial){
	n.params <- n.ind + 				#nuggets
				spatial*k*(3+time*1) +	#cov.pars, w/ or w/o time
				k + 					#cluster means
				n.ind * (k-1) + 		#admix proportions
				1						#shared mean
	return(n.params)
}

calc.DIC <- function(dir, burnin=500){
	# recover()
	setwd(dir)
	load(list.files(pattern="data.list"))
	load(list.files(pattern="output"))
	sampled.gens <- burnin:(max(which(!is.na(super.list$output.list$likelihood))))
	n.clusters <- super.list$model.options$n.clusters
	no.IBD <- super.list$model.options$no.st
	n.ind <- data.list$n.ind
	avg.parameters <- get.avg.params(n.clusters,super.list,data.list,sampled.gens,no.IBD)
	likelihood.at.avg <- calculate.likelihood.2(data.list,avg.parameters$inverse,avg.parameters$determinant)   ### D()
	avg.likelihood <- mean(super.list$output.list$likelihood[sampled.gens])  ###  Bar(D())
	pD <- (-2)*avg.likelihood - (-2)*likelihood.at.avg
	pV <- n.params(n.ind=data.list$n.ind,k=super.list$model.options$n.clusters,time=super.list$model.options$temporal.sampling,spatial=!super.list$model.options$no.st)			#var((-2)*super.list$output.list$likelihood[sampled.gens])/2
	DIC <- (-2*likelihood.at.avg) + 2*pV  # 2*pD
	cat(avg.likelihood, likelihood.at.avg, pD, pV,"\n")
	return(c(DIC, pD, pV))
}

plot.model.comp <- function(dir,n.runs=4){
	# recover()
	#Spatial vs. Nonspatial LnL comparison
	prob.vec.spatial <- rep(NA,n.runs)
	prob.vec.nonspatial <- rep(NA,n.runs)
	DIC.spatial <- numeric()
	DIC.nonspatial <- numeric()
	setwd(dir)
	for(k in 1:n.runs){
		sp.sub.dir <- paste(dir,"/spatial/k_",k,sep="")
		nsp.sub.dir <- paste(dir,"/nonspatial/k_",k,sep="")
		setwd(sp.sub.dir)
			load(list.files(pattern="output"))
			prob.vec.spatial[k] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
			DIC.spatial <- rbind(DIC.spatial,calc.DIC(sp.sub.dir))
		setwd(nsp.sub.dir)
			load(list.files(pattern="output"))
			prob.vec.nonspatial[k] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
			DIC.nonspatial <- rbind(DIC.nonspatial,calc.DIC(nsp.sub.dir))
	}
	layout(t(c(1:3)))
		plot(prob.vec.nonspatial,pch=19,col="blue",ylab="Posterior Probability",
				ylim=c(min(prob.vec.nonspatial,prob.vec.spatial,na.rm=TRUE),
						max(prob.vec.nonspatial, prob.vec.spatial,na.rm=TRUE)),
				xaxt='n',cex=2,xlab="")
			axis(side=1,at=1:n.runs,labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
			points(prob.vec.spatial,pch=19,col="green",cex=2)
			legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("nonspatial","spatial"))
			box(lwd=2)
		all.dev <- c(DIC.nonspatial[,1],DIC.spatial[,1])
		min.dev <- min(all.dev)
		DIC.nonspatial[,1] <- DIC.nonspatial[,1] - min.dev
		DIC.spatial[,1] <- DIC.spatial[,1] - min.dev
		plot(DIC.nonspatial[,1],pch=19,col="blue",ylab="DIC",ylim=c(0,max(all.dev - min.dev)),xlim=c(0.9,n.runs+0.5),xaxt='n',cex=2,xlab="") 
			points(DIC.spatial[,1],pch=19,col="green",cex=2)
			text(x=1:length(DIC.spatial[,1]) +.3,y=DIC.spatial[,1], format(DIC.spatial[,3],dig=3))
			text(x=1:length(DIC.nonspatial[,1])+.3,y=DIC.nonspatial[,1], format(DIC.nonspatial[,3],dig=3))
			axis(side=1,at=c(1:n.runs),labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
			box(lwd=2)
		plot(DIC.spatial[,1],pch=19,col="green",ylab="DIC",xaxt='n',cex=2,xlab="",xlim=c(0.9,n.runs+0.5)) 
			points(which.min(DIC.spatial[,1]),min(DIC.spatial[,1]),col="red",cex=2)
			text(x=1:length(DIC.spatial)+.3,y=DIC.spatial[,1], format(DIC.spatial[,3],dig=3))
			axis(side=1,at=c(1:n.runs),labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
			box(lwd=2)
	return(invisible(0))
}

plot.K.comp <- function(sub.dir,n.runs=4){
	# recover()
	#Cross-Run LnL comparison
	prob.vec <- rep(NA,n.runs)
	DIC.vec <- numeric()
	setwd(sub.dir)
	for(k in 1:n.runs){
		sub.sub.dir <- paste(sub.dir,"/k_",k,sep="")
		setwd(sub.sub.dir)
			load(list.files(pattern="output"))
			prob.vec[k] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
			DIC.vec <- rbind(DIC.vec,calc.DIC(sub.sub.dir))
	}
	par(mfrow=c(1,2),mar=c(2.5,2.5,1,1))
		plot(prob.vec,pch=19,col="blue",ylab="Posterior Probability",
				xaxt='n',cex=2,xlab="")
			axis(side=1,at=1:n.runs,labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
			box(lwd=2)
		DIC.vec[,1] <- DIC.vec[,1] - min(DIC.vec[,1])
		plot(DIC.vec[,1],pch=19,col="blue",ylab="DIC",xlim=c(0.9,n.runs+0.5),xaxt='n',cex=2,xlab="")
			points(which.min(DIC.vec[,1]),min(DIC.vec[,1]),col="red",cex=2)
			text(x=1:length(DIC.vec[,1]) +.3,y=DIC.vec[,1], format(DIC.vec[,3],dig=3))
			axis(side=1,at=c(1:n.runs),labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
			box(lwd=2)
	return(invisible(0))
}

plot.admix.pie.maps <- function(sub.dir,all.colors,n.runs,map.xlim=NULL,map.ylim=NULL,pdf.size=NULL){
	# recover()
	require(maps)
	setwd(sub.dir)
	use.colors <- all.colors
	for(i in 1:n.runs){
		sub.sub.dir <- paste(sub.dir,"/k_",i,sep="")
		setwd(sub.sub.dir)
			load(list.files(pattern="output"))
			load(list.files(pattern="data.list"))
		if(i > 2){
			use.colors <- get.color.order(ref.ad.props,super.list,use.colors,all.colors)
		}
		ref.ad.props <- super.list$parameter.list$admix.proportions
		cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
		if(is.null(map.xlim)){
			map.xlim <- range(data.list$geo.coords[,1]) + diff(range(data.list$geo.coords[,1]))/25 * c(-1,1)
		}
		if(is.null(map.ylim)){
			map.ylim <- range(data.list$geo.coords[,2]) + diff(range(data.list$geo.coords[,2]))/25 * c(-1,1)
		}
		if(is.null(pdf.size)){
			scale <- diff(map.xlim)%/%10
			pdf.size <- c(diff(map.xlim),diff(map.ylim)) / c(scale,scale)
		}
		setwd("..")
		pdf(file=paste("pie.admix.map.k_",i,".pdf",sep=""),width=pdf.size[1],height=pdf.size[2])
			par(mar=c(1,1,1,1))
			map(database="world",xlim=map.xlim,ylim=map.ylim,col="darkgray")
			make.admix.pie.plot(super.list,data.list,use.colors,cluster.names,radii=2.7,add=TRUE,title=paste("K=",i,sep=""))
			box(lwd=2)
		dev.off()
	}
	return(invisible(0))
}

make.all.metaplots <- function(dir,output.dir,K,sample.order=NULL,sample.names=NULL,sort.by=NULL){
	setwd(dir)
	sp.sub.dir <- paste(dir,"/spatial",sep="")
	nsp.sub.dir <- paste(dir,"/nonspatial",sep="")
	all.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	pdf(file=paste(output.dir,"/spatial.admix.pie.grid.pdf",sep=""),width=10,height=3*ceiling(K/4))
		plot.pie.grid(dir=sp.sub.dir,all.colors,K)
	dev.off()
	if(file.exists(nsp.sub.dir)){
	pdf(file=paste(output.dir,"/nonspatial.admix.pie.grid.pdf",sep=""),width=10,height=3*ceiling(K/4))
		plot.pie.grid(dir=nsp.sub.dir,all.colors,K)
	dev.off()
	}
	pdf(file=paste(output.dir,"/spatial.str.plots.pdf",sep=""),width=10,height=2*K)
		plot.structure.plot.grid(dir=sp.sub.dir,K,all.colors,sample.order,sample.names,sort.by)
	dev.off()
	if(file.exists(nsp.sub.dir)){
	pdf(file=paste(output.dir,"/nonspatial.str.plots.pdf",sep=""),width=10,height=2*K)
		plot.structure.plot.grid(dir=nsp.sub.dir,K,all.colors,sample.order,sample.names,sort.by)
	dev.off()
	}
	if(file.exists(nsp.sub.dir)){
	pdf(file=paste(output.dir,"/model.comp.pdf",sep=""),width=10,height=5)
		plot.model.comp(dir,n.runs=K)
	dev.off()
	}
	pdf(file=paste(output.dir,"/spatial.K.comp.pdf",sep=""),width=10,height=5)
		plot.K.comp(sub.dir=sp.sub.dir,n.runs=K)
	dev.off()
	if(file.exists(nsp.sub.dir)){
	pdf(file=paste(output.dir,"/nonspatial.K.comp.pdf",sep=""),width=10,height=5)
		plot.K.comp(sub.dir=nsp.sub.dir,n.runs=K)
	dev.off()
	}
	plot.admix.pie.maps(sub.dir=sp.sub.dir,all.colors=all.colors,n.runs=K)
	if(file.exists(nsp.sub.dir)){
		plot.admix.pie.maps(sub.dir=nsp.sub.dir,all.colors=all.colors,n.runs=K)
	}
	return(invisible(0))
}

# make.all.metaplots("~/Desktop/new_runs","~/desktop",10)

make.all.metaplots("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/globe/globe_analyses","~/desktop",10)

make.all.metaplots("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/hgdp/hgdp_analyses","~/desktop",10)

make.all.metaplots("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers/analyses/new_runs","~/desktop",10)
