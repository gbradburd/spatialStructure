################################################################
################################################################
#	figs for spatialStructure presentation
################################################################
################################################################

calc.DIC<-function(data.list.file,super.list.file, burnin=500){
	
load(data.list.file)
load(super.list.file)
runs<-burnin:(max(which(!is.na(super.list$output.list$likelihood))))
	n.clusters<-super.list$model.options$n.clusters
	no.IBD<-super.list$model.options$no.st
	n.ind <- nrow(data.list$geo.dist)
	avg.parameters <- list("admix.proportions" = NULL,
						"nuggets" = NULL,
						"cluster.list" = generate.clusters(n.clusters),
						"admixed.covariance" = NULL)
						
		avg.parameters$shared.mean <-  mean(super.list$output.list$shared.mean[runs])   #min(data$sample.covariance)
		avg.parameters$admix.proportions <-   Reduce("+",super.list$output.list$admix.proportions[runs])/length(runs)   #gtools::rdirichlet(n = n.ind,alpha = rep(1,model.options$n.clusters))
		avg.parameters$nuggets <- rowMeans(super.list$output.list$nuggets[,runs])     #rexp(n.ind)
		for(i in 1:n.clusters){
			if(!no.IBD){
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par1 <- mean(super.list$output.list$cluster.params$cov.par1[i,runs])
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par2 <- mean(super.list$output.list$cluster.params$cov.par2[i,runs])
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par3 <- mean(super.list$output.list$cluster.params$cov.par3[i,runs])
				avg.parameters$cluster.list[[i]]$covariance <- cluster.covariance(data.list$geo.dist,data.list$time.dist,avg.parameters$cluster[[i]]$covariance.params)
			}else{
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par1 <- 0
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par2 <- 0
				avg.parameters$cluster.list[[i]]$covariance.params$cov.par3 <- 0
				avg.parameters$cluster.list[[i]]$covariance <- matrix(0,nrow=n.ind,ncol=n.ind)
			}
			avg.parameters$cluster.list[[i]]$cluster.mean <- mean(super.list$output.list$cluster.params$cluster.mean[i,runs])
			avg.parameters$cluster.list[[i]]$admix.prop.matrix<-  avg.parameters$admix.proportions[,i] %*% t(avg.parameters$admix.proportions[,i])

		}
		
	avg.parameters$admixed.covariance <- admixed.covariance(avg.parameters$cluster.list,n.clusters,avg.parameters$shared.mean,avg.parameters$nuggets)
	avg.parameters$inverse <- solve(avg.parameters$admixed.covariance)
	avg.parameters$determinant <- determinant(avg.parameters$admixed.covariance,logarithm=TRUE)$modulus
	likelihood.at.avg <- calculate.likelihood.2(data.list,avg.parameters$inverse,avg.parameters$determinant)   ### D()
	avg.likelihood <- mean(super.list$output.list$likelihood[runs])  ###  Bar(D())

	pD<-  (-2)*avg.likelihood - (-2)*likelihood.at.avg
	pV<-var((-2)*super.list$output.list$likelihood[runs])/2
	DIC <-  (-2*likelihood.at.avg) +  2*pV  # 2*pD
	cat(avg.likelihood, likelihood.at.avg, pD, pV,"\n")
return(c(DIC, pD,pV))
}


plot.pie.grid<-function(pdf_file,my.path,geo.coords,data.list,super.list){
	pdf(file=pdf_file,width=15,height=6,pointsize=12)
	all.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen")
	#pdf(file="~/desktop/lattice_pops_nonspatial_adprop_piemaps.pdf",width=8,height=8,pointsize=12)
	x.lim <- c(min(geo.coords[,1]) - 1, max(geo.coords[,1]) + 1)
	y.lim <- c(min(geo.coords[,2]) - 1, max(geo.coords[,2]) + 1)
	radius <- 3.5
	layout(t(1:3));  par(mar=c(1,1,1,1)) # par(mfrow=c(2,2),mar=c(1,1,1,1))
	
	require(caroline)
	#K_1
	# load("spatial/k_1/k_1_output.Robj")
	# sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
	# cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
	# color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
	# pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
	# names(pie.list) <- sample.names
	# pies(pie.list,x0=geo.coords[,1],
					# y0=geo.coords[,2],
					# color.table=color.tab,border="black",radii=radius,
					# xlab="",ylab="",main="K = 1",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim) 
	# box(lwd=2)
	#K_2
	load(paste(my.path,"/k_2/k_2_output.Robj",sep=""))
	sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
	cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
	color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
	pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
	names(pie.list) <- sample.names
	pies(pie.list,x0=geo.coords[,1],
					y0=geo.coords[,2],
					color.table=color.tab,border="black",radii=radius,
					xlab="",ylab="",main="K = 2",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim)
	box(lwd=2)
	#K_3
	load(paste(my.path,"/k_3/k_3_output.Robj",sep=""))
	sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
	cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
	color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
	pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
	names(pie.list) <- sample.names
	pies(pie.list,x0=geo.coords[,1],
					y0=geo.coords[,2],
					color.table=color.tab,border="black",radii=radius,
					xlab="",ylab="",main="K = 3",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim)
	box(lwd=2)
	#K_4
	load(paste(my.path,"/k_4/k_4_output.Robj",sep=""))
	sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
	cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
	color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
	pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
	names(pie.list) <- sample.names
	pies(pie.list,x0=geo.coords[,1],
					y0=geo.coords[,2],
					color.table=color.tab,border="black",radii=radius,
					xlab="",ylab="",main="K = 4",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim)
	box(lwd=2)
	dev.off()
}

plot.pie.map <- function(my.path,dir.name,geo.coords,fig.filename,map.x.lim=NULL,map.y.lim=NULL,pdf.width=6,pdf.height=5){
	# recover()
	require(caroline)
	require(maps)
	setwd(dir.name)
	load(list.files(pattern="data.list"))
	load(list.files(pattern="output"))
	setwd(my.path)
	all.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen")
	sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
	cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
	color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
	pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
	names(pie.list) <- sample.names
		if(is.null(map.x.lim)){
			map.x.lim <- range(geo.coords[,1]) + c(-0.1,0.1) * diff(range(geo.coords[,1]))
		}
		if(is.null(map.x.lim)){
			map.y.lim <- range(geo.coords[,2]) + c(-0.1,0.1) * diff(range(geo.coords[,2]))
		}
	# op <- par()$new
	pdf(file=fig.filename,width=pdf.width,height=pdf.height)
		map(database="world",xlim=map.x.lim,ylim=map.y.lim)
		par(new=TRUE)
		pies(pie.list,x0=geo.coords[,1],
					y0=geo.coords[,2],
					color.table=color.tab,border="black",radii=3,
					xlab="",ylab="",lty=1,density=NULL)
		title(main=paste("K = ",super.list$model.options$n.clusters,sep=""))
		box(lwd=2)
	dev.off()
		# par(op)
	return(invisible("done"))
}

plot.model.comp<-function(my.path,pdf.names,n.runs=4){
	# recover()
	#Spatial vs. Nonspatial LnL comparison
	prob.vec.spatial <- rep(NA,n.runs)
	prob.vec.nonspatial <- rep(NA,n.runs)
	DIC.spatial<-numeric()
	DIC.nonspatial<-numeric()
	for(k in 1:n.runs){
		output.output<-paste(my.path, "spatial/k_", k, "/k_", k, "_output.Robj", sep="")
	
		load(output.output)
	#	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/line_sim")
	#	setwd(paste("spatial/k_",i,sep=""))
	#		load(list.files(pattern="output"))
			prob.vec.spatial[k] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	
		output.data<-paste(my.path, "spatial/k_", k,"/data.list.Robj",sep="")
		DIC.spatial<-rbind(DIC.spatial,calc.DIC(output.data,output.output))
	
	#	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/line_sim")
	#	setwd(paste("nonspatial/k_",i,sep=""))
	#		load(list.files(pattern="output"))

			output.output<-paste(my.path, "nonspatial/k_", k, "/k_", k, "_output.Robj", sep="")
		tryCatch({
			load(output.output) ;
			prob.vec.nonspatial[k] <- max(super.list$output.list$posterior.prob,na.rm=TRUE) ;
			output.data<-paste(my.path, "nonspatial/k_", k,"/data.list.Robj",sep="") ;
			DIC.nonspatial<-rbind(DIC.nonspatial,calc.DIC(output.data,output.output))},error=function(e){
				cat("folder K=",k,"does not exist\n")
		})
	}

	pdf(file=pdf.names,width=8,height=3,pointsize=11)
	layout(t(1:3))
	plot(prob.vec.nonspatial,pch=19,col="blue",ylab="Posterior Probability",
			ylim=c(min(prob.vec.nonspatial,prob.vec.spatial,na.rm=TRUE),
					max(prob.vec.nonspatial, prob.vec.spatial,na.rm=TRUE)),
			xaxt='n',cex=2,xlab="")
		axis(side=1,at=1:n.runs,labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
		points(prob.vec.spatial,pch=19,col="green",cex=2)
		box(lwd=2)
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("nonspatial","spatial"))		
		all.dev<-c(DIC.nonspatial[,1],DIC.spatial[,1])
		min.dev<-min(all.dev)
		DIC.nonspatial[,1]<-DIC.nonspatial[,1] - min.dev
		DIC.spatial[,1]<-DIC.spatial[,1] - min.dev
		plot(DIC.nonspatial[,1],pch=19,col="blue",ylab="DIC",ylim=c(0,max(all.dev - min.dev)),xlim=c(0.9,n.runs+0.5),xaxt='n',cex=2,xlab="") 
		points(DIC.spatial[,1],pch=19,col="green",cex=2)
		text(x=1:length(DIC.spatial[,1]) +.3,y=DIC.spatial[,1], format(DIC.spatial[,3],dig=3))
		text(x=1:length(DIC.nonspatial[,1])+.3,y=DIC.nonspatial[,1], format(DIC.nonspatial[,3],dig=3))
		axis(side=1,at=c(1:n.runs),labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
		plot(DIC.spatial[,1],pch=19,col="green",ylab="DIC",xaxt='n',cex=2,xlab="",xlim=c(0.9,n.runs+0.5)) 
		text(x=1:length(DIC.spatial)+.3,y=DIC.spatial[,1], format(DIC.spatial[,3],dig=3))
		axis(side=1,at=c(1:n.runs),labels=unlist(lapply(seq_along(1:n.runs),function(i){paste("K=",i,sep="")})))
	#plot(diff(prob.vec.nonspatial),pch=19,col="blue",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	#	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	#	box(lwd=2)
	#plot(diff(prob.vec.spatial),pch=19,col="green",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	#	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	#	box(lwd=2)
	dev.off()
}

if(file.exists("~/Dropbox/Students/gideon/")){
	my.path <- "~/Dropbox/Students/gideon/"
} else {
	my.path <- "~/Desktop/Dropbox/InspectorSpaceTime/"
}

load(paste(my.path,"spatialStructure/datasets/sims/glacier_sim/spatial/k_1/ms_glac_dataset.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/glacier_sim/spatial/k_1/k_1_output.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/glacier_sim/spatial/k_1/data.list.Robj",sep=""))
geo.coords <- ms_glac_dataset$geo.coords
plot.pie.grid(pdf_file="~/Downloads/glacier_spatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/glacier_sim/spatial/",sep=""),
				geo.coords,data.list,super.list)
plot.pie.grid(pdf_file="~/Downloads/glacier_nonspatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/glacier_sim/nonspatial/",sep=""),
				geo.coords,data.list,super.list)

load(paste(my.path,"spatialStructure/datasets/sims/barrier_sim/spatial/k_1/ms.barrier.dataset.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/barrier_sim/spatial/k_1/k_1_output.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/barrier_sim/spatial/k_1/data.list.Robj",sep=""))
geo.coords <- spacemix.dataset$population.coordinates
plot.pie.grid(pdf_file="~/Downloads/barrier_spatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/barrier_sim/spatial/",sep=""),
				geo.coords,data.list,super.list)
plot.pie.grid(pdf_file="~/Downloads/barrier_nonspatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/barrier_sim/nonspatial/",sep=""),
				geo.coords,data.list,super.list)

load(paste(my.path,"spatialStructure/datasets/sims/lattice_sim/spatial/k_1/latticepops.ms.dataset.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/lattice_sim/spatial/k_1/k_1_output.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/lattice_sim/spatial/k_1/data.list.Robj",sep=""))
geo.coords <- spacemix.dataset$population.coordinates
plot.pie.grid(pdf_file="~/Downloads/lattice_spatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/lattice_sim/spatial/",sep=""),
				geo.coords,data.list,super.list)
plot.pie.grid(pdf_file="~/Downloads/lattice_nonspatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/lattice_sim/nonspatial/",sep=""),
				geo.coords,data.list,super.list)

load(paste(my.path,"spatialStructure/datasets/sims/line_sim/spatial/k_1/linepops.ms.dataset.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/line_sim/spatial/k_1/k_1_output.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/sims/line_sim/spatial/k_1/data.list.Robj",sep=""))
geo.coords <- spacemix.dataset$population.coordinates
plot.pie.grid(pdf_file="~/Downloads/line_spatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/line_sim/spatial/",sep=""),
				geo.coords,data.list,super.list)
plot.pie.grid(pdf_file="~/Downloads/line_nonspatial.pdf",
				my.path=paste(my.path,"spatialStructure/datasets/sims/line_sim/nonspatial/",sep=""),
				geo.coords,data.list,super.list)
paste(my.path,"spatialStructure/datasets/sims/line_sim/spatial/k_1/linepops.ms.dataset.Robj",sep="")

source(paste(my.path,"spatialStructure/code/SpaceTimeStructureMix.R",sep=""))
plot.model.comp(my.path=paste(my.path,"spatialStructure/datasets/sims/lattice_sim/",sep=""),pdf.names="~/Downloads/Lattice_comp.pdf")
plot.model.comp(my.path=paste(my.path,"spatialStructure/datasets/sims/line_sim/",sep=""),pdf.names="~/Downloads/Line_comp.pdf")
plot.model.comp(my.path=paste(my.path,"spatialStructure/datasets/sims/barrier_sim/",sep=""),pdf.names="~/Downloads/Barrier_comp.pdf")
plot.model.comp(my.path=paste(my.path,"spatialStructure/datasets/sims/glacier_sim/",sep=""),pdf.names="~/Downloads/Glacier_comp.pdf")
plot.model.comp(my.path=paste(my.path,"spatialStructure/datasets/warblers/analyses/",sep=""),pdf.names="~/Downloads/warblers_comp.pdf")


#Warblers
load(paste(my.path,"spatialStructure/datasets/warblers/warbler_ind_dataset.Robj",sep=""))
	geo.coords <- warbler.ind.coords
	#Spatial
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/spatial/k_2",
					geo.coords,fig.filename="spatialStructure/figs/warb_k2_s.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/spatial/k_3",
					geo.coords,fig.filename="spatialStructure/figs/warb_k3_s.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/spatial/k_4",
					geo.coords,fig.filename="spatialStructure/figs/warb_k4_s.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/spatial/k_5",
					geo.coords,fig.filename="spatialStructure/figs/warb_k5_s.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/spatial/k_6",
					geo.coords,fig.filename="spatialStructure/figs/warb_k6_s.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/spatial/k_7",
					geo.coords,fig.filename="spatialStructure/figs/warb_k7_s.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/spatial/k_8",
					geo.coords,fig.filename="spatialStructure/figs/warb_k8_s.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)

	#Nonspatial
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/nonspatial/k_2",
					geo.coords,fig.filename="spatialStructure/figs/warb_k2_ns.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)	
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/nonspatial/k_3",
					geo.coords,fig.filename="spatialStructure/figs/warb_k3_ns.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)	
	plot.pie.map(my.path=my.path,dir.name="spatialStructure/datasets/warblers/analyses/nonspatial/k_4",
					geo.coords,fig.filename="spatialStructure/figs/warb_k4_ns.pdf",
					map.x.lim=c(37,110),map.y.lim=c(25,59),pdf.width=6,pdf.height=5)	

plot.model.comp(my.path=paste(my.path,"spatialStructure/datasets/warblers/analyses/",sep=""),pdf.names="~/Downloads/warblers_comp.pdf",n.runs=8)

load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_2/data.list.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_2/k_2_output.Robj",sep=""))
	inds <- row.names(warbler.ind.coords)
		inds[11] <- "Vir-STvi1"
		inds[12] <- "Vir-STvi2"
		inds[13] <- "Vir-STvi3"
	ind.subspp <- unlist(strsplit(inds,"-"))[seq(1,190,2)]
	inds.col <- numeric(length(ind.subspp))
		inds.col[grepl("Vir",ind.subspp)] <- "dodgerblue2"
		inds.col[grepl("Ni",ind.subspp)] <- "slateblue4"
		inds.col[grepl("Lud",ind.subspp)] <- "mediumseagreen"
		inds.col[grepl("Tro",ind.subspp)] <- "gold"
		inds.col[grepl("Obs",ind.subspp)] <- "orange"
		inds.col[grepl("Plu",ind.subspp)] <- "red"

		plot.inds <- gsub(" ","",inds)
		plot.inds <- gsub("[[:digit:]]","",plot.inds)
		plot.inds <- gsub(c("Plu-"),"",plot.inds)
		plot.inds <- gsub(c("Vir-"),"",plot.inds)
		plot.inds <- gsub(c("Ni-"),"",plot.inds)
		plot.inds <- gsub(c("Lud-"),"",plot.inds)
		plot.inds <- gsub(c("Tro-"),"",plot.inds)
		plot.inds <- gsub(c("Obs-"),"",plot.inds)
		plot.inds <- gsub(c("vi"),"",plot.inds)


make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=NULL,
						sample.names=NULL,sort.by=NULL)
for(i in 1:length(inds)){
	axis(side=1,at=seq(0.5,length(inds)-0.5,1)[i],labels=plot.inds[i],col.axis=inds.col[i],las=2,cex.axis=0.6,font=2)
}

load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_3/k_3_output.Robj",sep=""))						
make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=NULL,
						sample.names=NULL,sort.by=NULL)
for(i in 1:length(inds)){
	axis(side=1,at=seq(0.5,length(inds)-0.5,1)[i],labels=plot.inds[i],col.axis=inds.col[i],las=2,cex.axis=0.6,font=2)
}

load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_4/k_4_output.Robj",sep=""))						
make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=NULL,
						sample.names=NULL,sort.by=NULL)
for(i in 1:length(inds)){
	axis(side=1,at=seq(0.5,length(inds)-0.5,1)[i],labels=plot.inds[i],col.axis=inds.col[i],las=2,cex.axis=0.6,font=2)
}

load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_5/k_5_output.Robj",sep=""))						
make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=NULL,
						sample.names=NULL,sort.by=NULL)
for(i in 1:length(inds)){
	axis(side=1,at=seq(0.5,length(inds)-0.5,1)[i],labels=plot.inds[i],col.axis=inds.col[i],las=2,cex.axis=0.6,font=2)
}


load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_6/k_6_output.Robj",sep=""))						
make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=NULL,
						sample.names=NULL,sort.by=NULL)
for(i in 1:length(inds)){
	axis(side=1,at=seq(0.5,length(inds)-0.5,1)[i],labels=plot.inds[i],col.axis=inds.col[i],las=2,cex.axis=0.6,font=2)
}


load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_7/k_7_output.Robj",sep=""))						
make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=NULL,
						sample.names=NULL,sort.by=NULL)
for(i in 1:length(inds)){
	axis(side=1,at=seq(0.5,length(inds)-0.5,1)[i],labels=plot.inds[i],col.axis=inds.col[i],las=2,cex.axis=0.6,font=2)
}


load(paste(my.path,"spatialStructure/datasets/warblers/analyses/spatial/k_8/k_8_output.Robj",sep=""))						
make.structure.plot(data.list,super.list,sample.order=NULL,cluster.order=NULL,
						sample.names=NULL,sort.by=NULL)
for(i in 1:length(inds)){
	axis(side=1,at=seq(0.5,length(inds)-0.5,1)[i],labels=plot.inds[i],col.axis=inds.col[i],las=2,cex.axis=0.6,font=2)
}

#HAAK
metadata <- read.table(paste(my.path,"spatialStructure/datasets/HumanData/analyses/human_sample_metadata.txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
source(paste(my.path,"spatialStructure/code/viz_funcs.R",sep=""))

#spatiotemporal
load(paste(my.path,"spatialStructure/datasets/HumanData/analyses/temporal/spatial/k_2/haak_k2_st_output.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/HumanData/analyses/temporal/spatial/k_2/data.list.Robj",sep=""))
make.structure.plot(data.list,super.list,
					sample.order=rev(order(metadata$time)),
					cluster.order=NULL,
					sample.names=row.names(data.list$sample.covariance),
					sort.by=NULL)


load(paste(my.path,"spatialStructure/datasets/HumanData/analyses/temporal/spatial/k_3/haak_k3_st_output.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/HumanData/analyses/temporal/spatial/k_3/data.list.Robj",sep=""))
make.structure.plot(data.list,super.list,
					sample.order=rev(order(metadata$time)),
					cluster.order=NULL,
					sample.names=row.names(data.list$sample.covariance),
					sort.by=NULL)


load(paste(my.path,"spatialStructure/datasets/HumanData/analyses/temporal/spatial/k_4/haak_k4_st_output.Robj",sep=""))
load(paste(my.path,"spatialStructure/datasets/HumanData/analyses/temporal/spatial/k_4/data.list.Robj",sep=""))
make.structure.plot(data.list,super.list,
					sample.order=rev(order(metadata$time)),
					cluster.order=NULL,
					sample.names=row.names(data.list$sample.covariance),
					sort.by=NULL)







