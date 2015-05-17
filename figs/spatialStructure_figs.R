################################################################
################################################################
#	figs for spatialStructure presentation
################################################################
################################################################

################################
#Line scenario
################################
setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/line_sim/")

#NONSPATIAL
load("nonspatial/k_1/linepops.ms.dataset.Robj")
load("nonspatial/k_1/data.list.Robj")

geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")

pdf(file="~/desktop/line_pops_nonspatial_samplecov.pdf",width=10,height=3,pointsize=12)
par(mfrow=c(1,4),oma=c(0.5,0.5,0.5,0.5))
#K_1
load("nonspatial/k_1/k_1_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="covariance",main="K = 1")
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
	box(lwd=2)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
#K_2
load("nonspatial/k_2/k_2_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 2")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)

mtext(text="geographic distance",side=1,adj=4,padj=3)

#K_3
load("nonspatial/k_3/k_3_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
#K_4
load("nonspatial/k_4/k_4_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
dev.off()



pdf(file="~/desktop/line_pops_nonspatial_adprop_piemaps.pdf",width=4,height=6,pointsize=12)
par(mfrow=c(4,1),mar=c(1,1,1,1))
require(caroline)
#K_1
load("nonspatial/k_1/k_1_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 1",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
#K_2
load("nonspatial/k_2/k_2_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 2",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
#K_3
load("nonspatial/k_3/k_3_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 3",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
#K_4
load("nonspatial/k_4/k_4_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 4",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
dev.off()




#SPATIAL
load("spatial/k_1/linepops.ms.dataset.Robj")
load("spatial/k_1/data.list.Robj")

geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")

pdf(file="~/desktop/line_pops_spatial_samplecov.pdf",width=10,height=3,pointsize=12)
par(mfrow=c(1,4),oma=c(0.5,0.5,0.5,0.5))
#K_1
load("spatial/k_1/k_1_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="covariance",main="K = 1")
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
	box(lwd=2)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
#K_2
load("spatial/k_2/k_2_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 2")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)

mtext(text="geographic distance",side=1,adj=4,padj=3)

#K_3
load("spatial/k_3/k_3_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
#K_4
load("spatial/k_4/k_4_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
dev.off()



pdf(file="~/desktop/line_pops_spatial_adprop_piemaps.pdf",width=4,height=6,pointsize=12)
par(mfrow=c(4,1),mar=c(1,1,1,1))
require(caroline)
#K_1
load("spatial/k_1/k_1_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 1",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
#K_2
load("spatial/k_2/k_2_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 2",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
#K_3
load("spatial/k_3/k_3_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 3",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
#K_4
load("spatial/k_4/k_4_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3.5,
				xlab="",ylab="",main="K = 4",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=c(-0.5,20.5))
box(lwd=2)
dev.off()




#Spatial vs. Nonspatial LnL comparison
prob.vec.spatial <- numeric(4)
prob.vec.nonspatial <- numeric(4)
for(i in 1:4){
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/line_sim")
	setwd(paste("spatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.spatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/line_sim")
	setwd(paste("nonspatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.nonspatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/line_sim")
}

pdf(file="~/desktop/line_pops_spatial_vs_nonspatial.pdf",width=8,height=3,pointsize=11)
par(mfrow=c(1,3))
plot(prob.vec.nonspatial,pch=19,col="blue",ylab="Posterior Probability",
		ylim=c(min(prob.vec.nonspatial,prob.vec.spatial),
				max(prob.vec.nonspatial, prob.vec.spatial)),
		xaxt='n',cex=2,xlab="")
	axis(side=1,at=c(1,2,3,4),labels=c("K=1","K=2","K=3","K=4"))
	points(prob.vec.spatial,pch=19,col="green",cex=2)
	box(lwd=2)
	legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("nonspatial","spatial"))
plot(diff(prob.vec.nonspatial),pch=19,col="blue",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
plot(diff(prob.vec.spatial),pch=19,col="green",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
dev.off()

################################
#Lattice scenario
################################
setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/lattice_sim/")

#NONSPATIAL
load("nonspatial/k_1/latticepops.ms.dataset.Robj")
load("nonspatial/k_1/data.list.Robj")

geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")

pdf(file="~/desktop/lattice_pops_nonspatial_samplecov.pdf",width=10,height=3,pointsize=12)
par(mfrow=c(1,4),oma=c(0.5,0.5,0.5,0.5))
#K_1
load("nonspatial/k_1/k_1_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="covariance",main="K = 1")
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
	box(lwd=2)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
#K_2
load("nonspatial/k_2/k_2_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 2")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)

mtext(text="geographic distance",side=1,adj=4,padj=3)

#K_3
load("nonspatial/k_3/k_3_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
#K_4
load("nonspatial/k_4/k_4_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
dev.off()



pdf(file="~/desktop/lattice_pops_nonspatial_adprop_piemaps.pdf",width=8,height=8,pointsize=12)
x.lim <- c(min(geo.coords[,1]) - 1, max(geo.coords[,1]) + 1)
y.lim <- c(min(geo.coords[,2]) - 1, max(geo.coords[,2]) + 1)
radius <- 3.5
par(mfrow=c(2,2),mar=c(1,1,1,1))
require(caroline)
#K_1
load("nonspatial/k_1/k_1_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=radius,
				xlab="",ylab="",main="K = 1",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim) 
box(lwd=2)
#K_2
load("nonspatial/k_2/k_2_output.Robj")
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
load("nonspatial/k_3/k_3_output.Robj")
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
load("nonspatial/k_4/k_4_output.Robj")
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




#SPATIAL
load("spatial/k_1/latticepops.ms.dataset.Robj")
load("spatial/k_1/data.list.Robj")

geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")

pdf(file="~/desktop/lattice_pops_spatial_samplecov.pdf",width=10,height=3,pointsize=12)
par(mfrow=c(1,4),oma=c(0.5,0.5,0.5,0.5))
#K_1
load("spatial/k_1/k_1_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="covariance",main="K = 1")
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
	box(lwd=2)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
#K_2
load("spatial/k_2/k_2_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 2")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)

mtext(text="geographic distance",side=1,adj=4,padj=3)

#K_3
load("spatial/k_3/k_3_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
#K_4
load("spatial/k_4/k_4_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
dev.off()



pdf(file="~/desktop/lattice_pops_spatial_adprop_piemaps.pdf",width=8,height=8,pointsize=12)
x.lim <- c(min(geo.coords[,1]) - 1, max(geo.coords[,1]) + 1)
y.lim <- c(min(geo.coords[,2]) - 1, max(geo.coords[,2]) + 1)
radius <- 3.5
par(mfrow=c(2,2),mar=c(1,1,1,1))
require(caroline)
#K_1
load("spatial/k_1/k_1_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=radius,
				xlab="",ylab="",main="K = 1",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim) 
box(lwd=2)
#K_2
load("spatial/k_2/k_2_output.Robj")
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
load("spatial/k_3/k_3_output.Robj")
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
load("spatial/k_4/k_4_output.Robj")
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



#Spatial vs. Nonspatial LnL comparison
prob.vec.spatial <- numeric(4)
prob.vec.nonspatial <- numeric(4)
for(i in 1:4){
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/lattice_sim")
	setwd(paste("spatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.spatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/lattice_sim")
	setwd(paste("nonspatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.nonspatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/lattice_sim")
}

pdf(file="~/desktop/lattice_pops_spatial_vs_nonspatial.pdf",width=8,height=3,pointsize=11)
par(mfrow=c(1,3))
plot(prob.vec.nonspatial,pch=19,col="blue",ylab="Posterior Probability",
		ylim=c(min(prob.vec.nonspatial,prob.vec.spatial),
				max(prob.vec.nonspatial, prob.vec.spatial)),
		xaxt='n',cex=2,xlab="")
	axis(side=1,at=c(1,2,3,4),labels=c("K=1","K=2","K=3","K=4"))
	points(prob.vec.spatial,pch=19,col="green",cex=2)
	box(lwd=2)
	legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("nonspatial","spatial"))
plot(diff(prob.vec.nonspatial),pch=19,col="blue",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
plot(diff(prob.vec.spatial),pch=19,col="green",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
dev.off()

################################
#Barrier scenario
################################
setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim/")

#NONSPATIAL
load("nonspatial/k_1/ms.barrier.dataset.Robj")
load("nonspatial/k_1/data.list.Robj")

geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")

pdf(file="~/desktop/barrier_pops_nonspatial_samplecov.pdf",width=10,height=3,pointsize=12)
par(mfrow=c(1,4),oma=c(0.5,0.5,0.5,0.5))
#K_1
load("nonspatial/k_1/k_1_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="covariance",main="K = 1")
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
	box(lwd=2)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
#K_2
load("nonspatial/k_2/k_2_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 2")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)

mtext(text="geographic distance",side=1,adj=4,padj=3)

#K_3
load("nonspatial/k_3/k_3_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
#K_4
load("nonspatial/k_4/k_4_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
dev.off()



pdf(file="~/desktop/barrier_pops_nonspatial_adprop_piemaps.pdf",width=8,height=8,pointsize=12)
x.lim <- c(min(geo.coords[,1]) - 1, max(geo.coords[,1]) + 1)
y.lim <- c(min(geo.coords[,2]) - 1, max(geo.coords[,2]) + 1)
radius <- 3.5
par(mfrow=c(2,2),mar=c(1,1,1,1))
require(caroline)
#K_1
load("nonspatial/k_1/k_1_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=radius,
				xlab="",ylab="",main="K = 1",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim) 
box(lwd=2)
#K_2
load("nonspatial/k_2/k_2_output.Robj")
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
load("nonspatial/k_3/k_3_output.Robj")
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
load("nonspatial/k_4/k_4_output.Robj")
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




#SPATIAL
load("spatial/k_1/ms.barrier.dataset.Robj")
load("spatial/k_1/data.list.Robj")

geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")

pdf(file="~/desktop/barrier_pops_spatial_samplecov.pdf",width=10,height=3,pointsize=12)
par(mfrow=c(1,4),oma=c(0.5,0.5,0.5,0.5))
#K_1
load("spatial/k_1/k_1_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="covariance",main="K = 1")
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
	box(lwd=2)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
#K_2
load("spatial/k_2/k_2_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 2")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)

mtext(text="geographic distance",side=1,adj=4,padj=3)

#K_3
load("spatial/k_3/k_3_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
#K_4
load("spatial/k_4/k_4_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
dev.off()



pdf(file="~/desktop/barrier_pops_spatial_adprop_piemaps.pdf",width=8,height=8,pointsize=12)
x.lim <- c(min(geo.coords[,1]) - 1, max(geo.coords[,1]) + 1)
y.lim <- c(min(geo.coords[,2]) - 1, max(geo.coords[,2]) + 1)
radius <- 3.5
par(mfrow=c(2,2),mar=c(1,1,1,1))
require(caroline)
#K_1
load("spatial/k_1/k_1_output.Robj")
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=radius,
				xlab="",ylab="",main="K = 1",lty=1,density=NULL,xaxt='n',yaxt='n',xlim=x.lim,ylim=y.lim) 
box(lwd=2)
#K_2
load("spatial/k_2/k_2_output.Robj")
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
load("spatial/k_3/k_3_output.Robj")
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
load("spatial/k_4/k_4_output.Robj")
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



#Spatial vs. Nonspatial LnL comparison
prob.vec.spatial <- numeric(4)
prob.vec.nonspatial <- numeric(4)
for(i in 1:4){
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim")
	setwd(paste("spatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.spatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim")
	setwd(paste("nonspatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.nonspatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim")
}

pdf(file="~/desktop/barrier_pops_spatial_vs_nonspatial.pdf",width=8,height=3,pointsize=11)
par(mfrow=c(1,3))
plot(prob.vec.nonspatial,pch=19,col="blue",ylab="Posterior Probability",
		ylim=c(min(prob.vec.nonspatial,prob.vec.spatial),
				max(prob.vec.nonspatial, prob.vec.spatial)),
		xaxt='n',cex=2,xlab="")
	axis(side=1,at=c(1,2,3,4),labels=c("K=1","K=2","K=3","K=4"))
	points(prob.vec.spatial,pch=19,col="green",cex=2)
	box(lwd=2)
	legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("nonspatial","spatial"))
plot(diff(prob.vec.nonspatial),pch=19,col="blue",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
plot(diff(prob.vec.spatial),pch=19,col="green",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
dev.off()


#SPATIAL with ONE IBD RATE
load("barrier_sim_1ibd/k_2/data.list.Robj")
load("barrier_sim_1ibd/k_2/ms.barrier.dataset.Robj")
geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")

pdf(file="~/desktop/barrier_1_IBD_pops_spatial_samplecov.pdf",width=7.5,height=3,pointsize=12)
par(mfrow=c(1,3),oma=c(0.5,0.5,0.5,0.5))
#K_2
load("barrier_sim_1ibd/k_2/k_2_1ibd_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 2")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)

mtext(text="geographic distance",side=1,adj=17,padj=3.5)

#K_3
load("barrier_sim_1ibd/k_3/k_3_1ibd_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 3")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
#K_4
load("barrier_sim_1ibd/k_4/k_4_1ibd_output.Robj")
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="",
		ylab="",main="K = 4")
	box(lwd=2)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
dev.off()



pdf(file="~/desktop/barrier_pops_spatial_adprop_piemaps.pdf",width=8,height=8,pointsize=12)
x.lim <- c(min(geo.coords[,1]) - 1, max(geo.coords[,1]) + 1)
y.lim <- c(min(geo.coords[,2]) - 1, max(geo.coords[,2]) + 1)
radius <- 3.5
par(mfrow=c(2,2),mar=c(1,1,1,1))
require(caroline)

#K_2
load("barrier_sim_1ibd/k_2/k_2_1ibd_output.Robj")
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
load("barrier_sim_1ibd/k_3/k_3_1ibd_output.Robj")
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
load("barrier_sim_1ibd/k_4/k_4_1ibd_output.Robj")
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



#Spatial vs. Nonspatial LnL comparison
prob.vec.spatial <- numeric(4)
prob.vec.spatial.1ibd <- numeric(4)
prob.vec.nonspatial <- numeric(4)
for(i in 1:4){
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim")
	setwd(paste("spatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.spatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim")
	setwd(paste("nonspatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.nonspatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim")
	if(i > 1){
	setwd(paste("barrier_sim_1ibd/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.spatial.1ibd[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	}
	setwd("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/barrier_sim")
}

pdf(file="~/desktop/barrier_pops_spatial_vs_nonspatial_vs1ratespatial.pdf",width=8,height=3,pointsize=11)
par(mfrow=c(1,4))
plot(prob.vec.nonspatial,pch=19,col="blue",ylab="Posterior Probability",
		ylim=c(min(prob.vec.nonspatial,prob.vec.spatial),
				max(prob.vec.nonspatial, prob.vec.spatial)),
		xaxt='n',cex=2,xlab="")
	axis(side=1,at=c(1,2,3,4),labels=c("K=1","K=2","K=3","K=4"))
	points(prob.vec.spatial,pch=19,col="green",cex=2)
	points(prob.vec.spatial.1ibd,pch=19,col="purple",cex=2)
	box(lwd=2)
	legend(x="bottomright",pch=19,col=c("blue","green","purple"),legend=c("nonspatial","spatial","one-rate spatial"))
plot(diff(prob.vec.nonspatial),pch=19,col="blue",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
plot(diff(prob.vec.spatial),pch=19,col="green",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(1,2,3),labels=c("K=2 vs. K=1","K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
plot(c(1,2,3),c(NA,diff(prob.vec.spatial.1ibd[2:4])),pch=19,col="purple",ylab="Difference in Posterior Probability",cex=2,xaxt='n',xlab="")
	axis(side=1,at=c(2,3),labels=c("K=3 vs. K=2","K=4 vs. K=3"))
	box(lwd=2)
dev.off()