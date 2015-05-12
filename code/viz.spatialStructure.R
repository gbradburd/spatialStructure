#MCMC Plotting Functions
if(FALSE){
	load(list.files(pattern="output"))
	load("data.list.Robj")
}
geo.coords <- spacemix.dataset$population.coordinates
all.colors <- c("blue","red","green","yellow","purple","brown")
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))

#EMPIRICAL
quartz(width=11,height=5)
par(mfrow=c(1,3))
plot(data.list$geo.dist,data.list$sample.covariance,
		xlab="geographic distance",
		ylab="covariance")
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red",cex=0.8)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("sample covariance","parametric estimate"))
plot(data.list$time.dist,data.list$sample.covariance,
		xlab="temporal distance",
		ylab="covariance")
	points(data.list$time.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red")
plot(data.list$sample.covariance,super.list$parameter.list$admixed.covariance)
	abline(0,1,col="red")

quartz(width=8,height=5)

par(mfrow=c(2,length(super.list$parameter.list$cluster.list)))
for(i in 1:length(super.list$parameter.list$cluster.list)){
	plot(data.list$geo.dist,
			super.list$parameter.list$cluster.list[[i]]$covariance + 
			super.list$parameter.list$cluster.list[[i]]$cluster.mean,
			pch=20,col=all.colors[i],xlab="geo distance",ylab="covariance")
}
for(i in 1:length(super.list$parameter.list$cluster.list)){
	plot(data.list$time.dist,
			super.list$parameter.list$cluster.list[[i]]$covariance + 
			super.list$parameter.list$cluster.list[[i]]$cluster.mean,
			pch=20,col=all.colors[i],xlab="time distance",ylab="covariance")
}

plot(super.list$output.list$likelihood[300:length(which(!is.na(super.list$output.list$likelihood)))],xlab="",ylab="",type='l',
	ylim=c(min(super.list$output.list$likelihood[300:length(which(!is.na(super.list$output.list$likelihood)))],
				super.list$output.list$posterior.prob[300:length(which(!is.na(super.list$output.list$likelihood)))]),
			max(super.list$output.list$likelihood[300:length(which(!is.na(super.list$output.list$likelihood)))],
				super.list$output.list$posterior.prob[300:length(which(!is.na(super.list$output.list$likelihood)))])))
points(super.list$output.list$posterior.prob[300:length(which(!is.na(super.list$output.list$likelihood)))],type='l',col="red")

quartz(width=11,height=4)
par(mfrow=c(1,4))
matplot(t(super.list$output.list$cluster.params$cov.par1),type='l',lty=1,ylab="cov.par1",col=all.colors[1:length(super.list$parameter.list$cluster.list)])
	mtext("cov.par1",side=1,padj=-11)
matplot(t(super.list$output.list$cluster.params$cov.par2),type='l',lty=1,ylab="cov.par2",col=all.colors[1:length(super.list$parameter.list$cluster.list)])
	mtext("cov.par2",side=1,padj=-11)
matplot(t(super.list$output.list$cluster.params$cov.par3),type='l',lty=1,ylab="cov.par3",col=all.colors[1:length(super.list$parameter.list$cluster.list)])
	mtext("cov.par3",side=1,padj=-11)
matplot(t(super.list$output.list$cluster.params$cluster.mean),type='l',lty=1,ylab="cluster.mean",col=all.colors[1:length(super.list$parameter.list$cluster.list)])
	mtext("cluster.mean",side=1,padj=-11)
	legend(x="topright",col=all.colors[1:length(super.list$parameter.list$cluster.list)],lty=1,legend=cluster.names)

# make.admix.prop.mat <- function(super.list){
	# # recover()
	# x <- length(which(!is.na(super.list$output.list$posterior.prob)))
	# admix.prop.mat <- matrix(NA,nrow=nrow(super.list$output.list$admix.proportions[[1]]),ncol=x)
	# for(i in 1:x){
		# admix.prop.mat[,i] <- super.list$output.list$admix.proportions[[i]][,1]
	# }
	# return(admix.prop.mat)
# }
# admix.prop.mat <- make.admix.prop.mat(super.list)

par(mfrow=c(3,1),mar=c(1,1,1,1))
matplot(t(super.list$output.list$nuggets),type='l')
# matplot(t(admix.prop.mat),type='l')
	# abline(h=sim.param.list$sim.admix.props[,1],lwd=0.5)
plot(super.list$output.list$shared.mean,type='l',
	xlab="sampled mcmc generations",
	ylab="shared mean parameter estimate")
	
	
quartz(width=11,height=4)
par(mfrow=c(1,4),mar=c(1,1,1,1))
matplot(t(super.list$output.list$acceptance.rates$cov.par1),type='l',xlab="",ylab="",col=all.colors[1:length(super.list$parameter.list$cluster.list)])
matplot(t(super.list$output.list$acceptance.rates$cov.par2),type='l',xlab="",ylab="",col=all.colors[1:length(super.list$parameter.list$cluster.list)])
matplot(t(super.list$output.list$acceptance.rates$cov.par3),type='l',xlab="",ylab="",col=all.colors[1:length(super.list$parameter.list$cluster.list)])
matplot(t(super.list$output.list$acceptance.rates$cluster.mean),type='l',xlab="",ylab="",col=all.colors[1:length(super.list$parameter.list$cluster.list)])

par(mfrow=c(3,1),mar=c(1,1,1,1))
matplot(t(super.list$output.list$acceptance.rates$admix.proportions),type='l',xlab="",ylab="")
matplot(t(super.list$output.list$acceptance.rates$nuggets),type='l',xlab="",ylab="")
plot(super.list$output.list$acceptance.rates$shared.mean,type='l',xlab="",ylab="")


require(caroline)
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names

quartz(width=5,height=5)
# par(new=FALSE)
pies(pie.list,x0=geo.coords[,1],
				y0=geo.coords[,2],
				color.table=color.tab,border="black",radii=3,
				xlab="",ylab="",main="estimated",lty=1,density=NULL)
	# points(x=50.70267, y=52.84767,pch=8,col=1)
box(lwd=2)
################################################################
################################################################
################################################################
################################################################
if(FALSE){
	load("test1.Robj")
	load("data.list.Robj")
	load("simdata.Robj")
}
#SIMS
quartz(width=11,height=5)
par(mfrow=c(1,3))
plot(sim.param.list$geo.dist,sim.param.list$sim.admixed.cov.mat,
		xlab="geographic distance",
		ylab="covariance")
	points(sim.param.list$geo.dist,super.list$parameter.list$admixed.covariance,col="red",pch=20,cex=0.8)
	legend(x="topright",pch=c(1,20),col=c(1,2),legend=c("truth","parametric estimate"))
plot(sim.param.list$sim.admix.props[,1],
		xlab="sample index",
		ylab="admixture proportion in Cluster 1 (of 2)")
	points(super.list$parameter.list$admix.proportions[,1],col="red",pch=20)
plot(sim.param.list$geo.dist, sim.param.list$sim.cluster.list$Cluster_1$covariance + sim.param.list$sim.cluster.list$Cluster_1$cluster.mean,col="darkblue",
		ylim=c(min(sim.param.list$sim.cluster.list$Cluster_1$covariance + sim.param.list$sim.cluster.list$Cluster_1$cluster.mean,
					sim.param.list$sim.cluster.list$Cluster_2$covariance + sim.param.list$sim.cluster.list$Cluster_2$cluster.mean),
					max(sim.param.list$sim.cluster.list$Cluster_1$covariance + sim.param.list$sim.cluster.list$Cluster_1$cluster.mean,
					sim.param.list$sim.cluster.list$Cluster_2$covariance + sim.param.list$sim.cluster.list$Cluster_2$cluster.mean)+0.5),
					xlab="geographic distance",
					ylab="spatial covariance within a cluster")
	points(sim.param.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_1$covariance + super.list$parameter.list$cluster.list$Cluster_1$cluster.mean,col="lightblue",pch=20)
points(sim.param.list$geo.dist,sim.param.list$sim.cluster.list$Cluster_2$covariance + sim.param.list$sim.cluster.list$Cluster_2$cluster.mean,col="darkgreen")
	points(sim.param.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_2$covariance + super.list$parameter.list$cluster.list$Cluster_2$cluster.mean,col="lightgreen",pch=20)
	legend(x="topright",pch=19,col=c("darkblue","lightblue","darkgreen","lightgreen"),legend=c("truth - Cluster1","parametric estimate - Cluster1","truth - Cluster2","parametric estimate - Cluster2"))

plot(sim.param.list$sim.admix.props[,1],
	super.list$parameter.list$admix.proportions[,1],
		xlab="sim admix proportions",
		ylab="admixture proportion in Cluster 1 (of 2)",
		xlim=c(0,1),ylim=c(0,1))
		abline(0,1,col="red",lty=2)
		abline(1,-1,col="red",lty=2)

plot(super.list$output.list$likelihood[300:length(which(!is.na(super.list$output.list$likelihood)))],xlab="",ylab="",type='l',
	ylim=c(min(super.list$output.list$likelihood[300:length(which(!is.na(super.list$output.list$likelihood)))],
				super.list$output.list$posterior.prob[300:length(which(!is.na(super.list$output.list$likelihood)))]),
			max(super.list$output.list$likelihood[300:length(which(!is.na(super.list$output.list$likelihood)))],
				super.list$output.list$posterior.prob[300:length(which(!is.na(super.list$output.list$likelihood)))])))
points(super.list$output.list$posterior.prob[300:length(which(!is.na(super.list$output.list$likelihood)))],type='l',col="red")


par(mfrow=c(2,2),mar=c(2,2,2,2))
matplot(t(super.list$output.list$cluster.params$cov.par1),type='l',lty=1,
			ylim=c(min(sim.param.list$sim.cov.par1_1,sim.param.list$sim.cov.par1_2,super.list$output.list$cluster.params$cov.par1,na.rm=TRUE),
					max(sim.param.list$sim.cov.par1_1,sim.param.list$sim.cov.par1_2,super.list$output.list$cluster.params$cov.par1,na.rm=TRUE)))
	mtext("cov.par1",side=1,padj=-11)
	abline(h=c(sim.param.list$sim.cov.par1_1,sim.param.list$sim.cov.par1_2),col=c(1,2),lty=2)
	legend(x="topleft",col=c(1,1,2,2),lty=c(1,2,1,2),legend=c("cluster1 - estimate","cluster1 - truth","cluster2 - estimate","cluster2 - truth"))
matplot(t(super.list$output.list$cluster.params$cov.par2),type='l',lty=1,
			ylim=c(min(sim.param.list$sim.cov.par2_1,sim.param.list$sim.cov.par2_2,super.list$output.list$cluster.params$cov.par2,na.rm=TRUE),
					max(sim.param.list$sim.cov.par2_1,sim.param.list$sim.cov.par2_2,super.list$output.list$cluster.params$cov.par2,na.rm=TRUE)))
	mtext("cov.par2",side=1,padj=-11)
	abline(h=c(sim.param.list$sim.cov.par2_1,sim.param.list$sim.cov.par2_2),col=c(1,2),lty=2)
matplot(t(super.list$output.list$cluster.params$cov.par3),type='l',lty=1,
			ylim=c(min(sim.param.list$sim.cov.par3_1,sim.param.list$sim.cov.par3_2,super.list$output.list$cluster.params$cov.par3,na.rm=TRUE),
					max(sim.param.list$sim.cov.par3_1,sim.param.list$sim.cov.par3_2,super.list$output.list$cluster.params$cov.par3,na.rm=TRUE)))
	mtext("cov.par3",side=1,padj=-11)
	abline(h=c(sim.param.list$sim.cov.par3_1,sim.param.list$sim.cov.par3_2),col=c(1,2),lty=2)
matplot(t(super.list$output.list$cluster.params$cluster.mean),type='l',lty=1,
			ylim=c(min(sim.param.list$sim.cluster.mean.1,sim.param.list$sim.cluster.mean.2,super.list$output.list$cluster.params$cluster.mean,na.rm=TRUE),
					max(sim.param.list$sim.cluster.mean.1,sim.param.list$sim.cluster.mean.2,super.list$output.list$cluster.params$cluster.mean,na.rm=TRUE)))
	mtext("cluster.mean",side=1,padj=-11)
	abline(h=c(sim.param.list$sim.cluster.mean.1,sim.param.list$sim.cluster.mean.2),col=c(1,2),lty=2)

make.admix.prop.mat <- function(super.list){
	# recover()
	x <- length(which(!is.na(super.list$output.list$posterior.prob)))
	admix.prop.mat <- matrix(NA,nrow=nrow(super.list$output.list$admix.proportions[[1]]),ncol=x)
	for(i in 1:x){
		admix.prop.mat[,i] <- super.list$output.list$admix.proportions[[i]][,1]
	}
	return(admix.prop.mat)
}
admix.prop.mat <- make.admix.prop.mat(super.list)

par(mfrow=c(3,1),mar=c(1,1,1,1))
matplot(t(super.list$output.list$nuggets),type='l')
	abline(h=sim.param.list$sim.nuggets,lwd=0.5)
matplot(t(admix.prop.mat),type='l')
	abline(h=sim.param.list$sim.admix.props[,1],lwd=0.5)
plot(super.list$output.list$shared.mean,type='l',
	xlab="sampled mcmc generations",
	ylab="shared mean parameter estimate")
	abline(h=sim.param.list$sim.shared.mean,lwd=0.5)
	
	
par(mfrow=c(3,1),mar=c(1,1,1,1))
matplot(t(super.list$output.list$acceptance.rates$cov.par1),type='l',xlab="",ylab="")
matplot(t(super.list$output.list$acceptance.rates$cov.par2),type='l',xlab="",ylab="")
matplot(t(super.list$output.list$acceptance.rates$cov.par3),type='l',xlab="",ylab="")

par(mfrow=c(3,1),mar=c(1,1,1,1))
matplot(t(super.list$output.list$acceptance.rates$admix.proportions),type='l',xlab="",ylab="")
matplot(t(super.list$output.list$acceptance.rates$nuggets),type='l',xlab="",ylab="")
matplot(t(super.list$output.list$acceptance.rates$cluster.mean),type='l',xlab="",ylab="")


require(caroline)
sample.names <- unlist(lapply(1:sim.param.list$k,function(i){paste("sample_",i,sep="")}))
pie.list <- lapply(1:sim.param.list$k,function(i){nv(sim.param.list$sim.admix.props[i,],c("Cluster_1","Cluster_2"))})
names(pie.list) <- sample.names
color.tab <- nv(c("blue","red"),c("Cluster_1","Cluster_2"))
pie.list2 <- lapply(1:sim.param.list$k,function(i){nv(super.list$parameter.list$admix.proportions[i,],c("Cluster_1","Cluster_2"))})
names(pie.list2) <- sample.names
#pdf(file="admix_prop_map.pdf",width=10,height=5)
quartz(width=10,height=5)
par(mfrow=c(1,2))
pies(pie.list,x0=sim.param.list$spatial.coords[,1],y0=sim.param.list$spatial.coords[,2],color.table=color.tab,border="black",radii=3,
		xlab="",ylab="",main="truth",lty=1,density=NULL)
pies(pie.list2,x0=sim.param.list$spatial.coords[,1],y0=sim.param.list$spatial.coords[,2],color.table=color.tab,border="black",radii=3,
		xlab="",ylab="",main="estimated",lty=1,density=NULL)
#dev.off()
