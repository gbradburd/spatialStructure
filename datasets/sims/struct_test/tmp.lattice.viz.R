plot(data.list$geo.dist,data.list$sample.covariance)
points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,col="red",pch=20)



plot(super.list$parameter.list$admix.proportions[,1])

require(caroline)
sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
all.colors <- c("blue","red","green","yellow","purple","brown")
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list2 <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list2) <- sample.names
#pdf(file="admix_prop_map.pdf",width=10,height=5)
quartz(width=6,height=6)
pies(pie.list2,x0=spacemix.dataset$population.coordinates[,1],
				y0=spacemix.dataset$population.coordinates[,2],
				color.table=color.tab,border="black",radii=3,
				xlab="",ylab="",main="estimated",lty=1,density=NULL)
				
plot(super.list$output.list$posterior.prob[200:1000])

par(mfrow=c(1,3))
plot(data.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_1$covariance)
plot(data.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_2$covariance)
plot(data.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_3$covariance)

prob.vec.spatial <- numeric(4)
prob.vec.nonspatial <- numeric(4)
for(i in 1:4){
	setwd(paste("~/Desktop/struct_test/lattice_sim/spatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.spatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
	setwd(paste("~/Desktop/struct_test/lattice_sim/nonspatial/k_",i,sep=""))
		load(list.files(pattern="output"))
		prob.vec.nonspatial[i] <- max(super.list$output.list$posterior.prob,na.rm=TRUE)
}

plot(prob.vec.nonspatial,pch=19,col=1,ylab="posterior.prob",
		ylim=c(min(prob.vec.nonspatial,prob.vec.spatial),
				max(prob.vec.nonspatial, prob.vec.spatial)))
	points(prob.vec.spatial,pch=19,col=2)

source('~/Desktop/struct_test/lattice_sim/spatial/k_3/SpaceTimeStructureMix.R')	


new.cluster.list <- super.list$parameter.list$cluster.list
new.cluster.list$Cluster_1$covariance.params$cov.par2 <- 1.5e4
new.cluster.list$Cluster_1$cluster.mean <- 0.01
new.cluster.list$Cluster_1$covariance <- cluster.covariance(data.list$geo.dist,data.list$time.dist,new.cluster.list$Cluster_1$covariance.params)
new.adcov <- admixed.covariance(new.cluster.list,3,super.list$parameter.list$shared.mean,super.list$parameter.list$nuggets)
quartz(width=11,height=6)
par(mfrow=c(1,2))
plot(super.list$parameter.list$admixed.covariance,new.adcov) ; abline(0,1,col="red")
plot(data.list$geo.dist,data.list$sample.covariance,pch=19)
points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,col="red",pch=20)
points(data.list$geo.dist, new.adcov,col="blue",pch=20)


