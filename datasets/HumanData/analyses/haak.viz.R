#visualize haak output

metadata <- read.table("human_sample_metadata.txt",header=TRUE,stringsAsFactors=FALSE)
load(list.files(pattern="data.list"))
load(list.files(pattern="output"))
all.colors <- c("blue","red","green","yellow","purple","brown")

par(mfrow=c(1,3))
plot(data.list$time.dist,data.list$sample.covariance)
	points(data.list$time.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red")

plot(data.list$geo.dist,data.list$sample.covariance)
	points(data.list$geo.dist,super.list$parameter.list$admixed.covariance,pch=20,col="red")

plot(data.list$sample.covariance,super.list$parameter.list$admixed.covariance)
	abline(0,1,col="red")

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

# plot3d(x=data.list$geo.dist,y=data.list$time.dist,z=super.list$parameter.list$cluster.list$Cluster_1$covariance)
# plot3d(x=data.list$geo.dist,y=data.list$time.dist,z=super.list$parameter.list$cluster.list$Cluster_2$covariance)


require(maps)
require(caroline)


sample.names <- unlist(lapply(1:data.list$n.ind,function(i){paste("sample_",i,sep="")}))
cluster.names <- unlist(lapply(1:super.list$model.options$n.clusters,function(i){paste("Cluster_",i,sep="")}))
color.tab <- nv(c(all.colors[1:super.list$model.options$n.clusters]),cluster.names)
pie.list <- lapply(1:data.list$n.ind,function(i){nv(super.list$parameter.list$admix.proportions[i,],cluster.names)})
names(pie.list) <- sample.names


quartz(width=5,height=5)
map(xlim=c(-30,60),ylim=c(5,70))
par(new=TRUE)
pies(pie.list,x0=metadata$lon,
				y0=metadata$lat,
				color.table=color.tab,border="black",radii=1.5,
				xlab="",ylab="",main="estimated",lty=1,density=NULL)
	# points(x=50.70267, y=52.84767,pch=8,col=1)
box(lwd=2)

source('../viz_funcs.R')
make.spatiotemporal.sampling.plot(metadata)
plot.sample.ellipses(metadata,scale=0.2,alpha=0.8)
pies3d(x=pie.list,n.clusters=length(pie.list[[1]]),locations=metadata, radii=1)


plot(super.list$parameter.list$nuggets,type='n')
	text(super.list$parameter.list$nuggets,row.names(data.list$sample.covariance))
make.structure.plot(data.list,super.list,cluster.order=c(2,1,3),sample.names=row.names(data.list$sample.covariance),sort.by=2)
