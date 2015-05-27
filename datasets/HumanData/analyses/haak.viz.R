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
	plot.sample.ellipses2(metadata,radii=1)
	# plot.sample.ellipses(metadata,0.2,0.9)
	view3d(theta=0,phi=0)
		M1 <- par3d("userMatrix")
	view3d(theta=0,phi=-86.1)
		M2 <- par3d("userMatrix")
	view3d(theta=0,phi=-57.5)
		M3 <- par3d("userMatrix")
	view3d(theta=0,phi=-34.5)
		M4 <- par3d("userMatrix")

	view.matrix.list1 <- c(list(M1),lapply(1:3,function(i){rotate3d(M1,-0.5*i,1,0,0)}))
	view.matrix.list2 <- c(list(view.matrix.list1[[4]]),lapply(1:2,function(i){rotate3d(M2,0.25*i,1,0,0)}))
	view.matrix.list3 <- c(list(view.matrix.list1[[4]]),lapply(1:2,function(i){rotate3d(M2,0.45*i,1,0,0)}))
spin.x1 <- par3dinterp(userMatrix=view.matrix.list1,method="linear",extrapolate="constant")
spin.x2 <- par3dinterp(userMatrix=view.matrix.list2,method="linear",extrapolate="constant")
spin.x3 <- par3dinterp(userMatrix=view.matrix.list3,method="linear",extrapolate="constant")
spin.z <- spin3d(axis=c(0,0,-1))
plotfunc1 <- function(time,base1,base2){
	# recover()
	# par3d(skipRedraw = TRUE)
	if(time < 3){
		userMatrix <- spin.x1(time)
	} else if (time >= 3 && time <= 3.5){
		userMatrix <- base1
	} else if (time > 3.5 && time <= 5.5){
		userMatrix <- spin.x2(time-3.5)
	} else if (time > 5.5){
		userMatrix <- spin.z(time-5.5,base=base2)
	}
}
plotfunc2 <- function(time,base1,base2){
	# recover()
	# par3d(skipRedraw = TRUE)
	if(time < 3){
		userMatrix <- spin.x1(time)
	} else if (time >= 3 && time <= 3.5){
		userMatrix <- base1
	} else if (time > 3.5 && time <= 5.5){
		userMatrix <- spin.x3(time-3.5)
	} else if (time > 5.5){
		userMatrix <- spin.z(time-5.5,base=base2)
	}
}

play3d(plotfunc1,duration=17.5,base1=M2,base2=M3)

# play3d(par3dinterp(userMatrix=view.matrix.list,method="linear",extrapolate="constant"),duration=20)
movie3d(plotfunc1,base1=M2,base2=M3,duration=17.5,dir="~/Desktop/test_spinny",clean=TRUE,convert=TRUE,fps=15)


tmp.metadata <- metadata
tmp.metadata$time <- tmp.metadata$time + 101
pies3d(x=pie.list,n.clusters=length(pie.list[[1]]),locations=tmp.metadata, radii=1)
play3d(plotfunc2,duration=17.5,base1=M2,base2=M4)
movie3d(plotfunc2,base1=M2,base2=M4,duration=17.5,dir="~/Desktop/test_spinny_pies",clean=TRUE,convert=TRUE,fps=15)

pies3d.2(x=pie.list,n.clusters=length(pie.list[[1]]),locations=metadata, radii=1)
view.matrix.list2 <- c(list(view.matrix.list1[[4]]),lapply(1:2,function(i){rotate3d(M2,0.25*i,1,0,0)}))
spin.x2.pies <- par3dinterp(userMatrix=view.matrix.list2,method="linear",extrapolate="constant")
spin.z <- spin3d(axis=c(0,0,-1))
plotfunc <- function(time,base1,base2){
	# recover()
	# par3d(skipRedraw = TRUE)
	if(time < 3){
		userMatrix <- spin.x1(time)
	} else if (time >= 3 && time <= 3.5){
		userMatrix <- base1
	} else if (time > 3.5 && time <= 5.5){
		userMatrix <- spin.x2(time-3.5)
	} else if (time > 5.5){
		userMatrix <- spin.z(time-5.5,base=base2)
	}
}
tmp.metadata <- metadata
tmp.metadata$time <- tmp.metadata$time + 101
pies3d(x=pie.list,n.clusters=length(pie.list[[1]]),locations=tmp.metadata, radii=1)

movie3d(plotfunc,base1=M2,base2=M3,duration=17.5,dir="~/Desktop/test_spinny_pies",clean=TRUE,convert=TRUE,fps=15)


make.structure.plot(data.list,super.list,cluster.order=c(2,1,3),sample.names=row.names(data.list$sample.covariance),sort.by=2)



if(FALSE){
	M3 <- M2
		M3[1,1] <- -M3[1,1]
		M3[2,2] <- -M3[2,2]
		M3[3,2] <- -M3[3,2]
	view.matrix.list <- c(list(M1),
							lapply(1:3,function(i){rotate3d(M1,-0.5*i,1,0,0)}),
							list(M2),
							lapply(1:8,function(i){rotate3d(M2,-(i*pi)/10,0,0,1)}),
							list(M3),
							lapply(1:10,function(i){rotate3d(M3,-(i*pi)/10,0,0,1)}))
}