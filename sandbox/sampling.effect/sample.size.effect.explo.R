var.mat <- matrix(NA,nrow=100,ncol=100)
sample.size.mat <- matrix(seq(1,100,length.out=100),nrow=100,ncol=100)
freq.mat <- matrix(seq(1,100,length.out=100)/100,nrow=100,ncol=100,byrow=TRUE)
for(i in 1:100){
	for(j in 1:100){
		tmp.counts <- rbinom(1e5,sample.size.mat[i,j],prob=freq.mat[i,j])
		var.mat[i,j] <- var(tmp.counts/sample.size.mat[i,j])
	}
	message(i)
}

test <- (freq.mat*(1-freq.mat)) / sample.size.mat
plot(c(test),c(var.mat)) ; abline(0,1,col="red")

var.vec <- numeric(100)
freq.mat <- matrix(NA,nrow=100,ncol=100000)
for(i in 1:100){
	freq.mat[i,] <- runif(100000)
	tmp.counts <- rbinom(100000,sample.size.mat[i,1],prob=freq.mat[i,])
	var.vec[i] <- var(tmp.counts/sample.size.mat[i,1])
}

test2 <- unlist(lapply(1:100,function(i){mean(freq.mat[i,] * (1-freq.mat[i,]))/sample.size.mat[i,1]}))
plot(test2+var(freq.mat[1,]),var.vec) ; abline(0,1,col="red")


load("/Users/gburd/Desktop/k_1b/k1b_tri_output.Robj")
load("/Users/gburd/Desktop/k_1b/data.list.Robj")
source('~/Desktop/k_1b/SpaceTimeStructureMix.R')
plot(data.list$geo.dist,data.list$sample.covariance)
plot(data.list$sample.covariance,super.list$parameter.list$admixed.covariance/mean(norms))
	abline(0,1,col="red")
#points(data.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_1$covariance,col="darkgreen")
#points(rep(1,length(super.list$parameter.list$nuggets)),super.list$parameter.list$nuggets,col="pink")
#abline(h=super.list$parameter.list$shared.mean,col="blue")
#abline(h=super.list$parameter.list$cluster.list$Cluster_1$cluster.mean,col="orange")
points(data.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_1$covariance + 
										super.list$parameter.list$shared.mean + 
										super.list$parameter.list$cluster.list$Cluster_1$cluster.mean + 
										diag(super.list$parameter.list$nuggets),
										col="darkgreen")
points(data.list$geo.dist,super.list$parameter.list$cluster.list$Cluster_1$covariance + 
										min(data.list$sample.covariance) + 
										super.list$parameter.list$cluster.list$Cluster_1$cluster.mean + 
										diag(super.list$parameter.list$nuggets),
										col="red")



test <- admixed.covariance(super.list$parameter.list$cluster.list,
							1,
							super.list$parameter.list$shared.mean,
							super.list$parameter.list$nuggets)
test2 <- super.list$parameter.list$cluster.list$Cluster_1$covariance + 
			super.list$parameter.list$cluster.list$Cluster_1$cluster.mean + 
			super.list$parameter.list$shared.mean + 
			diag(super.list$parameter.list$nuggets)
identical(test,super.list$parameter.list$admixed.covariance)


load("~/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/sampling.effect/poplar.data.Robj")
load("~/Desktop/k_1b/k1b_tri_output.Robj")
tri.pops <- which(test.data$drainage.sp == "Populus trichocarpa")
mean.freqs <- colMeans(test.data$drainage.freqs[tri.pops,])
normalizers <- mean.freqs*(1-mean.freqs)
mean(normalizers)

plot(var(mean.freqs) + mean(normalizers)/test.data$drainage.sample.sizes[tri.pops],
		unlist(lapply(1:length(tri.pops),function(i){var(test.data$drainage.freqs[tri.pops[i],])})))
abline(0,1,col="red")
	mean(normalizers)
