k <- 20
n.loci <- 1e4
epsilon <- rbeta(n.loci,10,10)
coords <- cbind(runif(20),runif(20))
geoDist <- fields::rdist(coords)
cov.func <- function(a0,a1,a2,geoDist){
	geoCov <- a0 * exp(-(a1*geoDist)^a2)
	return(geoCov)
}
geoCov <- cov.func(0.01,1.9,0.7,geoDist)
#plot(geoDist,geoCov)
require(MASS)
true.freqs <- matrix(NA,nrow=k,ncol=n.loci)
for(i in 1:n.loci){
	true.freqs[,i] <- mvrnorm(1,mu=rep(epsilon[i],k),Sigma=geoCov)
}
bounded.freqs <- true.freqs
bounded.freqs[which(bounded.freqs < 0)] <- 0
bounded.freqs[which(bounded.freqs > 1)] <- 1
sample.sizes <- sample(2:10,20,replace=TRUE)
sample.size.mat <- matrix(sample.sizes,nrow=k,ncol=n.loci)
obs.counts <- matrix(rbinom(k*n.loci,sample.size.mat,bounded.freqs),nrow=k,ncol=n.loci)
obs.freqs <- obs.counts/sample.size.mat
#plot(bounded.freqs,obs.freqs,pch=20,col=adjustcolor(1,0.01)) ; abline(0,1,col="red")

true.freq.cov <- cov(t(true.freqs))
bounded.freq.cov <- cov(t(bounded.freqs))
obs.cov <- cov(t(obs.freqs))

var.eps <- var(epsilon)
binVar <- mean(epsilon*(1-epsilon))

par(mfrow=c(1,3))
	plot(geoDist,geoCov,pch=20,ylim=c(range(c(geoCov,true.freq.cov,bounded.freq.cov,obs.cov))))
		points(geoDist,true.freq.cov,pch=20,col="red")
		points(geoDist,bounded.freq.cov,pch=20,col="purple")
		points(geoDist,obs.cov,pch=20,col="blue")
	legend(x="right",pch=19,col=c("black","red","purple","blue"),
			legend=c("geoCov","true freq cov","bounded freq cov","observed freq cov"))
	plot(geoDist,geoCov,pch=20,ylim=c(range(c(geoCov,true.freq.cov-var.eps,bounded.freq.cov-var.eps,obs.cov-var.eps-diag(binVar/sample.sizes)))))
		points(geoDist,true.freq.cov-var.eps,pch=20,col="red")
		points(geoDist,bounded.freq.cov-var.eps,pch=20,col="purple")
		points(geoDist,obs.cov-var.eps-diag(binVar/sample.sizes),pch=20,col="blue")
	plot(bounded.freq.cov-var.eps,obs.cov-var.eps-diag(binVar/sample.sizes))
		abline(0,1,col="red")
if(FALSE){
	plot(geoDist,geoCov,pch=20,ylim=c(range(c(geoCov,true.freq.cov-var.eps,bounded.freq.cov-var.eps,obs.cov-var.eps))))
		points(geoDist,true.freq.cov-var.eps,pch=20,col="red")
		points(geoDist,bounded.freq.cov-var.eps+2*(mean(epsilon*(1-epsilon)))^2,pch=20,col="purple")
		points(geoDist,obs.cov-var.eps+2*(mean(epsilon*(1-epsilon)))^2,pch=20,col="blue")	
}
	plot(geoDist,obs.cov,pch=20,ylim=range(c(obs.cov,c(var.eps+geoCov+diag(binVar/sample.sizes)))))
		points(geoDist,var.eps + (geoCov+diag(binVar/sample.sizes)),pch=20,col="red")
		points(geoDist,var.eps + binVar2*(geoCov+diag(1/sample.sizes)),pch=20,col="green")


	
	
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
freq.mat2 <- matrix(NA,nrow=100,ncol=100000)
for(i in 1:100){
	freq.mat2[i,] <- runif(100000)
	tmp.counts <- rbinom(100000,sample.size.mat[i,1],prob=freq.mat2[i,])
	var.vec[i] <- var(tmp.counts/sample.size.mat[i,1])
}

test2 <- unlist(lapply(1:100,function(i){mean(freq.mat2[i,] * (1-freq.mat2[i,]))/sample.size.mat[i,1]}))
plot(test2+var(freq.mat2[1,]),var.vec) ; abline(0,1,col="red")


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
