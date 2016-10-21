load("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/metadata/Poplar_freqs.Robj")

pop.metadata <- pop.metadata[-1,]
ind.freqs <- poplar.freqs
drainages <- unique(pop.metadata$Drainage.Location.name)
drainage.freqs <- matrix(NA,nrow=length(drainages),ncol=ncol(poplar.freqs))
drainage.long <- numeric(length(drainages))
drainage.lat <- numeric(length(drainages))
drainage.sp <- numeric(length(drainages))
drainage.sample.sizes <- numeric(length(drainages))

for(i in 1:length(drainages)){
	ind.drainage.indices <- which(pop.metadata$Drainage.Location.name==drainages[i])
	drainage.sample.sizes[i] <- 2 * mean(apply(ind.freqs[ind.drainage.indices,,drop=FALSE],2,function(x){length(which(!is.na(x)))})) #length(ind.drainage.indices) 
	drainage.freqs[i,] <- colSums(ind.freqs[ind.drainage.indices,,drop=FALSE],na.rm=TRUE)/apply(ind.freqs[ind.drainage.indices,,drop=FALSE],2,function(x){length(which(!is.na(x)))}) #colSums(ind.freqs[ind.drainage.indices,,drop=FALSE],na.rm=TRUE)/length(ind.drainage.indices)
	drainage.long[i] <- mean(as.numeric(pop.metadata$Longitude[ind.drainage.indices]))
	drainage.lat[i] <- mean(as.numeric(pop.metadata$Latitude[ind.drainage.indices]))
	if(length(unique(pop.metadata$Species[ind.drainage.indices])) < 2){
		drainage.sp[i] <- unique(pop.metadata$Species[ind.drainage.indices])
	}	else {
		drainage.sp[i] <- "mixed"
	}
}

drainage.coords <- cbind(drainage.long,drainage.lat)

poplar.data <- list("freqs" = drainage.freqs,
					"coords" = drainage.coords,
					"sp.ID" = drainage.sp,
					"sample.sizes" = drainage.sample.sizes)

trichocarpa.data <- list("freqs" = drainage.freqs[drainage.sp=="Populus trichocarpa",],
						 "coords" = drainage.coords[drainage.sp=="Populus trichocarpa",],
						 "sp.ID" = drainage.sp[drainage.sp=="Populus trichocarpa"],
						 "sample.sizes" = drainage.sample.sizes[drainage.sp=="Populus trichocarpa"])

save(poplar.data,file="poplar.spStr.dataset.Robj")
save(trichocarpa.data,file="trichocarpa.spStr.dataset.Robj")


if(FALSE){
	require(maps)
	map(xlim=range(drainage.coords[,1]),ylim=range(drainage.coords[,2]))
	points(drainage.coords)
}

random.switcharoo <- function(x){
	swap <- ifelse(runif(1) < 0.5,1,0)
	if(swap){
		x <- 1-x
	}
	return(x)
}

for(i in 1:ncol(poplar.freqs)){
	ind.freqs[,i] <- random.switcharoo(ind.freqs[,i])
	drainage.freqs[,i] <- random.switcharoo(drainage.freqs[,i])
}

drainage.cov <- cov(t(drainage.freqs),use="pairwise.complete.obs")
poplar.data <- list("geo.coords" = drainage.coords,
					"sample.cov" = drainage.cov,
					"n.loci" = ncol(poplar.freqs))

test.data <- list("drainage.coords" = drainage.coords,
					"drainage.cov" =  drainage.cov,
					"drainage.freqs" = drainage.freqs,
					"drainage.sample.sizes" = drainage.sample.sizes,
					"drainage.sp" = drainage.sp)
save(test.data,file="~/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/sampling.effect/poplar.data.Robj")

poplar.drainage.freq.data <- list("drainage.cov" = drainage.cov,
								  "geo.coords" = drainage.coords,
								  "drainage.freqs" = drainage.freqs,
								  "n.loci" = ncol(drainage.freqs))

mean.tri.freqs <- colMeans(drainage.freqs[which(drainage.sp=="Populus trichocarpa"),])
trichocarpa.data <- list("geo.coords" = drainage.coords[which(drainage.sp=="Populus trichocarpa"),],
							"sample.cov" = cov(t(drainage.freqs[which(drainage.sp=="Populus trichocarpa"),]),use="pairwise.complete.obs"),
							"n.loci" = ncol(poplar.freqs),
							"sample.sizes" = drainage.sample.sizes[which(drainage.sp=="Populus trichocarpa")],
							"binom.var" = mean(mean.tri.freqs*(1-mean.tri.freqs)))
save(poplar.drainage.freq.data,file="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/poplar.drainage.freq.data.Robj")
save(trichocarpa.data,file="~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/trichocarpa_dataset.Robj")
save(poplar.data,file="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/poplar.spStr.dataset.Robj")

if(FALSE){
ind.cov <- cov(t(ind.freqs),use="pairwise.complete.obs")
drainage.cov <- cov(t(drainage.freqs),use="pairwise.complete.obs")
ind.dist <- fields::rdist(cbind(pop.metadata$Longitude,pop.metadata$Latitude))
drainage.dist <- fields::rdist(drainage.coords)
ind.ind.mat <- upper.tri(ind.dist,diag=TRUE)
drainage.ind.mat <- upper.tri(drainage.dist,diag=TRUE)
quartz(width=10,height=5)
par(mfrow=c(1,2))
plot(drainage.dist[drainage.ind.mat], drainage.cov[drainage.ind.mat])
plot(ind.dist[ind.ind.mat], ind.cov[ind.ind.mat])

max.neighborhood.size <- 0.5
geo.coords <- matrix(as.numeric(geo.coords),nrow=nrow(geo.coords),ncol=2)
D <- fields::rdist.earth(geo.coords,miles=TRUE)

neighbors <- vector("list",length = nrow(geo.coords))
for(i in 1:nrow(geo.coords)){
	tmp.D <- fields::rdist.earth(geo.coords[i,,drop=FALSE],geo.coords,miles=TRUE)
	neighbors[[i]] <- which(fields::rdist(geo.coords[i,,drop=FALSE],geo.coords) < max.neighborhood.size)
}


flags <- c()
for(i in 1:length(neighbors)){
	tmp.neighbs <- 	neighbors[[i]]
	neighbs.set <- unique(unlist(lapply(seq_along(tmp.neighbs),function(i){neighbors[[tmp.neighbs[i]]]})))
	if(!identical(neighbs.set,tmp.neighbs)){
		flags <- c(flags,i)
	}	
}
flags
geo.coords[flags,]
max(fields::rdist(geo.coords[flags,]))

groups <- vector("list")
n.groups <- 1
for(i in 1:length(neighbors)){
	tmp.neighbs <- 	neighbors[[i]]
	neighbs.set <- unique(unlist(lapply(seq_along(tmp.neighbs),function(i){neighbors[[tmp.neighbs[i]]]})))
	if(is.null(unlist(lapply(groups,identical,neighbs.set))) | !any(unlist(lapply(groups,identical,neighbs.set)))){
		groups[[n.groups]] <- neighbs.set
		n.groups <- n.groups + 1
	}
}
} 







