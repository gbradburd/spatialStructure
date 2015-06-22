load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/metadata/Poplar_freqs.Robj")

pop.metadata <- pop.metadata[-1,]
ind.freqs <- poplar.freqs
drainages <- unique(pop.metadata$Drainage.Location.name)
drainage.freqs <- matrix(NA,nrow=length(drainages),ncol=ncol(poplar.freqs))
drainage.long <- numeric(length(drainages))
drainage.lat <- numeric(length(drainages))

for(i in 1:length(drainages)){
	ind.drainage.indices <- which(pop.metadata$Drainage.Location.name==drainages[i])
	drainage.freqs[i,] <- colSums(ind.freqs[ind.drainage.indices,,drop=FALSE],na.rm=TRUE)/length(ind.drainage.indices)
	drainage.long[i] <- mean(as.numeric(pop.metadata$Longitude[ind.drainage.indices]))
	drainage.lat[i] <- mean(as.numeric(pop.metadata$Latitude[ind.drainage.indices]))
}

drainage.coords <- cbind(drainage.long,drainage.lat)
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







