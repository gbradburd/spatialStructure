bear.geno <- read.table("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/bears/bear_data.raw",header=TRUE,stringsAsFactors=FALSE)
sample.names <- bear.geno[,2]
bear.geno <- as.matrix(bear.geno[,-c(1:6)])
bear.freqs <- bear.geno/2
sample.cov <- cov(t(bear.freqs),use="pairwise.complete.obs")

metadata.raw <- read.csv("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/bears/metadata/Puckett-UamerSampleCoordinates.csv",stringsAsFactors=FALSE)
metadata <- as.matrix(metadata.raw)[match(sample.names, metadata.raw$Sample_ID),]

pc.obj <- eigen(sample.cov)$vectors
geo.factor <- (as.numeric(metadata[,4]) + abs(min(as.numeric(metadata[,4]))))
geo.factor <- geo.factor/(2*max(geo.factor))
geo.col <- hsv(h=geo.factor,s=1,v=1)
gen.col <- hsv(h=pc.obj[,2] + abs(min(pc.obj[,2])),s=1,v=1)
coords <- cbind(as.numeric(metadata[,4]),as.numeric(metadata[,3]))
layout(matrix(c(1,2,3,3),ncol=2,nrow=2,byrow=TRUE))
map(regions=c("Canada","USA"),xlim=c(-170,-50),ylim=c(25,80)) ; points(metadata[,4],metadata[,3],pch=20,col=gen.col)
plot(pc.obj[,2],pc.obj[,3],pch=20,col=geo.col)
plot(fields::rdist.earth(coords),sample.cov,pch=20)

col.mat1 <- matrix(NA,nrow=length(geo.col),ncol=length(geo.col))
col.mat2 <- matrix(NA,nrow=length(geo.col),ncol=length(geo.col))
for(i in 1:length(geo.col)){
	for(j in 1:length(geo.col)){
		col.mat1[i,j] <- geo.col[i]
		col.mat2[i,j] <- geo.col[j]
	}
}
par(mfrow=c(1,2))
	map(regions=c("Canada","USA"),xlim=c(-170,-50),ylim=c(25,80)) ; points(metadata[,4],metadata[,3],pch=20,col=gen.col)
	plot(fields::rdist.earth(coords),sample.cov,pch=21,bg=col.mat1,col=col.mat2)

bear.dataset <- list("sample.coords" = coords,
					 "sample.freqs" = bear.freqs,
					 "sample.sizes" = rep(2,nrow(coords)))

save(bear.dataset,file="bear.dataset.Robj")

# tmp.metadata <- metadata.raw
# tmp.metadata[["match_in_geno"]] <- rep("no.match",nrow(metadata.raw))
# in.geno.not.meta <- c()
# for(i in 1:length(sample.names)){
	# if(length(which(tmp.metadata$Sample_ID==sample.names[i]))==0){
		# in.geno.not.meta <- c(in.geno.not.meta,sample.names[i])	
	# }
	# tmp.metadata$match_in_geno[which(tmp.metadata$Sample_ID==sample.names[i])] <- sample.names[i]
# }

