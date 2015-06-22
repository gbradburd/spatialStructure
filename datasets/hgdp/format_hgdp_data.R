load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/hgdp/hgdp.himalaya.bedassle.Robj")
setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/hgdp")

hgdp.coords <- cbind(DATA$Latitude,DATA$Longitude)

random.switcharoo <- function(x){
	swap <- ifelse(runif(1) < 0.5,1,0)
	if(swap){
		x <- 1-x
	}
	return(x)
}

hgdp.freqs <- counts/sample_sizes
for(i in 1:ncol(hgdp.freqs)){
	hgdp.freqs[,i] <- random.switcharoo(hgdp.freqs[,i])
}

hgdp.sample.covariance <- cov(t(hgdp.freqs),use="pairwise.complete.obs")
row.names(hgdp.coords) <- DATA$Population
hgdp.data <- list("hgdp.sample.covariance" = hgdp.sample.covariance,
					"hgdp.coords" = hgdp.coords)
save(hgdp.data,file="hgdp.spStr.dataset.Robj")


#to fix Latitude/Longitude switch-mistake
for(i in 1:10){
	dir <- paste("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/hgdp/hgdp_analyses/spatial/k_",i,sep="")
	setwd(dir)
	load("data.list.Robj")
		data.list$geo.coords <- cbind(data.list$geo.coords[,2],data.list$geo.coords[,1])
		save(data.list,file="data.list.Robj")
	load("hgdp.spStr.dataset.Robj")
		hgdp.data$hgdp.coords <- cbind(hgdp.data$hgdp.coords[,2],hgdp.data$hgdp.coords[,1])
		save(hgdp.data,file="hgdp.spStr.dataset.Robj")
}
#nonspatial
for(i in 1:10){
	dir <- paste("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/hgdp/hgdp_analyses/nonspatial/k_",i,sep="")
	setwd(dir)
	load("data.list.Robj")
		data.list$geo.coords <- cbind(data.list$geo.coords[,2],data.list$geo.coords[,1])
		save(data.list,file="data.list.Robj")
	load("hgdp.spStr.dataset.Robj")
		hgdp.data$hgdp.coords <- cbind(hgdp.data$hgdp.coords[,2],hgdp.data$hgdp.coords[,1])
		save(hgdp.data,file="hgdp.spStr.dataset.Robj")
}