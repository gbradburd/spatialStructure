load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/globe/globetrotter_data.Robj")
setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/globe")

globe.coords <- cbind(globetrotter.long,globetrotter.lat)

random.switcharoo <- function(x){
	swap <- ifelse(runif(1) < 0.5,1,0)
	if(swap){
		x <- 1-x
	}
	return(x)
}

globe.freqs <- globetrotter.counts/globetrotter.sample.sizes
for(i in 1:ncol(globe.freqs)){
	globe.freqs[,i] <- random.switcharoo(globe.freqs[,i])
}

globe.sample.covariance <- cov(t(globe.freqs),use="pairwise.complete.obs")
row.names(globe.coords) <- row.names(globe.sample.covariance)
globe.data <- list("globe.sample.covariance" = globe.sample.covariance,
					"globe.coords" = globe.coords)
save(globe.data,file="globe.spStr.dataset.Robj")