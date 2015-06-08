load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers/warbler_ind_dataset.Robj")
setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers")

warb.counts <- warbler.ind.allele.counts
warb.sample.sizes <- warbler.ind.sample.sizes

random.switcharoo <- function(x){
	swap <- ifelse(runif(1) < 0.5,
					1,
					0)
	if(swap){
		x <- 1-x
	}
	return(x)
}

for(i in 1:ncol(warb.counts)){
	warb.counts[,i] <- random.switcharoo(warb.counts[,i])
}

warb.sample.covariance <- cov(t(warb.counts/warb.sample.sizes))

warb.data <- list("warb.sample.covariance" = warb.sample.covariance,
					"warbler.ind.coords" = warbler.ind.coords)
save(warb.data,file="warbler.ind.dataset.Robj")