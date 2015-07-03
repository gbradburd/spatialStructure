setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses")
load("human_sample_covariance.Robj")
metadata <- read.table("human_sample_metadata.txt",header=TRUE,stringsAsFactors=FALSE)
haak.dataset <- list("sample.cov" = sample.cov,
					 "geo.coords" = cbind(metadata$lon,metadata$lat),
					 "time.coords" = metadata$time)
save(haak.dataset,file="haak.dataset.Robj")