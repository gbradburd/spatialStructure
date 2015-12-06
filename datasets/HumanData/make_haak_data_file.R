setwd("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/metadata_prep")
load("human_sample_covariance.Robj")
load("human_sample_bin_var.Robj")
metadata <- read.table("human_sample_metadata.txt",header=TRUE,stringsAsFactors=FALSE)
haak.dataset <- list("sample.cov" = sample.cov,
					 "geo.coords" = cbind(metadata$lon,metadata$lat),
					 "time.coords" = metadata$time,
					 "sample.sizes" = metadata$sample.size,
					 "bin.var" = bin.var)
save(haak.dataset,file="haak.dataset.Robj")