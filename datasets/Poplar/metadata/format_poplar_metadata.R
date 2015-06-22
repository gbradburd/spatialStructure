load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/Poplar_freqs.Robj")

geo.coords <- cbind(pop.metadata$Longitude,pop.metadata$Latitude)[2:435,]
save(poplar.cov,geo.coords,file="poplar_data.Robj")