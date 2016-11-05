require(geoStructure)
load(list.files(pattern="dataset"))
geoStructure(spatial = TRUE,
			 K = 2,
			 freqs = bear.dataset$sample.freqs,
			 D = fields::rdist.earth(bear.dataset$sample.coords),
			 coords = bear.dataset$sample.coords,
			 sample.sizes = bear.dataset$sample.sizes,
			 prefix = "bears_K=2")