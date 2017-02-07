require(geoStructure)
load(list.files(pattern="dataset.Robj"))
						
geoStructure(spatial = TRUE,
					  K = as.numeric(tail(strsplit(getwd(),"=")[[1]],1)),
					  freqs = sim.dataset$data.list$allele.freqs,
					  D = fields::rdist(sim.dataset$data.list$coords),
					  coords = sim.dataset$data.list$coords,
					  sample.sizes = rowMeans(sim.dataset$data.list$sample.sizes),
					  prefix = tail(strsplit(getwd(),"/")[[1]],1),
					  n.iter = 1e4)