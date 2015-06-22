source("SpaceTimeStructureMix.R")

load(list.files(pattern="dataset.Robj"))

sample.cov <- poplar.data$sample.cov
geo.coords <- poplar.data$geo.coords

sim.data <- list("geo.coords" = geo.coords,
				"time.coords" = matrix(0,nrow=nrow(geo.coords),ncol=2),
				"sample.covariance" = sample.cov,
				"n.loci" = 1000)
model.options = list("round.earth" = TRUE,
						"n.clusters" = 4,
						"temporal.sampling"=FALSE,
						"no.st" = TRUE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name" = "k4_poplar_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
