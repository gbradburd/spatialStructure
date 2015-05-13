source("SpaceTimeStructureMix.R")

load("poplar_data.Robj")

sim.data <- list("geo.coords" = geo.coords,
				"time.coords" = rep(0,nrow(geo.coords)),
				"sample.covariance" = poplar.cov,
				"n.loci" = 33070)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 2,
						"temporal.sampling"=FALSE,
						no.st=FALSE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e5,
					"output.file.name"="poplar_k2_s_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
