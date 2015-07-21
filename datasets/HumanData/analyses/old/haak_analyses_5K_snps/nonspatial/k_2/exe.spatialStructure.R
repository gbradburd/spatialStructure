source("SpaceTimeStructureMix.R")
load("haak.dataset.Robj")
sim.data <- list("geo.coords" = haak.dataset$geo.coords,
				"time.coords" = haak.dataset$time.coords,
				"sample.covariance" = haak.dataset$sample.cov,
				"n.loci" = 5e3)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 2,
						"temporal.sampling" = FALSE,
						"no.st" = TRUE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 5e6,
					"output.file.name" = "k2_haak_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
