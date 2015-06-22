source("SpaceTimeStructureMix.R")

load(list.files(pattern="dataset.Robj"))

sample.cov <- warb.data$warb.sample.covariance
warbler.ind.coords <- warb.data$warbler.ind.coords

sim.data <- list("geo.coords" = warbler.ind.coords,
				"time.coords" = matrix(0,nrow=nrow(warbler.ind.coords),ncol=2),
				"sample.covariance" = sample.cov,
				"n.loci" = 2247)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 1,
						"temporal.sampling"=FALSE,
						"no.st" = TRUE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name" = "k1_warb_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
