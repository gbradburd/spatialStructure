source("SpaceTimeStructureMix.R")

load(list.files(pattern="dataset.Robj"))

sample.cov <- globe.data$globe.sample.covariance
globe.coords <- globe.data$globe.coords

sim.data <- list("geo.coords" = globe.coords,
				"time.coords" = matrix(0,nrow=nrow(globe.coords),ncol=2),
				"sample.covariance" = sample.cov,
				"n.loci" = 10000)
model.options = list("round.earth" = TRUE,
						"n.clusters" = 2,
						"temporal.sampling"=FALSE,
						"no.st" = TRUE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name" = "k2_globe_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
