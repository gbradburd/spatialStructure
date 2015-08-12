source("SpaceTimeStructureMix.R")
load("haak.dataset.Robj")
# load(list.files(pattern="init"))
sim.data <- list("geo.coords" = haak.dataset$geo.coords,
				"time.coords" = haak.dataset$time.coords,
				"sample.covariance" = haak.dataset$sample.cov,
				"n.loci" = 87158)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 6,
						"temporal.sampling" = FALSE,
						"no.st" = TRUE)
mcmc.options = list("ngen" = 3e7,
					"samplefreq" = 3e4,
					"printfreq" = 1e3,
					"savefreq" = 5e6,
					"output.file.name" = "k6_haak_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options)
#initial.parameters=initial.parameters)
