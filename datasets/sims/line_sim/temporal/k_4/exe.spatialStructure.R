source("SpaceTimeStructureMix.R")

load(list.files(pattern="dataset.Robj"))
sim.sample.covariance <- cov(t(spacemix.dataset$allele.counts/spacemix.dataset$sample.sizes))
sim.data <- list("geo.coords" = matrix(0,nrow=nrow(spacemix.dataset$population.coordinates),ncol=2),
				"time.coords" = spacemix.dataset$population.coordinates[,1],
				"sample.covariance" = sim.sample.covariance,
				"n.loci" = ncol(spacemix.dataset$allele.counts))
model.options = list("round.earth" = FALSE,
						"n.clusters" = 4,
						"temporal.sampling"=TRUE,
						no.st=FALSE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name"="k_4t_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
