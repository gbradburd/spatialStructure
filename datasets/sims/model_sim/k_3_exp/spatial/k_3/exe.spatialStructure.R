source("SpaceTimeStructureMix.R")

load(list.files(pattern="dataset"))

# freqs <- switcharoo.data(spacemix.dataset$allele.counts/spacemix.dataset$sample.sizes)

# sim.sample.covariance <- cov(t(freqs))
sim.data <- list("geo.coords" = sim.data$sampling.data$geo.coords,
				"time.coords" = sim.data$sampling.data$temporal.coords,
				"sample.covariance" = sim.data$sample.covariance,
				"n.loci" = sim.data$n.loci)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 3,
						"temporal.sampling"=TRUE,
						no.st=FALSE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name"="k_3exp_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)