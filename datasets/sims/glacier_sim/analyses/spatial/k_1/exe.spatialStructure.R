source("SpaceTimeStructureMix.R")

load(list.files(pattern="dataset.Robj"))

sim.data <- list("geo.coords" = ms_glac_dataset$geo.coords,
				"time.coords" = matrix(0,nrow=nrow(ms_glac_dataset$geo.coords),ncol=2),
				"sample.covariance" = ms_glac_dataset$sample.cov,
				"n.loci" = ms_glac_dataset$n.loci)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 1,
						"temporal.sampling"=FALSE,
						no.st=FALSE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name"="k_1_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
