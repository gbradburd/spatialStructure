source("SpaceTimeStructureMix.R")

load(list.files(pattern="dataset.Robj"))

sample.cov <- trichocarpa.data$sample.cov
geo.coords <- trichocarpa.data$geo.coords

data <- list("geo.coords" = geo.coords,
				"time.coords" = matrix(0,nrow=nrow(geo.coords),ncol=2),
				"sample.covariance" = sample.cov,
				"n.loci" = 1000)
model.options = list("round.earth" = TRUE,
						"n.clusters" = 9,
						"temporal.sampling" = FALSE,
						"no.st" = TRUE)
mcmc.options = list("ngen" = 1e7,
					"samplefreq" = 1e4,
					"printfreq" = 1e3,
					"savefreq" = 1e6,
					"output.file.name" = "k9_tri_output.Robj")
MCMC.gid(data,model.options,mcmc.options,initial.parameters=NULL)
