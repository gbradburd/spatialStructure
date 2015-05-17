source("SpaceTimeStructureMix.R")

load("human_sample_covariance.Robj")
metadata <- read.table("human_sample_metadata.txt",header=TRUE,stringsAsFactors=FALSE)
sim.data <- list("geo.coords" = cbind(metadata$lon,metadata$lat),
				"time.coords" = metadata$time,
				"sample.covariance" = sample.cov,
				"n.loci" = 87158)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 4,
						"temporal.sampling"=FALSE,
						no.st=FALSE)
mcmc.options = list("ngen" = 5e7,
					"samplefreq" = 5e4,
					"printfreq" = 1e3,
					"savefreq" = 5e6,
					"output.file.name"="haak_k4_s_output.Robj")
MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)
