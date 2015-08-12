source("SpaceTimeStructureMix.R")
load("haak.dataset.Robj")
sim.data <- list("geo.coords" = haak.dataset$geo.coords,
				"time.coords" = haak.dataset$time.coords,
				"sample.covariance" = haak.dataset$sample.cov,
				"n.loci" = 5e3)
model.options = list("round.earth" = FALSE,
						"n.clusters" = 3,
						"temporal.sampling" = TRUE,
						"no.st" = FALSE)
mcmc.options = list("ngen" = 1e6,
					"samplefreq" = 1e3,
					"printfreq" = 1e3,
					"savefreq" = 1e5,
					"output.file.name" = "k3_haak_output.Robj")
Rprof(file="spStr_Rprof.out")

MCMC.gid(sim.data,model.options,mcmc.options,initial.parameters=NULL)


Rprof(NULL)

require(proftools)
pdf(file="~/Desktop/spStr_Rprof/spStr_rprof_callgraph.pdf",width=10,height=7)
	plotProfileCallGraph(readProfileData("~/Desktop/spStr_Rprof/spStr_Rprof.out"),score="self")
dev.off()

pdf(file="~/Desktop/spStr_Rprof/spStr_rprof_callgraph_total.pdf",width=10,height=7)
	plotProfileCallGraph(readProfileData("~/Desktop/spStr_Rprof/spStr_Rprof.out"),score="total")
dev.off()