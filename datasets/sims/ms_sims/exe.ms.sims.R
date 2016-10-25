source("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/ms_sims/ms.spStr.sim.R")

ms.command.line <- generate.ms.command.line.values(diploid.population.size = 1e3,
														locus.size = 1e6, 
														per.bp.mu = 1e-8, 
														migration.fraction = 0.0005,
														deep.split = 1e4,
														shallow.split = 5e3)
														
sim.dataset <- 	generate.spStr.dataset(n.loci = 1e3,
										K = 2,
										sampling.coords = as.matrix(expand.grid(1:4,1:4,stringsAsFactors=FALSE)),
										n.chromo = 20,
										theta = ms.command.line$theta,
										migration.rate = ms.command.line$m,
										shallow.split = ms.command.line$shallow.split,
										deep.split = ms.command.line$deep.split,
										admix.list = list("sources" = c(1,2),
														  "targets" = c(17,18),
														  "admixture.proportions" = c(0.3,0.9),
														  "time.point" = 0.0001),
										drop="admix.sources")


save(sim.dataset,file="sim.dataset.Robj")