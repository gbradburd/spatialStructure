source("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/ms_sims/ms.spStr.sim.R")

ms.command.line <- generate.ms.command.line.values(diploid.population.size = 1e3,
														locus.size = 1e6, 
														per.bp.mu = 1e-8, 
														migration.fraction = 0.0005,
														deep.split = 1e4,
														shallow.split = 5e3)
pop1.sources <- c(1,2,6,7,11,12,16,17,21,22)
pop2.sources <- c(29,30,34,35,39,40,44,45,49,50)
pop1.targets <- pop2.sources - 25
pop2.targets <- pop1.sources + 25
admixers.s <- c(28,33,38,43,48)
admixers.t <- admixers.s - 25
admix.sources <- c(pop1.sources,pop2.sources,admixers.s)
admix.targets <- c(pop2.targets,pop1.targets,admixers.t)
admix.props <- c(rep(1,20),gtools::rdirichlet(length(admixers.s ),c(1.5,1.5))[,1])

sim.dataset <- 	generate.spStr.dataset(n.loci = 1e4,
										K = 2,
										sampling.coords = as.matrix(expand.grid(1:5,1:5,stringsAsFactors=FALSE)),
										n.chromo = 2,
										theta = ms.command.line$theta,
										migration.rate = ms.command.line$m,
										shallow.split = ms.command.line$shallow.split,
										deep.split = ms.command.line$deep.split,
										admix.list = list("sources" = admix.sources,
														  "targets" = admix.targets,
														  "admixture.proportions" = admix.props,
														  "time.point" = 0.0001),
										drop="admix.sources")
save(sim.dataset,file="sim.dataset.Robj")