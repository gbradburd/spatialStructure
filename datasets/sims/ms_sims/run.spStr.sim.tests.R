n.sims <- 10

for(i in 1:n.sims){
	dir.name <- paste0("sim_",i)
	setwd(dir.name)
		source("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/ms_sims/exe.spStr.R")
	setwd("..")
}