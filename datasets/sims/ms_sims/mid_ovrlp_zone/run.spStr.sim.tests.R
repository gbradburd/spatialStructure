deep.split.times <- seq(5e3,5e2,length.out=10)
for(i in 1:length(deep.split.times)){
	dir.name <- paste0("sim_", deep.split.times[i])
	setwd(dir.name)
		for(j in 1:10){
			dir.name <- paste0("sim_", j)
			setwd(dir.name)
				source("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/ms_sims/exe.spStr.R")
			setwd("..")
		}
	setwd("..")
}