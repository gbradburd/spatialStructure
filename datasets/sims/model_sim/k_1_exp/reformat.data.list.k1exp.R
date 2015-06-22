load("/Users/gburd/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/model_sim/k_1_exp/k1_exp_dataset.Robj")
geo.coords <- sim.data$sampling.data$geo.coords
time.coords <- sim.data$sampling.data$temporal.coords


sp.path <- "~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/model_sim/k_1_exp/spatial/"
for(i in 1:4){
	setwd(paste(sp.path,"k_",i,sep=""))
	load(list.files(pattern="data.list"))
	data.list$geo.coords <- geo.coords
	data.list$time.coords <- time.coords
	save(data.list,file="data.list.Robj")
	setwd("..")
}

nsp.path <- "~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/model_sim/k_1_exp/nonspatial/"
for(i in 1:4){
	setwd(paste(nsp.path,"k_",i,sep=""))
	load(list.files(pattern="data.list"))
	data.list$geo.coords <- geo.coords
	data.list$time.coords <- time.coords
	save(data.list,file="data.list.Robj")
	setwd("..")
}
