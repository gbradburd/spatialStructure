load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers/warbler.ind.dataset.Robj")
geo.coords <- warb.data$warbler.ind.coords
time.coords <- rep(0,nrow(warb.data$warbler.ind.coords))


sp.path <- "~/Desktop/new_runs/spatial/"
for(i in 1:8){
	setwd(paste(sp.path,"k_",i,sep=""))
	load(list.files(pattern="data.list"))
	data.list$geo.coords <- geo.coords
	data.list$time.coords <- time.coords
	save(data.list,file="data.list.Robj")
	setwd("..")
}

nsp.path <- "~/Desktop/new_runs/nonspatial/"
for(i in 1:8){
	setwd(paste(nsp.path,"k_",i,sep=""))
	load(list.files(pattern="data.list"))
	data.list$geo.coords <- geo.coords
	data.list$time.coords <- time.coords
	save(data.list,file="data.list.Robj")
	setwd("..")
}

prob.vec <- numeric(8)
sp.path <- "~/Desktop/new_runs/spatial/"
for(i in 1:8){
	setwd(paste(sp.path,"k_",i,sep=""))
	load(list.files(pattern="output"))
	prob.vec[i] <- mean(super.list$output.list$posterior.prob)
	setwd("..")
}