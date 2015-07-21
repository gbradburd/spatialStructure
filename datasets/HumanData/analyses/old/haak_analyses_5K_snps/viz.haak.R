source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/spStr.metaplots.R")
source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/viz_funcs.R")
require(maps)

big.dir <- "~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/haak_analyses"
n.K <- 6

for(i in 1:n.K){
	dir <- paste(big.dir,"/spatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
}
#nonspatial
for(i in 1:n.K){
	dir <- paste(big.dir,"/nonspatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
}


make.all.metaplots(big.dir,big.dir,6)

