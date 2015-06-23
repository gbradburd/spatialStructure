source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/spStr.metaplots.R")
source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/viz_funcs.R")
require(maps)

big.dir <- "~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/poplar/poplar_analyses"

for(i in 1:10){
	dir <- paste(big.dir,"/spatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
}
#nonspatial
for(i in 1:10){
	dir <- paste(big.dir,"/nonspatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
}


make.all.metaplots(big.dir,big.dir,10)

