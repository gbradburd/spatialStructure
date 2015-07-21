source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/make.spstr.run.dirs.R")
make.spstr.run.dirs(spacemix.code="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/SpaceTimeStructureMix.R",
					exe.code="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/exe.spatialStructure.R",
					dataset="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/haak.dataset.Robj",
					K=6,
					path="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/haak_analyses",
					output.file.name = "haak_output")


source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/SpaceTimeStructureMix.R")
sp.dir <- "~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/old/haak_analyses/spatial"
setwd(sp.dir)
for(i in 1:6){
	make.starting.params.list(dir=paste("k_",i,sep=""),
							init.par.file.name=paste("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/haak_analyses/spatial/k_",i,"/k_",i,"_init.params.list.Robj",sep=""))
	setwd(sp.dir)
}


nsp.dir <- "~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/old/haak_analyses/nonspatial"
setwd(nsp.dir)
for(i in 1:6){
	make.starting.params.list(dir=paste("k_",i,sep=""),
							init.par.file.name=paste("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/haak_analyses/nonspatial/k_",i,"/k_",i,"_init.params.list.Robj",sep=""))
	setwd(nsp.dir)
}


