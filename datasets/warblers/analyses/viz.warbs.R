source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/spStr.metaplots.R")
source("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/viz_funcs.R")
require(maps)

for(i in 1:10){
	dir <- paste("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers/analyses/new_runs/spatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
	pdf(file=paste("sp_admix_map_k=",i,".pdf",sep=""),width=10,height=5)
		map()
		make.admix.pie.plot(super.list,data.list,all.colors=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray"),
							cluster.names=paste(1:i),radii=2,add=TRUE,title=NULL,xlim=NULL,ylim=NULL)
	dev.off()
}
#nonspatial
for(i in 1:10){
	dir <- paste("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers/analyses/new_runs/nonspatial/k_",i,sep="")
	make.all.the.plots(dir,dir,i)
	load(list.files(pattern="output"))
	load(list.files(pattern="data.list"))
	pdf(file=paste("nsp_admix_map_k=",i,".pdf",sep=""),width=10,height=5)
		map()
		make.admix.pie.plot(super.list,data.list,all.colors=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray"),
							cluster.names=paste(1:i),radii=2,add=TRUE,title=NULL,xlim=NULL,ylim=NULL)
	dev.off()
}


make.all.metaplots("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers/analyses/new_runs","~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/warblers/analyses/new_runs",10)

