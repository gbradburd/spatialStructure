get.tabs <- function(x,pattern){
	n.tabs <- length(which(strsplit(x[grep(pattern,x)],"")[[1]] == "\t"))
	tabs <- paste(rep("\t",n.tabs),collapse="")
	return(tabs)
}

make.spstr.run.dirs <- function(spacemix.code,exe.code,dataset,K,path,output.file.name){
	# recover()
	if(file.exists(path)){
		stop("directory already exists\n")
	} else {
		dir.create(path)
	}
	setwd(path)
	dir.create("spatial")
	dir.create("nonspatial")
	for(k in 1:K){
		dir.name <- paste("k_",k,sep="")
		dir.create(paste("spatial/",dir.name,sep=""))
		dir.create(paste("nonspatial/",dir.name,sep=""))
		file.copy(from=spacemix.code,to=paste("spatial/",dir.name,sep=""))
		file.copy(from=spacemix.code,to=paste("nonspatial/",dir.name,sep=""))
		file.copy(from=dataset,to=paste("spatial/",dir.name,sep=""))
		file.copy(from=dataset,to=paste("nonspatial/",dir.name,sep=""))
		x <- readLines(exe.code,warn=FALSE)
		y <- x
		y[grep("n.clusters",x)] <- paste(get.tabs(x,"n.clusters"),'"n.clusters" = ',k,",",sep="")
		y[grep("output",x)] <- paste(get.tabs(x,"output"),'"output.file.name" = \"k',k,"_",output.file.name,".Robj\")",sep="")
		y.spatial <- y
		y.nonspatial <- y
		y.spatial[grep("no.st",y.spatial)] <- paste(get.tabs(y.spatial,"no.st"),'"no.st" = FALSE)',sep="")
		y.nonspatial[grep("no.st",y.spatial)] <- paste(get.tabs(y.nonspatial,"no.st"),'"no.st" = TRUE)',sep="")
		cat(y.spatial,file=paste("spatial/",dir.name,"/exe.spatialStructure.R",sep=""),sep="\n")
		cat(y.nonspatial,file=paste("nonspatial/",dir.name,"/exe.spatialStructure.R",sep=""),sep="\n")
	}
	return(invisible("rock&roll"))
}

if(FALSE){
make.spstr.run.dirs(spacemix.code="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/code/SpaceTimeStructureMix.R",
					exe.code="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/globe/exe.spatialStructure.R",
					dataset="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/globe/globe.spStr.dataset.Robj",
					K=10,
					path="~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/globe/analyses",
					output.file.name = "globe_output")
}
