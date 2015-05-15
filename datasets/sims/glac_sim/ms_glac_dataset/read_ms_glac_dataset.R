setwd("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/glac_sim/ms_glac_dataset")
dirs <- list.dirs()
dirs <- dirs[-which(dirs==".")]

hap.snp.dataset <- matrix(0,nrow=60,ncol=length(dirs))
for(i in 1:ncol(hap.snp.dataset)){
	setwd(dirs[i])
	text <- scan("msoutput.txt",what=character(0), sep="\n", quiet=TRUE)
	if(length(text) > 0){
		marker <- grep("positions",text)
		hap.snp.dataset[,i] <- as.numeric(text[(marker+1):length(text)])
	}
	if(file.exists("msarg.txt")){
		file.remove(file="msarg.txt")
	}
	setwd("..")
}

ind.snp.dataset <- matrix(0,nrow=30,ncol=ncol(hap.snp.dataset))
for(i in 1:nrow(ind.snp.dataset)){
	ind.snp.dataset[i,] <- colSums(hap.snp.dataset[c(i*2-1,i*2),])
}

freq.mat <- ind.snp.dataset/2
load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/glac_sim/ms_glac_dataset/ms_glac_scenario.Robj")
pmat <- diag(30) - 1/30
eig.covmat <- eigen(pmat %*% cov(t(freq.mat)) %*% pmat)
plot.col <- rev(rainbow(30,start=4/6,end=6/6))[as.numeric(cut(eig.covmat$vectors[,1],30))]

geo.dist <- fields::rdist(ms_scenario$sample.config[,1:2])
plot(geo.dist,cov(t(ind.snp.dataset)))
plot(geo.dist,cov(t(freq.mat)))
plot(ms_scenario$sample.config[,c(2,1)],pch=19,col=plot.col,cex=2)

plot(eig.covmat$vectors[,1], eig.covmat$vectors[,2],pch=19,col=plot.col)

ms_glac_dataset <- list("sample.cov" = cov(t(freq.mat)),
						"geo.coords" = ms_scenario$sample.config[,c(2,1)],
						"n.loci" = ncol(ind.snp.dataset))
save(ms_glac_dataset,file="ms_glac_dataset.Robj")