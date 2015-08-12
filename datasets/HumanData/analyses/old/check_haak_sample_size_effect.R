metadata <- read.table("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/metadata_prep/human_sample_metadata.txt",
						header=TRUE,stringsAsFactors=FALSE)
						
sample.sizes <- metadata$sample.size
time.cols <- rainbow(length(unique(metadata$time)),start=4/6,end=6/6)[as.numeric(cut(metadata$time,length(unique(metadata$time))))]

load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/old/haak_analyses/spatial/k_3/k3_haak_output.Robj")

plot(1/sample.sizes,super.list$parameter.list$nuggets,xlab="1/sample size",ylab="nugget",type='n')
	text(1/sample.sizes,super.list$parameter.list$nuggets,labels=metadata$Population,col=time.cols)


plot(sample.sizes,super.list$parameter.list$admix.proportions[,2],
		xlab="sample size",ylab="admix.prop.Clst2",type='n',xlim=c(-200,1000))
	text(sample.sizes,super.list$parameter.list$admix.proportions[,2],
		labels=metadata$Population,col=time.cols)
		
load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/analyses/old/haak_analyses/spatial/k_3/data.list.Robj")
load("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/metadata_prep/haak_snp_chr_1.Robj")
vars <- unlist(lapply(1:nrow(chr.data.snps),function(i){var(chr.data.snps[i,])}))
to.drop <- c(grep("other",row.names(chr.data.snps)),grep("Europe",row.names(chr.data.snps)))
vars <- vars[-to.drop]
epsilon <- colMeans(chr.data.snps[-to.drop,])
globe.mean.var <- mean(epsilon*(1-epsilon))

plot(min(data.list$sample.covariance)*(1+1/sample.sizes),vars,xlab="1/sample.size",ylab="variance in allele freqs",col=time.cols) ; abline(0,1,col="red")

plot(diag(data.list$sample.covariance),vars) ; abline(0,1,col="red")

