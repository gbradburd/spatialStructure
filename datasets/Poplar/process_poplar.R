
##downloaded from http://datadryad.org/resource/doi:10.5061/dryad.7s848
poplar<-read.table("~/Downloads/FOR_DRYAD_TandB_fwd_june2013_434.txt",as.is=TRUE,sep="\t")
##Supp tables from the paper
pop.metadata<-read.csv(file="~/Dropbox/Students/gideon/spatialStructure/datasets/Poplar/Populus_metadata.csv",as.is=TRUE,skip=1)
pop.samples<-poplar[10,]

poplar.genos<-poplar[11:nrow(poplar),2:ncol(poplar)]
poplar.freqs<-apply(poplar.genos[1:1000,],1,function(geno){
	genotype<-sapply(as.character(geno),function(x){strsplit(x,"/|")[[1]][1:2]})
	bases<-unique(c(genotype[1,],genotype[2,]))
	my.bases<-bases[bases %in% c("A","T","C","G")]
	my.bases<-sample(my.bases)
	allele.1<-(genotype[1,] == my.bases[1]) + (genotype[2,] == my.bases[1]) 
	allele.2<-(genotype[1,] == my.bases[2]) + (genotype[2,] == my.bases[2]) 
	allele.1 / (allele.1+allele.2)
})

poplar.cov<-cov(t(poplar.freqs),use="pairwise.complete.obs")

map<-match(pop.samples,pop.metadata$Accession)

pop.metadata<-pop.metadata[map,]

save(file="~/Dropbox/Students/gideon/spatialStructure/datasets/Poplar/Poplar_freqs.Robj",poplar.freqs,poplar.cov,pop.metadata)


species.points<-c(19,22)
names(species.points)<-c("Populus trichocarpa","Populus balsamifera")

drainage.cols<-rainbow(length(unique(pop.metadata$Drainage.Location.name)))
drainage.cols<-sample(drainage.cols)
names(drainage.cols)<-unique(pop.metadata$Drainage.Location.name)

PCAs.pop<-eigen(poplar.cov)
plot(PCAs.pop$vectors[,2],PCAs.pop$vectors[,3],pch=species.points[pop.metadata$Species],col=drainage.cols[pop.metadata$Drainage.Location.name])
