setwd("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/wolves/data")

raw.geno.data <- as.matrix(read.table("nonAdmix_nacanids_94indiv_unrel_noYNP_42Ksnps_wEcotypes.tped",stringsAsFactors=FALSE))
snp.chromos <- raw.geno.data[,1]
snp.ID <- raw.geno.data[,2]
raw.snp.data <- raw.geno.data[,5:ncol(raw.geno.data)]

ind.count.data <- (raw.snp.data[,seq(1,ncol(raw.snp.data),2)] == 1) +
		   		  (raw.snp.data[,seq(2,ncol(raw.snp.data),2)] == 1)

ind.sample.sizes <- (raw.snp.data[,seq(1,ncol(raw.snp.data),2)] != "N") +
		   		    (raw.snp.data[,seq(2,ncol(raw.snp.data),2)] != "N")

ind.freqs <- t(ind.count.data/ind.sample.sizes)
ind.cov <- cov(t(ind.freqs),use="pairwise.complete.obs")
image(ind.cov)

metadata <- read.csv("~/Desktop/wolves/EnvironmentalData_6pops_94indiv.csv",header=TRUE,stringsAsFactors=FALSE)
require(maps)
quartz(width=8,height=4)
map(database="world",
		xlim=range(metadata$X.Longitude)+c(-8,3),
		ylim=range(metadata$Y.Latitude)+c(5,5),
		mar=c(1,1,1,1))
points(metadata$X.Longitude,metadata$Y.Latitude,pch=20,col="blue")
box(lwd=2)

geoDist <- fields::rdist.earth(cbind(metadata$X.Longitude,metadata$Y.Latitude))
plot(geoDist,ind.cov,ylim=c(0,0.07))

set.seed(123)
freq.data <- switcharoo.data(ind.freqs)
obsSigma <- cov(t(freq.data),use="pairwise.complete.obs")

wolf.data <- list("geoCoords" = cbind(metadata$X.Longitude,metadata$Y.Latitude),
					"geoDist" = geoDist,
					"obsSigma" = obsSigma,
					"n.loci" = ncol(freq.data),
					"n.ind" = nrow(obsSigma),
					"sample.sizes" = colMeans(ind.sample.sizes))

save(wolf.data,file="wolf.dataset.Robj")