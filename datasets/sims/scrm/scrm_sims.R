#SpatialStructure Simulations
setwd("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/scrm")
require(scrm)
require(ape)

total.grid <- as.matrix(expand.grid(1:100,1:100,stringsAsFactors=FALSE))
random.pts <- total.grid[sample(1:10000,20),]
plot(total.grid)
	points(random.pts,col="red",pch=20)



sim.coords <- cbind(runif(10),runif(10))
sample.dist <- fields::rdist(sim.coords)
plot(sample.dist,dexp(sample.dist,rate=migration.rate))

vignette("scrm-Arguments")

# first number is number of haplotypes
# second number is number of loci
# -r is followed by the recombination rate (4*Ne*r)
# 	and the total length of all loci
# -t is followed by the mutation rate (4*Ne*mu)
test <- scrm("5 1 -r 5 1000 -t 5 -L ")
plot(read.tree(text=test$trees[[1]]))
?scrm

