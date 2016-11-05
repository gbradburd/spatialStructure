require(BEDASSLE)
require(geoStructure)
setwd("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/sims/ms_sims/mid_ovrlp_zone")
deep.split.times <- seq(5e3,5e2,length.out=10)


get.lr.fst <- function(sim.dataset){
	c1 <- sim.dataset$data.list$allele.freqs[sim.dataset$data.list$coords[,1] < 3,] * 
								  sim.dataset$data.list$sample.sizes[sim.dataset$data.list$coords[,1] < 3,]
	s1 <- sim.dataset$data.list$sample.sizes[sim.dataset$data.list$coords[,1] < 3,]
	c2 <- sim.dataset$data.list$allele.freqs[sim.dataset$data.list$coords[,1] > 3,] * 
								  sim.dataset$data.list$sample.sizes[sim.dataset$data.list$coords[,1] > 3,]
	s2 <- sim.dataset$data.list$sample.sizes[sim.dataset$data.list$coords[,1] > 3,]
	fst <- calculate.pairwise.Fst(rbind(colSums(c1),colSums(c2)),
									rbind(colSums(s1),colSums(s2)))
	return(fst)
}

get.ad.prop.fit <- function(sim.dataset,geoStr.results){
	cluster.order <- get.cluster.order(K = 2,
										admix.props = geoStr.results[[1]]$post$admix.proportions,
										ref.admix.props = sim.dataset$par.list$admix.props)
	sum.sq.dif <- sum(sqrt((geoStr.results[[1]]$point$admix.proportions[,cluster.order[,2]][,1] -
					 	sim.dataset$par.list$admix.props[,1])^2))
	return(sum.sq.dif)
}

x.pop.fst <- vector("list",length=length(deep.split.times))
for(i in 1:length(deep.split.times)){
	setwd(paste0("sim_",deep.split.times[i]))
	x.pop.fst[[i]] <- numeric(10)
		for(j in 1:10){
			setwd(paste0("sim_",j))
				load("sim.dataset.Robj")
					x.pop.fst[[i]][j] <- get.lr.fst(sim.dataset)
			setwd("..")
		}
	setwd("..")	
}

ad.prop.fit <- vector("list",length=length(deep.split.times))
for(i in 1:length(deep.split.times)){
	setwd(paste0("sim_",deep.split.times[i]))
	ad.prop.fit[[i]] <- numeric(10)
		for(j in 1:10){
			setwd(paste0("sim_",j))
				load("sim.dataset.Robj")
				load(list.files(pattern="geoStr.results"))
					ad.prop.fit[[i]][j] <- get.ad.prop.fit(sim.dataset,geoStr.results)
			setwd("..")
		}
	setwd("..")	
}

pdf(file="sim.meta.figs.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(5,5.5,4,1))
plot(deep.split.times[unlist(lapply(1:10,function(x){rep(x,10)}))],
		unlist(x.pop.fst),pch=19,col=adjustcolor(1,0.5),
		xlab = "deep population split time",
		ylab = "pairwise Fst between \nclusters 1 and 2")
plot(deep.split.times[unlist(lapply(1:10,function(x){rep(x,10)}))],
		unlist(ad.prop.fit)/25,pch=19,col=adjustcolor(1,0.5),
		xlab = "deep population split time",
		ylab = "avg. sum sq. difference between \ntrue and esimated admix proportions")
dev.off()