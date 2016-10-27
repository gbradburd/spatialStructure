require(geoStructure)
load(list.files(pattern="dataset.Robj"))
						
geoStructure(spatial = TRUE,
			 K = 2,
			 freqs = sim.dataset$data.list$allele.freqs,
			 D = fields::rdist(sim.dataset$data.list$coords),
			 coords = fields::rdist(sim.dataset$data.list$coords),
			 sample.sizes = rowMeans(sim.dataset$data.list$sample.sizes),
			 prefix = tail(strsplit(getwd(),"/")[[1]],1),
			 n.iter = 5e3)
			 
load(list.files(pattern="geoStr.results"))
load(list.files(pattern="data.block"))

cluster.order <- get.cluster.order(data.block$K,geoStr.results[[1]]$post$admix.proportions,ref.admix.props=sim.dataset$par.list$admix.props)

plot.admix.prop.post <- function(K,N,admix.props,sim.admix.props,K.cols,cluster.order){
#	recover()
	qs <- lapply(1:K,function(i){
				apply(admix.props[,,i],2,
					function(x){quantile(x,c(0.025,0.975))})
			})
	mean.props <- lapply(1:K,function(i){
						apply(admix.props[,,i],2,
							function(x){mean(x)})
					})
	cluster.cols <- K.cols[cluster.order[,2]]
	par(mfrow=c(1,2))
	plot(0,ylim=c(0,1),type='n',xlim=c(0,N),ylab="admix proportion",xlab="Sample")
		lapply(1:K,function(i){segments(x0=1:N,x1=1:N,
									    y0=apply(qs[[i]],2,min),
									    y1=apply(qs[[i]],2,max),
										col=cluster.cols[i],lwd=2)})
		lapply(1:K,function(i){points(mean.props[[i]],col=cluster.cols[i])})
		lapply(1:K,function(i){points(sim.admix.props[,i],col=K.cols[i],pch=18)})
	plot(0,xlim=c(0,1),ylim=c(0,1))
		abline(0,1,col=1)
		lapply(1:K,function(k){
				lapply(seq(1,dim(geoStr.results[[1]]$post$admix.proportions)[1],length.out=100),function(i){
						points(sim.admix.props[,cluster.order[k,2]],
								admix.props[i,,k],
								col=adjustcolor(K.cols[cluster.order[k,2]],0.1),
								pch=20)
					})
			})
	return(invisible("plotted"))
}
	

pdf(file="admix.prop.fit.pdf",width=9,height=5)
	plot.admix.prop.post(data.block$K,data.block$N,geoStr.results[[1]]$post$admix.proportions,sim.dataset$par.list$admix.props,c("blue","red"),cluster.order=cluster.order)
dev.off()

