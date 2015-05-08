big.table <- read.table("~/Desktop/Haak_POPRES_freqs.frq.strat.gz",header=TRUE,stringsAsFactors=FALSE)
n.lines.total <- nrow(big.table)
chr.breaks <- matrix(0,nrow=22,ncol=2)
chr.breaks[,1] <- unlist(lapply(1:max(unique(big.table$CHR)),match,table=big.table$CHR))
chr.breaks[,2] <- cumsum(c(chr.breaks[2:22,1],n.lines.total+1) - chr.breaks[,1])
write.table(chr.breaks,file="haak_chr_breaks.txt",row.names=FALSE,col.names=FALSE)
rm(big.table) ; gc()

chr.breaks <- read.table("haak_chr_breaks.txt",header=FALSE,stringsAsFactors=FALSE)

random.switcharoo <- function(x){
	swap <- ifelse(runif(1) < 0.5,
					1,
					0)
	if(swap){
		x <- 1-x
	}
	return(x)
}

get.snp.data <- function(chr.data,snp.name){
	random.switcharoo(chr.data$MAF[which(!is.na(match(chr.data$SNP,snp.names[j])))])
}

for(i in 1:nrow(chr.breaks)){
	chr.data <- read.table("~/Desktop/Haak_POPRES_freqs.frq.strat.gz",
							col.names=c("CHR","SNP","CLST","A1","A2","MAF","MAC","NCHROBS"),
							nrows=chr.breaks[i,2]-chr.breaks[i,1]+1,skip=chr.breaks[i,1],stringsAsFactors=FALSE)
	pops <- unique(chr.data$CLST)
	pops[is.na(pops)] <- "other2"
	snp.names <- unique(chr.data$SNP)
	chr.data.snps <- matrix(0,nrow=length(pops),ncol=length(snp.names))
	row.names(chr.data.snps) <- pops
		for(j in 1:length(snp.names)){
			chr.data.snps[,j] <- get.snp.data(chr.data,snp.name[j])
			if(j %% 100){
				cat(j,"\t")
			}
		}
	save(chr.data.snps,file=paste("haak_snp_chr_",i,".Robj",sep=""))
}
# for(i in 1:nrow(chr.breaks)){
	# chr.data <- read.table("~/Desktop/Haak_POPRES_freqs.frq.strat.gz",
							# col.names=c("CHR","SNP","CLST","A1","A2","MAF","MAC","NCHROBS"),
							# nrows=chr.breaks[i,2]-chr.breaks[i,1]+1,skip=chr.breaks[i,1],stringsAsFactors=FALSE)
	# cat(unique(chr.data$CHR),"\n")
# }
file.list <- list.files(pattern="haak_snp_chr")
chr.cov.list <- list(length(file.list))
for(i in 1:length(file.list)){
	load(file.list[i])
	chr.cov.list[[i]] <- list("cov" = cov(t(chr.data.snps)), 
								"n.loci" = ncol(chr.data.snps))
	cat(i,"\t")
}

total.loci <- sum(unlist(lapply(chr.cov.list,"[[","n.loci")))
loci.proportions <- lapply(1:length(chr.cov.list),function(i){chr.cov.list[[i]]$n.loci/total.loci})
weighted.cov.list <- lapply(1:length(chr.cov.list),function(i){loci.proportions[[i]] * chr.cov.list[[i]]$cov})
total.cov <- Reduce("+",weighted.cov.list)


pops <- scan("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/haak_pop_names.txt",what="character")
mc.mat <- diag(nrow(total.cov)) - 1/nrow(total.cov)
eig.covmat <- eigen(mc.mat %*% total.cov %*% t(mc.mat))
par(mfrow=c(1,2))
plot(eig.covmat$vectors[,1],eig.covmat$vectors[,2],type='n')
	text(eig.covmat$vectors[,1],eig.covmat$vectors[,2],labels=pops,cex=0.5)
plot(eig.covmat$vectors[,1],eig.covmat$vectors[,2],type='n',xlim=c(-0.02,0.05),ylim=c(-0.1,0.1))
	text(eig.covmat$vectors[,1],eig.covmat$vectors[,2],labels=pops,cex=0.5)
plot(eig.covmat$vectors[,1],eig.covmat$vectors[,2],type='n',xlim=c(0.01,0.05),ylim=c(-0.1,0.05))
	text(eig.covmat$vectors[,1],eig.covmat$vectors[,2],labels=pops,cex=0.5)





plot(eig.covmat$vectors[,2],eig.covmat$vectors[,3],type='n',xlim=c(-0.07,0.05),ylim=c(-0.02,0.0001))
	text(eig.covmat$vectors[,2],eig.covmat$vectors[,3],labels=pops,cex=0.5)
