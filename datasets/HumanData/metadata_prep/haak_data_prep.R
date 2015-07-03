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
sample.cov <- total.cov
to.drop <- c(grep("other",row.names(sample.cov)),grep("Europe",row.names(sample.cov)))
sample.cov <- sample.cov[-to.drop,-to.drop]
save(sample.cov,file="human_sample_covariance.Robj")



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


#metadata
pops <- scan("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/metadata_prep/haak_pop_names.txt",what="character")
haak.metadata <- read.csv("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/metadata_prep/Haak_table.csv",header=TRUE,skip=1,stringsAsFactors=FALSE)
laz.metadata <- read.csv("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/metadata_prep/Lazaridis_table_poplatlong.csv",header=TRUE,stringsAsFactors=FALSE)
popres.metadata <- read.table("~/Desktop/Dropbox/InspectorSpaceTime/spatialStructure/datasets/HumanData/metadata_prep/popres_geo_data.txt",header=TRUE,stringsAsFactors=FALSE)
	popres.metadata$COUNTRY_SELF[2] <- "Swiss-French"
	popres.metadata$COUNTRY_SELF[4] <- "Czech-Republic"
	popres.metadata$COUNTRY_SELF[5] <- "Swiss-German"
	popres.metadata$COUNTRY_SELF[15] <- "United-Kingdom"
sample.metadata <- as.data.frame(matrix(cbind(pops,0,0,0,0),nrow=length(pops),ncol=5),stringsAsFactors=FALSE)
	colnames(sample.metadata) <- c("Population","test.pop","lon","lat","time")

matches <- match(sample.metadata$Population,laz.metadata$Population)
matches2 <- match(laz.metadata[matches[which(!is.na(matches))],]$Population,sample.metadata$Population)
sample.metadata[matches2,2] <- laz.metadata$Population[matches[which(!is.na(matches))]]
sample.metadata[matches2,3] <- laz.metadata$Long[matches[which(!is.na(matches))]]
sample.metadata[matches2,4] <- laz.metadata$Lat[matches[which(!is.na(matches))]]


samp.in.haak <- match(haak.metadata$Culture.ID,sample.metadata$Population)
	samp.in.haak <- unique(samp.in.haak[which(!is.na(samp.in.haak))])
haak.samp.variance <- vector("list",length=length(samp.in.haak))

process.sample.times <- function(element){
	# recover()
	time <- unlist(lapply(1:length(element$date),function(i){gsub(" calBCE","",element$date[i])}))
	time <- unlist(lapply(1:length(time),function(i){gsub(" cal BCE","",time[i])}))
	time <- unlist(lapply(1:length(time),function(i){gsub(" BCE","", time[i])}))
	multiples <- which(unlist(lapply(lapply(time,grep,pattern=","),sum))==1)
	if(sum(multiples) > 0){
		for(z in 1:length(multiples)){
			tmp.times <- lapply(1:length(multiples[z]),
									function(i){
										unlist(strsplit(time[multiples[i]]," , "))})
			tmp.times2 <- strsplit(tmp.times[[1]],",")
			time[multiples[z]] <- mean(unlist(lapply(1:length(tmp.times2[[1]]),function(i){mean(as.numeric(strsplit(tmp.times2[[1]][i],"-")[[1]]))})))
		}
	}
	mean.times <- unlist(lapply(1:length(time),function(i){mean(as.numeric(strsplit(time[i],"-")[[1]]))}))
	return(mean.times)
}

return.mean.spatime.locations <- function(element){
	# recover()
	loc.mat <- colMeans(cbind(as.numeric(element$Lon),as.numeric(element$Lat)))
	mean.times <- mean(process.sample.times(element))
	return(list("loc"=loc.mat,"time" = mean.times))
}

for(i in 1:length(samp.in.haak)){
	haak.samp <- match(haak.metadata$Culture.ID,sample.metadata$Population[samp.in.haak[i]])
	haak.samp <- which(!is.na(haak.samp))
	haak.samp.variance[[i]] <- list("date" = haak.metadata$Date[haak.samp],
									"Lon" = haak.metadata$Longitude[haak.samp],
									"Lat" = haak.metadata$Latitude[haak.samp])
	haak.samp.data.i <- return.mean.spatime.locations(haak.samp.variance[[i]])
	sample.metadata$test.pop[samp.in.haak[i]] <- unique(haak.metadata$Culture.ID[haak.samp])
	sample.metadata$lon[samp.in.haak[i]] <- haak.samp.data.i$loc[1]
	sample.metadata$lat[samp.in.haak[i]] <- haak.samp.data.i$loc[2]
	sample.metadata$time[samp.in.haak[i]] <- -haak.samp.data.i$time
}

popres.in.sample.matches <- match(sample.metadata$Population,popres.metadata$COUNTRY_SELF)
popres.in.sample.matches <- popres.in.sample.matches[which(!is.na(popres.in.sample.matches))]
sample.in.popres.matches <- match(popres.metadata$COUNTRY_SELF,sample.metadata$Population)
sample.in.popres.matches <- sample.in.popres.matches[which(!is.na(sample.in.popres.matches))]
sample.metadata$test.pop[sample.in.popres.matches] <- popres.metadata$COUNTRY_SELF[popres.in.sample.matches]
sample.metadata$lon[sample.in.popres.matches] <- popres.metadata$long[popres.in.sample.matches]
sample.metadata$lat[sample.in.popres.matches] <- popres.metadata$lat[popres.in.sample.matches]
to.drop <- c(grep("other",sample.metadata$Population),grep("Europe",sample.metadata$Population))
sample.metadata <- sample.metadata[-to.drop,]
write.table(sample.metadata[,c(1,3,4,5)],file="human_sample_metadata.txt")


# lapply(haak.samp.variance,check.unique)
# lapply(haak.samp.variance,return.unique)
# return.spatemp.locations(haak.samp.variance[[12]])
# match(haak.metadata$Culture.ID,sample.metadata$Population[samp.in.haak[i]])
# spatemp.data <- lapply(haak.samp.variance, return.mean.spatime.locations)
# sample.metadata

#require(maps)
# map(database="world",xlim=c(-20,70),ylim=c(35,65))
# # color.vec <- rainbow(18)
	# lapply(1:length(spatemp.data),function(i){map(database="world",xlim=c(-20,70),ylim=c(35,65)) ; points(spatemp.data[[i]],pch=i,col=color.vec[i]) ; Sys.sleep(1)})
# check.unique <- function(element){
	# unique <- TRUE
	# un.Date <- length(unique(element$date))
	# un.Lon <- length(unique(element$Lon))
	# un.Lat <- length(unique(element$Lat))
	# if(any(c(un.Date,un.Lon,un.Lat) > 1)){
		# unique <- FALSE
	# }
	# return(unique)
# }

# return.unique <- function(element,check,return.unique){
	# # recover()
	# un.Date <- unique(element$date)
	# un.Lon <- unique(element$Lon)
	# un.Lat <- unique(element$Lat)
	# unique <- list("Date" = un.Date,"Lon" = un.Lon,"Lat" = un.Lat)
	# return(unique)
# }

# return.spa.locations <- function(element){
	# # recover()
	# loc.mat <- cbind(as.numeric(element$Lon),as.numeric(element$Lat))
	# return(loc.mat)
# }



