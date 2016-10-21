##################################################################
##################################################################
##	run ms on a lattice
##################################################################
##################################################################

require(geoStructure)

within.cluster.m <- function(sample.inds,migration.index,pairwise.migration.matrix){
	#recover()
	migration.rate.vector <- c()
	for(i in 1:nrow(migration.index)){
		migration.rate.vector <- c(migration.rate.vector,
									sprintf("-m %s %s %s",
										sample.inds[migration.index[i,1]],
										sample.inds[migration.index[i,2]],
										pairwise.migration.matrix[migration.index[i,1],migration.index[i,2]]))
	}
	return(migration.rate.vector)
}

within.cluster.merge <- function(sample.inds,shallow.split){
	merge.vector <- unlist(lapply(sample.inds[2]:sample.inds[length(sample.inds)],
						function(i){
							sprintf("-ej %s %s %s",
									shallow.split,
									i,
									sample.inds[1])}))
	return(merge.vector)
}

btwn.cluster.merge <- function(K,n.pops,deep.split){
	if(K > 1){
		merge.vector <- unlist(lapply(2:K,
						function(k){
							sprintf("-ej %s %s %s",
								deep.split,
								1+(k-1)*n.pops,
								1)
						}))
	} else {
		merge.vector <- c()
	}
	return(merge.vector)
}

write.migration.rates <- function(K,sampling.coords,migration.rate,shallow.split,deep.split){
	# recover()
	pop.dist <- fields::rdist(sampling.coords)
	pop.dist[which(pop.dist > sqrt(2) | pop.dist < 1)] <- Inf
	pairwise.migration.matrix <- lapply(1:K,function(k){migration.rate[k]/pop.dist})
	migration.index <- lapply(pairwise.migration.matrix,function(x){which(x != 0,arr.ind=TRUE)})
	n.pops <- nrow(sampling.coords)
	sample.inds <- lapply(1:K,function(k){1:n.pops + (k-1) * n.pops})
	migration.rate.vector <- c()
	for(k in 1:K){
		migration.rate.vector <- c(migration.rate.vector,
									within.cluster.m(sample.inds[[k]],migration.index[[k]],pairwise.migration.matrix[[k]]))
	}
	for(k in 1:K){
		migration.rate.vector <- c(migration.rate.vector,
									within.cluster.merge(sample.inds[[k]],shallow.split[k]))
	}
	migration.rate.vector <- c(migration.rate.vector,
									btwn.cluster.merge(K,n.pops,deep.split))
	return(migration.rate.vector)
}

# code cannibalized from Dan Denison that reads ms output into R
read.ms.haplotype.matrices <- function(nsam, ndraws, ms.output.file) {
    txt <- scan(file=ms.output.file, what=character(0), sep="\n", quiet=TRUE)
    h <- list()
    ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
    marker <- grep("segsites", txt)
    stopifnot(length(marker) == ndraws)
    ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
    segsites <- sapply(strsplit(txt[marker], split=":"), function(vec) as.integer(vec[2]))
    for(draw in seq(along=marker)) {
        if(!(draw %% 100)) cat(draw, " ")
        if(segsites[draw] > 0) {
            haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
            haplotypes <- strsplit(haplotypes, split="")
            h[[draw]] <- sapply(haplotypes, function(el) c(as.integer(el)))
            ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
            if(segsites[draw] == 1) h[[draw]] <- as.matrix(h[[draw]])
            ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
            else h[[draw]] <- t(h[[draw]])
        }
        else h[[draw]] <- matrix(nrow=nsam, ncol=0)
        stopifnot(all(dim(h[[draw]]) == c(nsam, segsites[draw])))  
    }
    cat("\n")
    h
}

# function that actually calls ms with the call put together by all these R functions
ms <- function(sampling.coords,K,n.chromo,theta,migration.rate,shallow.split,deep.split){
		ms.output.file <- "ms_output"
		random.seeds <- c(sample(1:100000,3,replace=TRUE))
		sampled.pops <- nrow(sampling.coords)
		call <- paste(
					sprintf(
						"/Applications/ms.folder/msdir/ms %s 1 -t %s -s 1 -I %s %s 0.0-m %s -seeds %s %s %s",
							K*sampled.pops*n.chromo,
							theta,
							K*sampled.pops,
							paste0(rep(1*n.chromo,K*sampled.pops),collapse=" "),
							paste0(write.migration.rates(K,
														 sampling.coords,
														 migration.rate,
														 shallow.split,
														 deep.split),
														 collapse=" "),
							random.seeds[1],
							random.seeds[2],
							random.seeds[3]
					),
					"-T", ">", ms.output.file
				)
		cat(call,file="ms_call.txt")
		system(call)
		read.ms.haplotype.matrices(nsam=K*sampled.pops*n.chromo,ndraws=1,ms.output.file=ms.output.file)
}

# generate times and migration rates in ms units
generate.ms.command.line.values <- function(diploid.population.size,locus.size,per.bp.mu,migration.fraction,deep.split,shallow.split){
	ms.command.line.values <- vector("list",length=4)
		names(ms.command.line.values) <- c("theta","m","deep.split","shallow.split")
			ms.command.line.values$theta <- 4*diploid.population.size*per.bp.mu*locus.size
			ms.command.line.values$m <- 4*diploid.population.size*migration.fraction
			ms.command.line.values$deep.split <- deep.split/(4*diploid.population.size)
			ms.command.line.values$shallow.split <- shallow.split/(4*diploid.population.size)
			# ms.command.line.values$time <- 4*diploid.population.size*generations.ago
			#admixture on time scale more recent than 1/(4Nm k) won't spread
	return(ms.command.line.values)
}

check.param <- function(param,K){
	if(length(param != K)){
		param <- rep(param,K)
	}
	return(param)
}

sim.spCor.admix.props <- function(N,K,geoDist){
	#recover()
	spCor <- 5*exp(-geoDist/0.7)
	t.ad.props <- matrix(NA,nrow=N,ncol=K)
	t.ad.props <- apply(t.ad.props,2,function(x){MASS::mvrnorm(1,rep(0,N),Sigma=spCor)})
	t.ad.props <- apply(t.ad.props,2,function(x){x + abs(min(x))})
	lt.ad.props <- exp(t.ad.props)
	sum.ad.props <- rowSums(lt.ad.props)
	ad.props <- apply(lt.ad.props,2,function(x){x/sum.ad.props})
	return(ad.props)
}

admix.freqs <- function(K,admix.props,allele.counts,sample.sizes){
	sample.freqs <- lapply(1:K,function(k){allele.counts[[k]]/sample.sizes[[k]]})
	weighted.freqs <- lapply(1:K,
						function(k){
							apply(sample.freqs[[k]],2,
									function(f){
										f * admix.props[,k]
									})
						})
	admixed.freqs <- Reduce("+",weighted.freqs)
	return(admixed.freqs)	
}

# make dataset for use by, e.g., spatialStructure
generate.spStr.dataset <- function(n.loci,K,sampling.coords,n.chromo,sim.sample.sizes,theta,migration.rate,shallow.split,deep.split){
	#recover()
	#Allele Counts
		migration.rate <- check.param(migration.rate,K)
		shallow.split <- check.param(shallow.split,K)
		data.matrix <- do.call(cbind,
							replicate(n.loci,
								ms(sampling.coords,K,n.chromo,theta,migration.rate,shallow.split,deep.split)))
		sampled.pops <- K*nrow(sampling.coords)
		population.membership <- c()
			for(i in 1:sampled.pops){
				population.membership <- c(population.membership,rep(i,n.chromo))
			}
		allele.counts <- matrix(0,nrow=sampled.pops,ncol=n.loci)
		for(i in 1:sampled.pops){
			allele.counts[i,] <- colSums(data.matrix[which(population.membership==i),,drop=FALSE])
		}	
	#Sample sizes
		sample.sizes <- matrix(n.chromo,nrow=sampled.pops,ncol=n.loci)
	#Broken up by cluster
		cluster.indices <- lapply(1:K,function(k){1:nrow(sampling.coords) + (k-1)*nrow(sampling.coords)})
		cluster.allele.counts <- lapply(cluster.indices,function(x){allele.counts[x,]})
		cluster.sample.sizes <- lapply(cluster.indices,function(x){sample.sizes[x,]})
	#Create admixed sample
		D <- fields::rdist(sampling.coords)
		admix.props <- sim.spCor.admix.props(nrow(sampling.coords),K,D)
		admix.freqs <- admix.freqs(K,admix.props,cluster.allele.counts,cluster.sample.sizes)
		sample.admix.allele.counts <- matrix(rbinom(nrow(sampling.coords)*n.loci,sim.sample.sizes,admix.freqs),nrow=nrow(sampling.coords),ncol=n.loci)
		sample.admix.freqs <- sample.admix.allele.counts/sim.sample.sizes
	#Return sim dataset
#		freqs <- switcharoo.data(sample.admix.allele.counts/sim.sample.sizes)
		freqs <- sample.admix.freqs
		sample.cov <- cov(t(freqs))
		spStr.dataset <- list("freqs.list" = list(
									"cluster.allele.counts" = cluster.allele.counts,
									"cluster.sample.sizes" = cluster.sample.sizes,
									"admix.freqs" = admix.freqs,
									"sample.admix.allele.counts" = sample.admix.allele.counts),
							  "par.list" = list(
							  		"n.loci" = n.loci,
							  		"admix.props" = admix.props,
							  		"n.chromo" = n.chromo,
							  		"theta" = theta,
							  		"migration.rate" = migration.rate,
							  		"shallow.split" = shallow.split,
							  		"deep.split" = deep.split,
							  		"sim.sample.sizes" = sim.sample.sizes),
							  "data.list" = list(
							  		"allele.freqs" = freqs,	
							  		"sample.cov" = sample.cov,
							  		"D" = D,
							  		"sampling.coords" = sampling.coords))
	return(spStr.dataset)
}