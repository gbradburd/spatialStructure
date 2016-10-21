# Load package
require(msarg)

# change path so that I can call ms
Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/Applications/ms.folder/msdir",sep=":"))

# function for reading ms output
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

# function for generating n.loci X segsites to create a full dataset
create.dataset <- function(ga,sample.config,n.loci){
	n.chromos <- nrow(sample.config)
	genotypes <- matrix(NA,nrow=n.chromos,ncol=n.loci)
	for(i in 1:n.loci){
		run_ms(ga, outdir="ms_output", nsamp=sample.config, theta=0.01, nreps=1)
		genos <- read.ms.haplotype.matrices(nrow(sample.config)*2,1,"ms_output/msoutput.txt")[[1]]
		genotypes[,i] <- genos[seq(1,n.chromos*2,by=2)] + genos[seq(2,n.chromos*2,by=2)]
	}
	dataset <- list("genos" = genotypes,
					"n.loci" = n.loci,
					"sample.config" = sample.config)
	return(dataset)
}

# specify simulation parameters
nrow <- 16               # width of grid
ncol <- 16                # height of grid
overlap.size <- 2        # with of the overlap (in the middle)
admix.rate <- 0.05
mig.rate <- 1

ga <- grid_array(nlayers=2,nrow=nrow,ncol=ncol,N=1,mig.rate=mig.rate,admix.rate=admix.rate) 

mask.1 <- mask.2 <- matrix(1,nrow=nrow(ga),ncol=ncol(ga))
mask.1[row(mask.1) > nrow(ga)/2+overlap.size/2] <- 0
mask.2[row(mask.1) < nrow(ga)/2-overlap.size/2] <- 0
ga <- modify_grid_layer( ga, layer=1, dN=mask.1 )
ga <- modify_grid_layer( ga, layer=2, dN=mask.2 )
#plot(ga)

n.loci <- 1e4

sample.config <- sample_locations( ga, n=30, each=2 )
sample.cols <- rainbow(nrow(sample.config))

# simulate dataset
sim.dataset <- create.dataset(ga,sample.config,n.loci)
save(sim.dataset,file="sim.dataset.Robj")

# plot output
sample.freqs <- sim.dataset$genos/2
sample.cov <- cov(t(sample.freqs))
coords <- sim.dataset$sample.config[,1:2]
geoDist <- fields::rdist(coords)
pc.obj <- eigen(sample.cov)$vectors
pc1 <- pc.obj[,2] + abs(min(pc.obj[,2]))
cols <- hsv(h=pc1,s=1,v=1)

pdf("sim_check_plots.pdf",width=9,height=5)
par(mfrow=c(1,2))
plot(geoDist,sample.cov)
plot(coords,col=cols,pch=19)
dev.off()