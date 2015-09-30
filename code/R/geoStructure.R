validate.list <- function(data.block){
	if(!"space" %in% names(data.block)){
		stop("\nUser must specify a \"space\" option\n\n")
	}
	if(!"time" %in% names(data.block)){
		stop("\nUser must specify a \"time\" option\n\n")
	}
	if(!"N" %in% names(data.block)){
		stop("\nUser must specify a \"N\" option\n\n")
	}
	if(!"K" %in% names(data.block)){
		stop("\nUser must specify a \"K\" option\n\n")
	}
	if(!"L" %in% names(data.block)){
		stop("\nUser must specify a \"L\" option\n\n")
	}
	if(!"obsSigma" %in% names(data.block)){
		stop("\nUser must specify a \"obsSigma\" option\n\n")
	}
	if(!"geoDist" %in% names(data.block)){
		stop("\nUser must specify a \"geoDist\" option\n\n")
	}
	if(!"timeDist" %in% names(data.block)){
		stop("\nUser must specify a \"timeDist\" option\n\n")
	}
	if(!"DirichAlpha" %in% names(data.block)){
		stop("\nUser must specify a \"DirichAlpha\" option\n\n")
	}
	if(!"sampleSize" %in% names(data.block)){
		stop("\nUser must specify a \"sampleSize\" option\n\n")
	}
	if(!"binVar" %in% names(data.block)){
		stop("\nUser must specify a \"binVar\" option\n\n")
	}
	return(invisible("list elements validated"))	
}

validate.n.samples <- function(data.block){
	n.samples <- data.block$N
	n.samples <- c(data.block$N,nrow(data.block$obsSigma))
	n.samples <- c(n.samples,length(data.block$sampleSize))
	if(!is.null(data.block$geoDist)){
		n.samples <- c(n.samples,nrow(data.block$geoDist))
	}
	if(!is.null(data.block$timeDist)){
		n.samples <- c(n.samples,nrow(data.block$timeDist))
	}
	if(length(unique(n.samples)) > 1){
		stop("\nthe number of samples is not consistent 
				across entries in the data.block\n\n")
	}
	return(invisible("n.samples validated"))
}

validate.model <- function(data.block){
	if(data.block$space){
		if(is.null(data.block$geoDist)){
			stop("\nyou have specified a spatial model,
				  but you have not specified a matrix 
				  of pairwise geographic distances")
		}
	}
	if(data.block$time){
		if(is.null(data.block$timeDist)){
			stop("\nyou have specified a spatial model,
				  but you have not specified a matrix 
				  of pairwise geographic distances")
		}
	}
	n.clusters <- data.block$K
	if(n.clusters > 1){
		if(n.clusters != length(data.block$DirichAlpha)){
			stop("\nthe vector of Dirichlet concentration prior parameters must 
					have a length equal to the number of specified clusters\n\n")
		}
	}
	return(invisible("model validated"))
}


make.data.block.S3 <- function(data.block){
	data.block <- data.block
	class(data.block) <- "data.block"
	return(data.block)
}

print.data.block <- function(data.block){
	print(str(data.block,max.level=1))
}

validate.data.block <- function(data.block){
	message("\nchecking data.block\n")
	tmp <- validate.list(data.block)
	tmp <- validate.n.samples(data.block)
	message(sprintf("reading %s samples",data.block$N))
	message(sprintf("reading %s loci",data.block$L))
	message("\nchecking specified model\n")
	tmp <- validate.model(data.block)
	message(sprintf("reading %s cluster(s)",data.block$K))
	if(data.block$space & !data.block$time){
		message(sprintf("user has specified a spatial model"))
	}
	if(!data.block$space & data.block$time){
		message(sprintf("user has specified a temporal model"))
	}
	if(data.block$space & data.block$time){
		message(sprintf("user has specified a spatiotemporal model"))
	}
	if(!data.block$space & !data.block$time){
		message(sprintf("user has specified a purely discrete model"))
	}
	data.block <- make.data.block.S3(data.block)
	return(data.block)
}

make.stan.code.block <- function(space,time,n.clusters){
	stan.code.block.name <- "stan.block"
	if(n.clusters == 1){
		stan.code.block.name<- paste("oneK.",stan.code.block.name,sep="")
	}
	if(n.clusters > 1){
		stan.code.block.name<- paste("multiK.",stan.code.block.name,sep="")	
	}
	if(time){
		stan.code.block.name <- paste("time.",stan.code.block.name,sep="")
	}
	if(space){
		stan.code.block.name <- paste("space.",stan.code.block.name,sep="")
	}
	return(get(stan.code.block.name))
}

geoStructure <- function(data.block,n.chains,n.iter,prefix){
	#validate on data block
	data.block <- validate.data.block(data.block)
		save(data.block,file=paste(prefix,"_data.block.Robj"))
	#validate on model specification
	#make.stan.code.block
	stan.block <- make.stan.code.block(data.block$space,data.block$time,data.block$K)
		#write stan block to file
	#run model
	#put stan in tryCatch, email me
	require(rstan)
	model.fit <- stan(model_code = stan.block,
						data = data.block,
						iter = n.iter,
						chains = n.chains)
	#save fit obj
	save(model.fit,file=paste(prefix,"_model.fit.Robj"))
	#return fit obj
	return(model.fit)
}