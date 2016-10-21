load("/Users/gburd/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/norm_lnL/test/spStr.dataset.Robj")
source("~/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/norm_lnL/geoStructure.R")
model.block <- paste0(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/sandbox/norm_lnL/stan_model_blocks/space_multiK.txt",warn=FALSE),collapse="\n")
std.list <- standardize.freqs(spStr.dataset$data.list$allele.freqs)

data.block <- list( "spatial" = TRUE,
					"K" = 2,
					"N" = nrow(spStr.dataset$data.list$sampling.coords),
					"L" = std.list$n.loci,
					"obsSigma" = std.list$sample.cov,
					"geoDist" = spStr.dataset$data.list$D,
					"stdErrs" = std.list$std.errs,
					"sampleSize" = rep(10,nrow(spStr.dataset$data.list$sampling.coords)))
	require(rstan)
	# fit <- stan(model_code = model.block, data =data.block, iter=1e4, chains = 1,thin=1e4/500)
	
	ml.fits <- lapply(1:10,function(i){
					optimizing(object = stan_model(model_code= model.block),
								data = data.block,
								iter=1e2)})

	fit.ml <- optimizing(object = stan_model(model_code= model.block),
							data = data.block,
							iter=1e5)
	save(spStr.dataset,data.block,fit,file=paste0(dir,"/test_output.Robj",collapse=""))
	
	par.cov <- get.ML.par.cov(fit.ml,data.block)
	geoStr.results <- get.geoStructure.ML.results(fit.ml,data.block)
	
	plot(data.block$obsSigma,geoStr.results$point$par.cov) ; abline(0,1,col="red")
	ml.results <- lapply(ml.fits,function(x){get.geoStructure.ML.results(x,data.block)})
	par(mfrow=c(2,5))
	lapply(ml.results,function(x){plot(data.block$obsSigma,x$point$par.cov) ; abline(0,1,col="red")})