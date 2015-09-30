


load("~/Dropbox/InspectorSpaceTime/spatialStructure/datasets/Poplar/trichocarpa_dataset.Robj")
data.block <- list( "space" = TRUE,
					"time" = FALSE,
					"K" = 2,
					"N" = nrow(trichocarpa.data$geo.coords),
					"L" = trichocarpa.data$n.loci,
					"obsSigma" = trichocarpa.data$sample.cov,
					"geoDist" = fields::rdist(trichocarpa.data$geo.coords),
					"timeDist" = NULL,
					"DirichAlpha" = rep(0.1,2),
					"sampleSize" = trichocarpa.data$sample.sizes,
					"binVar" = trichocarpa.data$binom.var)

test2 <- geoStructure(data.block=data.block,n.chains=2,n.iter=10000,prefix="~/desktop/STANTEST")

ad.cov <- apply(extract(test2,"Sigma")$Sigma,c(2,3),mean)
plot(data.block$geoDist,data.block$obsSigma)
points(data.block$geoDist,ad.cov,col="red",pch=20)
plot(data.block$obsSigma,ad.cov) ; abline(0,1,col="red")
