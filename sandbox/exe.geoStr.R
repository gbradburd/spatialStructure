require(geoStructure)
load(list.files(pattern="dataset"))
data.block <- list( "space" = TRUE,
					"time" = FALSE,
					"K" = 3,
					"N" = nrow(haak.dataset$geo.coords),
					"L" = 87158,
					"obsSigma" = haak.dataset$sample.cov,
					"geoCoords" = haak.dataset$geo.coords,
					"geoDist" = fields::rdist(haak.dataset$geo.coords),
					"timeCoords" = NULL,
					"timeDist" = NULL,
					"DirichAlpha" = rep(0.01,3),
					"sampleSize" = haak.dataset$sample.size,
					"binVar" = haak.dataset$bin.var)

model.fit <- geoStructure(data.block=data.block,n.chains=1,n.iter=10000,prefix="haak_mod_s_K3")
my.dir <- getwd()
make.all.the.plots(my.dir,my.dir)