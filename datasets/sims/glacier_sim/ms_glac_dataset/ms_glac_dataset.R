source("msarg.R")

nrow <- 100               # width of grid
ncol <- 40                # height of grid
glacier.width <- 70        # with of the glacier (in the middle)
mig.rate <- 1             # migration rate between neighbors in the grid
growth.time <- 200          # units of time (in Ne) it took to re-expand post-glacier
glacier.end <- 500          # units of time ago (in Ne) that glaciation ended
glacier.begin <- 1000        # units of time ago (in Ne) that glaciation began

dem <- demography( grid_array(nlayers=1,nrow=nrow,ncol=ncol,N=1,mig.rate=mig.rate) )
plot(dem[[1]])

mask <- matrix(1,nrow=nrow(dem),ncol=ncol(dem))
mask[ abs(row(mask)-nrow(dem)/2) < glacier.width/2] <- 0
dem <- add_to_demography( dem, tnew=glacier.end, fn=modify_grid_layer, layer=1, dN=mask )
plot(dem[[2]])

speed <- 1.2 * ( glacier.width ) / growth.time 
width <- 10
dem <- logistic_interpolation( dem, 
    t.end=glacier.end-growth.time,   # time ago expansion ended
    t.begin=glacier.end,             # time ago exapnsion began
    nsteps=10,                        # number of time steps to pass to ms
    speed=speed, width=width )       # parameters; could also pass sigma= and r= .


dem <- add_to_demography( dem, tnew=glacier.begin, pop=1 )
plot(dem,do.arrows=FALSE)

sample.config <- sample_locations( dem, n=30, each=2 )
sample.cols <- rainbow(nrow(sample.config))

n.loci <- 1e3
locus_names <- lapply(seq_along(1:n.loci),function(i){
						ms.output <- run_ms( dem, nsamp=sample.config, segsites=1, theta=0.001 )
						})
ms_scenario <- list("dem" = dem,"sample.config"=sample.config,"locus_names"=locus_names)
save(ms_scenario,file="ms_glac_scenario.Robj")

