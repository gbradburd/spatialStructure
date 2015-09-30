setwd("~/Dropbox/InspectorSpaceTime/spatialStructure/code")
require(devtools)
require(rmarkdown)
if(FALSE){
setup() #ONLY RUN setup() once!
# edit description file
# put R code in R file
# add dependencies
oneK.stan.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/basic.txt"),collapse="\n")
time.oneK.stan.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/time.txt"),collapse="\n")
space.oneK.stan.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/space.txt"),collapse="\n")
space.time.oneK.stan.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/space_time.txt"),collapse="\n")
multiK.stan.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/multiK.txt"),collapse="\n")
time.multiK.stan.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/time_multiK.txt"),collapse="\n")
space.multiK.stan.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/space_multiK.txt"),collapse="\n")
space.time.multiK.block <- paste(readLines("~/Dropbox/InspectorSpaceTime/spatialStructure/code/stan_model_blocks/space_time_multiK.txt"),collapse="\n")
use_data(oneK.stan.block,
 		 time.oneK.stan.block,
		 space.oneK.stan.block,
		 space.time.oneK.stan.block,
		 multiK.stan.block,
		 time.multiK.stan.block,
		 space.multiK.stan.block,
		 space.time.multiK.block,
		 internal=TRUE)
}

# edit NAMESPACE
load_all()
# write documentation
document()
# after updating documentation, run document() again, then run
# R CMD Rdconv -t html manfile.Rd -o out.html
#have to add texbin to R's path so that it can find pdflatex
#	to build the vignette
Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/usr/texbin",sep=":"))
check(build_args="--resave-data")
build(path="~/Dropbox/InspectorSpaceTime/spatialStructure/code",args="--resave-data")
install()



devtools::install_github("gbradburd/SpaceMix",build_vignettes=TRUE)
vignette("spacemix_vignette")


help(make.spacemix.map.list)

help(run.spacemix.analysis)

help(fade.admixture.source.points)

help(spacemix.procrustes)