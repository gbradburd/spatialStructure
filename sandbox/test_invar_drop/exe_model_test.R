source("model_test_lib.R")
model.block.file <- "stan_block.txt"

for(i in 1:100){
	dir <- sprintf("testrun_%s",i)
	dir.create(dir)
	run.test(20,1e4,model.block.file,dir,n.iter=1e4)
}
