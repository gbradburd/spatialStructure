################################
################################
#	notes on beta-distributed epsilon
################################
################################


beta.pars <- seq(1e-3,10,length.out=1000)
n.loci <- 1e5
eps <- lapply(beta.pars,function(x){
				rbeta(n.loci,x,x)
				})
				
var.eps <- lapply(eps,function(x){var(x)})
binVar <- lapply(eps,function(x){mean(x*(1-x))})

par(mfrow=c(1,3))
plot(beta.pars,unlist(var.eps))
	points(beta.pars,beta.var(beta.pars),col="red",pch=20)
	points(beta.pars,beta.var2(beta.pars),col="green",pch=20)
plot(beta.pars,unlist(binVar))
	points(beta.pars,0.25-beta.var(beta.pars),col="blue",pch=20)
plot(unlist(var.eps),unlist(binVar))

beta.var <- function(a){
	beta.var <- (a^2)/((2*a)^2*(2*a+1))
	return(beta.var)
}

beta.var2 <- function(a){
	beta.var <- 1/(8*a+4)
	return(beta.var)
}

# mean(x*(1-x)) == mean(x - x^2) == mean(x) - mean(x^2)
# ==mean(x) -(Var(x) + mean(x)*mean(x))
# 0.5 - (Var(x) + 0.25)

