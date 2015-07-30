
###Simple two param covar matrix model.

library("rstan")
library("MASS")

ex_model_code <-
'
data {
  int<lower=2> K; 		
  int<lower=K+2> L;
  matrix[K,K] obsSigma; // observed covar
}
parameters {
real<lower=0>  alpha;
real<lower=0>  beta;
}
transformed parameters {
  matrix[K,K] Sigma;
  for (i in 1:K){
  	 for (j in 1:K){
            Sigma[i, j] <- alpha;
         }
         }
    for (i in 1:K){
       Sigma[i, i]<- Sigma[i, i] + beta;
	}
}
model {
alpha ~ exponential(1);
beta ~ exponential(1);
(L*obsSigma) ~ wishart(L,Sigma);
}
'

myData <- list(K = 2L,L=100L,
                    obsSigma = matrix(c(2.5,1.5,1.5,2.5),byrow=TRUE,nrow=2)
                    )

fit <- stan(model_code = ex_model_code, data =myData, iter=1000, chains = 4)

#####  Simple Single IBD population






sp_model_code <-
'
data {
  int<lower=2> N; 		
  int<lower=N+2> L;
  matrix[N,N] obsSigma; // observed covar
  matrix[N, N] myDists;  
}
parameters {
	real<lower=0>  alpha0;
	real<lower=0>  alpha1;
	real<lower=0>  alpha2;
  real<lower=0> nugget[N]; // s.e. of effect estimates 
	
}
transformed parameters {
  matrix[N,N] Sigma;
  for (i in 1:N){
  	 for (j in 1:N){
            Sigma[i, j] <- alpha0 * exp(-(alpha1* myDists[i, j] )^alpha2 );
         }
         }
    for (i in 1:N){
      Sigma[i, i]<- Sigma[i, i] + nugget[i];
	}
}
model {
	alpha0 ~ exponential(1);
	alpha1 ~ exponential(1);
	alpha2 ~ exponential(1);
	nugget ~ exponential(1);
(L*obsSigma) ~ wishart(L,Sigma);
}
'

library("fields")
N<-20; L=1000
x<-runif(N)
y<-runif(N)
my.dists<-rdist(cbind(x,y))

alpha2<-1; alpha0<-1;alpha1<-2; 
my.covar<- alpha0 *exp(-(alpha1*my.dists)^alpha2)

nuggets<-rexp(N)
diag(my.covar)<- diag(my.covar) + nuggets

mv.data<-mvrnorm(n = L, mu=rep(0,N), Sigma=my.covar)

myData <- list(N=N,L=L,                                      #K = 10L,L=1000L,
                    obsSigma = cov(mv.data),
                    myDists = my.dists
                    )
fit <- stan(model_code = sp_model_code, data =myData, iter=1000, chains = 4)

my.params<-extract(fit,pars=c("alpha0","alpha1","alpha2","nugget"),permuted=FALSE)

my.model <- stan_model(model_code = sp_model_code)
map.fit <- optimizing( my.model , data =myData)


###### Spatial  STRUCTURE code  

library("fields")
N<-50; L=1000
x<-runif(N)
y<-runif(N)
my.dists<-rdist(cbind(x,y))

alpha2<-1; alpha0<-1;alpha1<-2; 
my.covar.1<- alpha0 *exp(-(alpha1*my.dists)^alpha2)

my.ws<-rbeta(N,0.9,0.9)
my.covar<- (my.ws %*% t(my.ws)) * my.covar.1  +  ((1-my.ws) %*% t(1-my.ws)) * my.covar.1

nuggets<-rexp(N,rate=10)
diag(my.covar)<- diag(my.covar) + nuggets

plot(c(my.dists),c(my.covar))
mv.data<-mvrnorm(n = L, mu=rep(0,N), Sigma=my.covar)

spStruct_model_code <-
'
data {
	int<lower=1> K;
	int<lower=2> N; 		
	int<lower=N+2> L;
  matrix[N,N] obsSigma; // observed covar
  matrix[N, N] myDists;  
  vector[K]  DirichAlpha; 
}
parameters {
	real<lower=0>  alpha0[K];
	real<lower=0>  alpha1[K];
	real<lower=0, upper=2>  alpha2[K];
	real<lower=0> mu[K];
	real<lower=0> gamma;
  real<lower=0> nugget[N]; // s.e. of effect estimates 
	simplex[K]  w[N];    //every individual (N in tot) has a K simplex (i.e. K clusters) 
}
transformed parameters {
  matrix[N,N] Sigma;
  
	 for (i in 1:N){
	  	 for (j in 1:N){
	  	 	  Sigma[i, j] <- gamma;
			  for(k in 1:K){
	            Sigma[i, j] <- Sigma[i, j]  +  w[i,k] * w[j,k] * ( alpha0[k] * exp(-(alpha1[k]* myDists[i, j] )^alpha2[k] ) + mu[k]);  
	         }
	         if(i==j){
		         Sigma[i, i]<- Sigma[i, i] + nugget[i];
	         }
	         }
	   }
}
model {
	alpha0 ~ exponential(1);
	alpha1 ~ exponential(1);
	alpha2 ~ exponential(1);
	nugget ~ exponential(1);
	mu ~ exponential(1);
	gamma ~ exponential(1);
	for(i in 1:N) w[i] ~ dirichlet(DirichAlpha);
	(L*obsSigma) ~ wishart(L,Sigma);
}
'
K<-2
myData <- list(K=K, N=N,L=L,                                      #K = 10L,L=1000L,
                    obsSigma = cov(mv.data),
                    myDists = my.dists,
                    DirichAlpha = rep(0.1,K)
                    )
                    
                    
fit <- stan(model_code = spStruct_model_code, data =myData, iter=1000, chains = 4)
my.params<-extract(fit,par = c("gamma","mu"),permuted=FALSE)
pairs(my.params[,1,])
hist(my.params[,4,2])


##optimizing code
my.model <- stan_model(model_code = spStruct_model_code)
x<-rbeta(N,1,1)
my.init<-list(alpha0=rep(1,K),alpha1=rep(1,K), alpha2=rep(1,K), mu=rep(0.001,K),gamma=0.001, w=matrix(c(x,1-x),nrow=N,ncol=K), nugget=rep(0.001,N))
map.fit <- optimizing( my.model , data =myData,init=my.init)
map.fit$par[grep("w",names(map.fit$par))]