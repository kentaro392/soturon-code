#lambda使わない

n=150
true_beta = 2
true_alpha = 1.5
x_total <- c()
for (i in 1:n){
  l <- rgamma(1,true_alpha,1)
  x <- rexp(1,l/true_beta)
  x_total <- append(x_total,x)
}


#メトロポリスアルゴリズム

post <- #事後分布
  function(beta, alpha) 
    ((alpha^(n-1/2))*beta^(-n-1))/((alpha + 1)*(alpha + 2)^(1/2))*(prod(1+x_total/beta))^(-(alpha + 1))

nsteps <- 31000 #mcmcサンプリングの数
beta_mcmc <- rep(0,nsteps) #要素が欠損しているnsteps次元ベクトルを生成
alpha_mcmc <- rep(0,nsteps)
beta_mcmc[1] <- 5
alpha_mcmc[1] <- 5
burn_in <- 1000


for (i in 1:(nsteps-1)){
  propose_beta <- abs(rnorm(1,beta_mcmc[i],1))
  propose_alpha <- abs(rnorm(1,alpha_mcmc[i],1))
  postOdds <- post(propose_beta, propose_alpha) / post(beta_mcmc[i], alpha_mcmc[i])
  if (min(postOdds,1) > runif(1,0,1)) {
    beta_mcmc[i+1] <- propose_beta
    alpha_mcmc[i+1] <- propose_alpha	
  }else {
    beta_mcmc[i+1] <- beta_mcmc[i]
    alpha_mcmc[i+1] <- alpha_mcmc[i]
  }
}

length(beta_mcmc)
length(alpha_mcmc)
plot(beta_mcmc,type="l", col="blue")
plot(alpha_mcmc,type="l", col="blue")


#burn-in
burned_beta <- beta_mcmc[-(1:burn_in)]
burned_alpha <- alpha_mcmc[-(1:burn_in)]
length(burned_beta)
length(burned_alpha)

#thin
thinned_beta <- c()
thinned_alpha <- c()

for (i in 1:(nsteps-burn_in)){
  if(i%%30 == 0){
    thinned_beta <- append(thinned_beta, burned_beta[i])
    thinned_alpha <- append(thinned_alpha, burned_alpha[i])
  }
}

plot(thinned_beta,type="l")
plot(thinned_alpha,type="l")

mean(thinned_beta)
mean(thinned_alpha)

#mcmcへ変換
library(coda)
before_beta <- mcmc(beta_mcmc)
before_alpha <- mcmc(alpha_mcmc)
data_beta <- mcmc(thinned_beta)
data_alpha <- mcmc(thinned_alpha)

summary(data_beta)
summary(data_alpha)

plot(data_beta)
plot(data_alpha)

rejectionRate(before_beta)
autocorr.plot(data_beta)
autocorr.plot(data_alpha)

geweke.diag(data_beta)
geweke.diag(data_alpha) 





#----------------------------------------------------------

  