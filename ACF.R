#論文の手法
rm(list=ls())
library("univariateML")
library("extraDistr")

#Lomax sampling
n=150
true_beta = 2
true_alpha = 1.5
x_total <- c()
l_total <- c()
for (i in 1:n){
  l <- rgamma(1,true_alpha,1)
  x <- rexp(1,l/true_beta)
  x_total <- append(x_total,x)
  l_total <- append(l_total,l)
}

#メトロポリスアルゴリズム

post_alpha <- #alphaの事後分布
  function(alpha) (1/((alpha+1)*alpha^(1/2)*(alpha+2)^(1/2)*(gamma(alpha))^n))*
  (prod(l_total))^(alpha - 1)


Metropolis <- function(current) {
  propose <- abs(rnorm(1,current,1))
  postOdds <- post_alpha(propose) / post_alpha(current)
  if(is.nan(postOdds)) current
  else if(min(postOdds,1) > runif(1,0,1)) propose else current
}


nsteps <- 15500 #mcmcサンプリングの数
mcmcsample <- rep(0,nsteps) #要素が欠損しているnsteps次元ベクトルを生成
mcmcsample[1] <- 5
burn_in <- 500

for (i in 1:(nsteps-1)){
  mcmcsample[i+1] <- Metropolis(mcmcsample[i])
}

#確認
length(mcmcsample)
plot(mcmcsample,type="l",col="blue")

#burn-in
mcmc_burned <- mcmcsample[-(1:burn_in)]
length(mcmc_burned)

#thin
mcmc_thinned <- c()

for (i in 1:(nsteps-burn_in)){
  if(i%%15 == 0){
    mcmc_thinned <- append(mcmc_thinned, mcmc_burned[i])
  }
}
#mcmc_thinned
#length(mcmc_thinned)

#plot(mcmc_burned,type="l",col="blue")
#plot(mcmc_thinned,type="l",col="blue")

#mcmcへ変換
library(coda)
data_before <- mcmc(mcmcsample)
data_burned <- mcmc(mcmc_burned)
data_thinned <- mcmc(mcmc_thinned) 

#plot(data)
plot(data_before)
plot(data_burned)

pdf("acf.pdf")
acf(mcmc_burned,ci=0,main="",20)
dev.off()
#-------------------------------------------------------------------------
#notlambda


n=150
true_beta = 2
true_alpha = 1.5
x_total <- c()
for (i in 1:n){
  l <- rgamma(1,true_alpha,1)
  x <- rexp(1,l/true_beta)
  x_total <- append(x_total,x)
}


post <- #事後分布
  function(beta, alpha) 
    ((alpha^(n-1/2))*beta^(-n-1))/((alpha + 1)*(alpha + 2)^(1/2))*(prod(1+x_total/beta))^(-(alpha + 1))

nsteps <- 11000 #mcmcサンプリングの数
beta_mcmc <- rep(0,nsteps) #要素が欠損しているnsteps次元ベクトルを生成
alpha_mcmc <- rep(0,nsteps)
beta_mcmc[1] <- 10
alpha_mcmc[1] <- 10
burn_in <- 1000


for (i in 1:(nsteps-1)){
  propose_beta <- abs(rnorm(1,beta_mcmc[i],1))
  propose_alpha <- abs(rnorm(1,alpha_mcmc[i],1))
  postOdds <- post(propose_beta, propose_alpha) / post(beta_mcmc[i], alpha_mcmc[i])
  if(is.nan(postOdds)){
    beta_mcmc[i+1] <- beta_mcmc[i]
    alpha_mcmc[i+1] <- alpha_mcmc[i]
  }
  else if(min(postOdds,1) > runif(1,0,1)){
    beta_mcmc[i+1] <- propose_beta
    alpha_mcmc[i+1] <- propose_alpha
  }else {
    beta_mcmc[i+1] <- beta_mcmc[i]
    alpha_mcmc[i+1] <- alpha_mcmc[i]
  }
}

#length(beta_mcmc)
#length(alpha_mcmc)
#plot(beta_mcmc,type="l", col="blue")
#plot(alpha_mcmc,type="l", col="blue")


#burn-in
burned_beta <- beta_mcmc[-(1:burn_in)]
burned_alpha <- alpha_mcmc[-(1:burn_in)]
length(burned_beta)
length(burned_alpha)

#thin
thinned_beta <- c()
thinned_alpha <- c()

for (i in 1:(nsteps-burn_in)){
  if(i%%10 == 0){
    thinned_beta <- append(thinned_beta, burned_beta[i])
    thinned_alpha <- append(thinned_alpha, burned_alpha[i])
  }
}

#mcmcへ変換
library(coda)
data_before_beta <- mcmc(beta_mcmc)
data_before_alpha <- mcmc(alpha_mcmc)
data_burned_beta <- mcmc(burned_beta)
data_burned_alpha <- mcmc(burned_alpha)
data_thinned_beta <- mcmc(thinned_beta)
data_thinned_alpha <- mcmc(thinned_alpha)

plot(data_before_alpha)
plot(data_before_beta)
plot(data_burned_alpha)
plot(data_burned_beta)

pdf("acf_notlambda.pdf")
par(mfrow = c(2,1))
acf(burned_beta,ci=0,main="",ylab="ACF(beta)",40)
acf(burned_alpha,ci=0,main="",ylab="ACF(alpha)",40)
dev.off()



#-------------------------------------------------------
#mle-bayes




library("univariateML")
library("extraDistr")
#Lomax sampling

n=150
true_beta = 2
true_alpha = 1.5
x_total <- c()
lambdas <- c()
for (i in 1:n){
  l <- rgamma(1,true_alpha,1)
  x <- rexp(1,l/true_beta)
  x_total <- append(x_total,x)
}


#最尤推定の結果
MLE_alpha <- mllomax(x_total)[2]
MLE_beta <- 1/mllomax(x_total)[1]

for (i in 1:n){
  lambdas[i] <- rgamma(1, MLE_alpha + 1, 1 + x_total[i]/MLE_beta)
}

#alphaの事後分布
post_alpha <- 
  function(alpha) 
    (1/((alpha+1)*alpha^(1/2)*(alpha+2)^(1/2)*(gamma(alpha))^n))*(prod(lambdas))^(alpha - 1)


Metropolis <- function(current) {
  propose <- abs(rnorm(1,current,1))
  postOdds <- post_alpha(propose) / post_alpha(current)
  if(is.nan(postOdds)) current
  else if(min(postOdds,1) > runif(1,0,1)) propose else current
}


nsteps <- 15500 #mcmcサンプリングの数
mcmcsample <- rep(0,nsteps) #要素が欠損しているnsteps次元ベクトルを生成
mcmcsample[1] <- 5
burn_in <- 500

for (i in 1:(nsteps-1)){
  mcmcsample[i+1] <- Metropolis(mcmcsample[i])
}

#burn-in
mcmc_burned <- mcmcsample[-(1:burn_in)]
length(mcmc_burned)

#thin
mcmc_thinned <- c()

for (i in 1:(nsteps-burn_in)){
  if(i%%15 == 0){
    mcmc_thinned <- append(mcmc_thinned, mcmc_burned[i])
  }
}
mcmc_thinned
length(mcmc_thinned)

#mcmcへ変換
library(coda)
data_before <- mcmc(mcmcsample)
data_burned <- mcmc(mcmc_burned)
data_thinned <- mcmc(mcmc_thinned) 

pdf("acf_mle-bayes.pdf")
acf(data_burned,ci=0,main="",20)
dev.off()