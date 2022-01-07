#simulation draft

library("univariateML")
library("extraDistr")
#Lomax sampling
	
n=250
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

  
nsteps <- 5500 #mcmcサンプリングの数
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
	if(i%%5 == 0){
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

#plot(data)
#plot(data_before)

#autocorr.plot(data_thinned)
autocorr.plot(data_burned)

rejectionRate(data)
geweke.diag(data,frac1=0.1,frac2=0.5)





-------------------------------------------------------------------------------------------
#以下nが1種類(n=100)の場合のアルゴリズム

library("univariateML")
library("extraDistr")

true_beta = 2
true_alpha = 1.5

#メトロポリスアルゴリズム
Metropolis <- function(current) {
  propose <- abs(rnorm(1,current,0.5))
  postOdds <- post_alpha(propose) / post_alpha(current)
  if(post_alpha(current) == 0) propose
  else if(min(postOdds,1) > runif(1,0,1)) propose else current
}



#1種類のサンプル数で、N回分の推定値を格納するベクトル
#これを平均したものがサンプル数nの時の推定値とする
bayes_mres_alpha <- c()
bayes_mses_alpha <- c()

bayes_mres_beta <- c()
bayes_mses_beta <- c()

MLE_mres_alpha <- c()
MLE_mses_alpha <- c()

MLE_mres_beta <- c()
MLE_mses_beta <- c()

N = 100 #1種類のサンプル数につき N回実験して平均を取る
for (j in 1:N){
  x_total <- c()
  l_total <- c()
  n=250　#サンプルサイズ
    
    #Lomaxに従う確率変数をn個生成
    for (k in 1:n){
    	l <- rgamma(1,true_alpha,1)
        x <- rexp(1,l/true_beta)
        x_total <- append(x_total,x)
        l_total <- append(l_total,l)
        }
        
    #対数尤度
    loglikelihood <- function(p,x) {
    	n*log(p[2])-n*log(p[1])-(p[2]+1)*sum(log((p[1]+x)/p[1]))
    	}
    
    #対数尤度最大化
    like_par <- optim(c(10,10),loglikelihood,x = x_total,method = "SANN",control = list(fnscale = -1))
    
    #alphaの事後分布
    post_alpha <- 
    function(alpha) 
    (1/((alpha+1)*alpha^(1/2)*(alpha+2)^(1/2)*(gamma(alpha))^n))*(prod(l_total))^(alpha - 1)


    nsteps <- 5500 #mcmcサンプリングの数
    alpha_mcmcsample <- rep(0,nsteps) #サンプルを格納するベクトル
    alpha_mcmcsample[1] <- 5 #初期値
    burn_in <- 500 
    
    for (i in 1:nsteps){
    	alpha_mcmcsample[i+1] <- Metropolis(alpha_mcmcsample[i])
    	}
    	
    #burn-in
    alpha_mcmc_burned <- alpha_mcmcsample[-(1:burn_in)]
    length(alpha_mcmc_burned)
    
    #thin
    alpha_mcmc_thinned <- c()
    thin = 5
    for (i in 1:(nsteps-burn_in)){
    	if(i%%thin == 0){
    		alpha_mcmc_thinned <- append(alpha_mcmc_thinned, alpha_mcmc_burned[i])
    		}
    	}
    	
    bayes_alpha_estimated <- mean(alpha_mcmc_thinned)
    #betaの事後分布はは逆ガンマ分布に従うので事後期待値は以下の形
    bayes_beta_estimated <- sum(x_total*l_total)/(n-1)    	
    MLE_alpha <- mllomax(x_total)[2]
    MLE_beta <- 1/mllomax(x_total)[1]

    
    bayes_mres_alpha <- append(bayes_mres_alpha, bayes_alpha_estimated*(1/true_alpha))
    bayes_mses_alpha <- append(bayes_mses_alpha, (bayes_alpha_estimated - true_alpha)^2)
        
    bayes_mres_beta <- append(bayes_mres_beta, bayes_beta_estimated*(1/true_beta))
    bayes_mses_beta <- append(bayes_mses_beta, (bayes_beta_estimated - true_beta)^2)
    
    
    MLE_mres_alpha <- append(MLE_mres_alpha, MLE_alpha*(1/true_alpha))
    MLE_mses_alpha <- append(MLE_mses_alpha, (MLE_alpha - true_alpha)^2)

    MLE_mres_beta <- append(MLE_mres_beta, MLE_beta*(1/true_beta))
    MLE_mses_beta <- append(MLE_mses_beta, (MLE_beta - true_beta)^2)

    }       
    
   
    
(1/N)*sum(bayes_mres_alpha)
(1/N)*sum(bayes_mses_alpha)

(1/N)*sum(bayes_mres_beta)
(1/N)*sum(bayes_mses_beta)

(1/N)*sum(MLE_mres_alpha)
(1/N)*sum(MLE_mses_alpha)

(1/N)*sum(MLE_mres_beta)
(1/N)*sum(MLE_mses_beta)

warnings()
 
#確認
plot(alpha_mcmcsample,type="l",col="blue")

library(coda)
alpha_mcmc <- mcmc(alpha_mcmc_thinned)
summary(alpha_mcmc)

plot(alpha_mcmc)



warnings()
-------------------------------------------------------------------------------------------
#以下最終的なプログラム

rm(list=ls())

library("univariateML")
library("extraDistr")
true_beta = 2
true_alpha = 1.5

#メトロポリスアルゴリズム
Metropolis <- function(current) {
  propose <- abs(rnorm(1,current,0.5))
  postOdds <- post_alpha(propose) / post_alpha(current)
  if(post_alpha(current) == 0) propose
  else if(min(postOdds,1) > runif(1,0,1)) propose else current
}


sample_size <- seq(100,200,by = 10)

#実際にplotする値を格納するベクトル
bayes_mre_alpha_all <- c()
bayes_mse_alpha_all <- c()

bayes_mre_beta_all <- c()
bayes_mse_beta_all <- c()

MLE_mre_alpha_all <- c()
MLE_mse_alpha_all <- c()

MLE_mre_beta_all <- c()
MLE_mse_beta_all <- c()



for (n in sample_size){
  #1種類のサンプル数につき N回実験して平均を取る
  N = 2000
  #後で平均をとるために推定値を格納するベクトル
  bayes_mres_alpha <- c()
  bayes_mses_alpha <- c()
  
  bayes_mres_beta <- c()
  bayes_mses_beta <- c()
  
  MLE_mres_alpha <- c()
  MLE_mses_alpha <- c()

　MLE_mres_beta <- c()
　MLE_mses_beta <- c()

  for (j in 1:N){
    x_total <- c()
    l_total <- c()
    
    #Lomaxに従う確率変数をn個生成
    for (k in 1:n){
      l <- rgamma(1,true_alpha,1)
      x <- rexp(1,l/true_beta)
      x_total <- append(x_total,x)
      l_total <- append(l_total,l)
    }
    
    #alphaの事後分布
    post_alpha <- 
      function(alpha) 
        (1/((alpha+1)*alpha^(1/2)*(alpha+2)^(1/2)*(gamma(alpha))^n))*(prod(l_total))^(alpha - 1)
    
    
    nsteps <- 5500 #mcmcサンプリングの数
    alpha_mcmcsample <- rep(0,nsteps) #サンプルを格納するベクトル
    alpha_mcmcsample[1] <- 5 #初期値
    burn_in <- 500 
    
    for (i in 1:nsteps){
      alpha_mcmcsample[i+1] <- Metropolis(alpha_mcmcsample[i])
    }
    
    #burn-in
    alpha_mcmc_burned <- alpha_mcmcsample[-(1:burn_in)]
    length(alpha_mcmc_burned)
    
    #thin
    alpha_mcmc_thinned <- c()
    thin = 5
    for (i in 1:(nsteps-burn_in)){
      if(i%%thin == 0){
        alpha_mcmc_thinned <- append(alpha_mcmc_thinned, alpha_mcmc_burned[i])
      }
    }
    
    #推定値
    bayes_alpha_estimated <- mean(alpha_mcmc_thinned)
    #betaの事後期待値は以下の形
    bayes_beta_estimated <- sum(x_total*l_total)/(n-1)
    MLE_alpha <- mllomax(x_total)[2]
    MLE_beta <- 1/mllomax(x_total)[1]
    
    
    
    bayes_mres_alpha <- append(bayes_mres_alpha, bayes_alpha_estimated*(1/true_alpha))
    bayes_mses_alpha <- append(bayes_mses_alpha, (bayes_alpha_estimated - true_alpha)^2)
    
    bayes_mres_beta <- append(bayes_mres_beta, bayes_beta_estimated*(1/true_beta))
    bayes_mses_beta <- append(bayes_mses_beta, (bayes_beta_estimated - true_beta)^2)
    
    MLE_mres_alpha <- append(MLE_mres_alpha, MLE_alpha*(1/true_alpha))
    MLE_mses_alpha <- append(MLE_mses_alpha, (MLE_alpha - true_alpha)^2)
    
    MLE_mres_beta <- append(MLE_mres_beta, MLE_beta*(1/true_beta))
    MLE_mses_beta <- append(MLE_mses_beta, (MLE_beta - true_beta)^2)
    
  }
  bayes_mre_alpha_all <- append(bayes_mre_alpha_all, (1/N)*sum(bayes_mres_alpha))
  bayes_mse_alpha_all <- append(bayes_mse_alpha_all, (1/N)*sum(bayes_mses_alpha))
  
  bayes_mre_beta_all <- append(bayes_mre_beta_all, (1/N)*sum(bayes_mres_beta))
  bayes_mse_beta_all <- append(bayes_mse_beta_all, (1/N)*sum(bayes_mses_beta))
  
  MLE_mre_alpha_all <- append(MLE_mre_alpha_all, (1/N)*sum(MLE_mres_alpha))
  MLE_mse_alpha_all <- append(MLE_mse_alpha_all, (1/N)*sum(MLE_mses_alpha))
  
  MLE_mre_beta_all <- append(MLE_mre_beta_all,(1/N)*sum(MLE_mres_alpha))
  MLE_mse_beta_all <- append(MLE_mse_beta_all,(1/N)*sum(MLE_mses_alpha))
}

par(mfrow=c(2,2)) 
plot(sample_size, bayes_mre_alpha_all, ylim = c(0.9, 1.2), ylab = "", type="l", col = "yellow")
par(new=T)
plot(sample_size, MLE_mre_alpha_all, ylim = c(0.9, 1.2), ylab = "MRE(alpha)", type="l", col = "blue")
abline(h=1,lty=2)

plot(sample_size, bayes_mse_alpha_all, ylim = c(0.0, 1.0), ylab = "", type="l", col = "yellow")
par(new=T)
plot(sample_size, MLE_mse_alpha_all, ylim = c(0.0, 1.0), ylab = "MSE(alpha)", type="l", col = "blue")
abline(h=0,lty=2)

plot(sample_size, bayes_mre_beta_all, ylim = c(0.95, 1.2),  ylab = "", type="l", col = "yellow") 
par(new=T)
plot(sample_size, MLE_mre_beta_all, ylim = c(0.95, 1.2),  ylab = "MRE(beta)", type="l", col = "blue")
abline(h=1,lty=2)

plot(sample_size, bayes_mse_beta_all, ylim = c(0.0, 1.0), ylab = "", type="l", col = "yellow")
par(new=T)
plot(sample_size, MLE_mse_beta_all, ylim = c(0.0, 1.0),  ylab = "MSE(beta)", type="l", col = "blue")
abline(h=0,lty=2)
