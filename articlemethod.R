#論文のアルゴリズム

rm(list=ls())

library("univariateML")
library("extraDistr")

true_beta = 2
true_alpha = 1.5

Metropolis <- function(current) {
  propose <- abs(rnorm(1,current,1))
  postOdds <- post_alpha(propose) / post_alpha(current)
  if(is.nan(postOdds)) current
  else if(min(postOdds,1) > runif(1,0,1)) propose else current
}

sample_size <- seq(100,250,by = 10)

bayes_mre_alpha_all <- c()
bayes_mse_alpha_all <- c()

bayes_mre_beta_all <- c()
bayes_mse_beta_all <- c()

MLE_mre_alpha_all <- c()
MLE_mse_alpha_all <- c()

MLE_mre_beta_all <- c()
MLE_mse_beta_all <- c()


for (n in sample_size){
  N = 1000 #1種類のサンプル数につき N回実験して平均を取る
  
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
    
    post_alpha <- 
      function(alpha) 
        (1/((alpha + 1)*alpha^(1/2)*(alpha + 2)^(1/2)*(gamma(alpha))^n))*(prod(l_total))^(alpha - 1)
    
    
    nsteps <- 15500 #mcmcサンプリングの数
    alpha_mcmcsample <- rep(0,nsteps) #サンプルを格納するベクトル
    alpha_mcmcsample[1] <- 5 #初期値
    burn_in <- 500 
    
    for (i in 1:nsteps){
      alpha_mcmcsample[i+1] <- Metropolis(alpha_mcmcsample[i])
    }
    
    #burn_in
    alpha_mcmc_burned <- alpha_mcmcsample[-(1:burn_in)]
    length(alpha_mcmc_burned)
    
    #thin
    alpha_mcmc_thinned <- c()
    thin = 15
    for (i in 1:(nsteps-burn_in)){
      if(i%%thin == 0){
        alpha_mcmc_thinned <- append(alpha_mcmc_thinned, alpha_mcmc_burned[i])
      }
    }
    
    
    bayes_alpha_estimated <- mean(alpha_mcmc_thinned)
    bayes_beta_estimated <- sum(x_total*l_total)/(n-1)
    MLE_beta <- 1/mllomax(x_total)[1]
    MLE_alpha <- mllomax(x_total)[2]
    
    
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
  
  MLE_mre_beta_all <- append(MLE_mre_beta_all,(1/N)*sum(MLE_mres_beta))
  MLE_mse_beta_all <- append(MLE_mse_beta_all,(1/N)*sum(MLE_mses_beta))
}

pdf("result.pdf")
par(mfrow=c(2,2)) 
plot(sample_size, bayes_mre_alpha_all, ylim = c(0.95, 1.15), ylab = "", type="l", lwd = 2, col = "red")
par(new=T)
plot(sample_size, MLE_mre_alpha_all, ylim = c(0.95, 1.15), ylab = "MRE(alpha)",lty=6, type="l", lwd = 2, col = "blue")
abline(h=1,lty=2)

plot(sample_size, bayes_mse_alpha_all, ylim = c(0.0, 0.5), ylab = "", type="l", lwd = 2, col = "red")
par(new=T)
plot(sample_size, MLE_mse_alpha_all, ylim = c(0.0, 0.5), ylab = "MSE(alpha)", lty=6, type="l", lwd = 2, col = "blue")
abline(h=0,lty=2)
legend("topright", legend=c("MLE","Bayes"), lty=c(6,1),lwd = 2, col=c("blue","red"))

plot(sample_size, bayes_mre_beta_all, ylim = c(0.95, 1.2),  ylab = "", type="l", lwd = 2, col = "red") 
par(new=T)
plot(sample_size, MLE_mre_beta_all, ylim = c(0.95, 1.2),  ylab = "MRE(beta)", lty=6, type="l", lwd = 2, col = "blue")
abline(h=1,lty=2)

plot(sample_size, bayes_mse_beta_all, ylim = c(0.0, 2.5), ylab = "", type="l", lwd = 2, col = "red")
par(new=T)
plot(sample_size, MLE_mse_beta_all, ylim = c(0.0, 2.5),  ylab = "MSE(beta)", lty=6, type="l", lwd = 2, col = "blue")
abline(h=0,lty=2)
dev.off()
