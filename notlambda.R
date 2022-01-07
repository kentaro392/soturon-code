#2変数メトロポリス

library("univariateML")
library("extraDistr")

true_beta = 2
true_alpha = 1.5

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
  #1種類のサンプル数につき N回実験して平均を取る
  N = 1000
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
    
    #Lomaxに従う確率変数をn個生成
    for (k in 1:n){
      l <- rgamma(1,true_alpha,1)
      x <- rexp(1,l/true_beta)
      x_total <- append(x_total,x)
    }
    
    post <- #事後分布
      function(beta, alpha) 
        ((alpha^(n-1/2))*beta^(-n-1))/((alpha + 1)*(alpha + 2)^(1/2))*(prod(1+x_total/beta))^(-(alpha + 1))
    
    nsteps <- 31000 #mcmcサンプリングの数
    beta_mcmc <- rep(0,nsteps) 
    alpha_mcmc <- rep(0,nsteps)
    beta_mcmc[1] <- 5
    alpha_mcmc[1] <- 5
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
    
    bayes_beta_estimated <- mean(thinned_beta)
    bayes_alpha_estimated <- mean(thinned_alpha)	
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

pdf("notlambda.pdf")
par(mfrow=c(2,2)) 
plot(sample_size, bayes_mre_alpha_all, ylim = c(0.95, 1.2), ylab = "", type="l", lwd = 2, col = "red")
par(new=T)
plot(sample_size, MLE_mre_alpha_all, ylim = c(0.95, 1.2), ylab = "MRE(alpha)", lty=6, type="l", lwd = 2, col = "blue")
abline(h=1, lty=2)

plot(sample_size, bayes_mse_alpha_all,ylim = c(0, 1.5), ylab = "", type="l", lwd = 2, col = "red") 
par(new=T)
plot(sample_size, MLE_mse_alpha_all,ylim = c(0, 1.5), ylab = "MSE(alhpa)", lty=6, type="l", lwd = 2, col = "blue")
abline(h=0, lty=2)
legend("topright", legend=c("MLE","Bayes"), lty=c(6,1), lwd = 2, col=c("blue","red"))

plot(sample_size, bayes_mre_beta_all, ylim = c(0.95, 1.5), ylab = "", type="l", lwd = 2, col = "red") 
par(new=T)
plot(sample_size, MLE_mre_beta_all, ylim = c(0.95, 1.5), ylab = "MRE(beta)", lty=6, type="l", lwd = 2, col = "blue")
abline(h=1, lty=2)

plot(sample_size, bayes_mse_beta_all, ylim = c(0, 3.0), ylab = "", type="l", lwd = 2, col = "red")
par(new=T)
plot(sample_size, MLE_mse_beta_all, ylim = c(0, 3.0), ylab = "MSE(beta)", lty=6, type="l", lwd = 2, col = "blue")
abline(h=0, lty=2)
dev.off()

