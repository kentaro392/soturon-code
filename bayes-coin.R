#コインを10回投げた結果
data <- c(0,1,1,1,1,0,1,1,0,1)

#尤度
LL_Bern <- function(x,q) {
  q^sum(x)*(1-q)^(length(x)-sum(x))
}
plot(seq(0,1,0.01),LL_Bern(data,seq(0,1,0.01)),type = "l")


#対数尤度
logLL_Bern <- function(x,q) {
  sum(x)*log(q)+(length(x)-sum(x))*log(1-q)
}
plot(seq(0,1,0.01),logLL_Bern(data, seq(0,1,0.01)),type = "l")

#最大化問題を解くとMLEが得られる
optimize(function(q) LL_Bern(data,q),c(0,1),maximum = TRUE) 

optimize(function(q) logLL_Bern(data,q),c(0,1),maximum = TRUE) 


#事後分布の推定
#ここでは事前分布はBeta(1,1)とする

prior_beta <- function(q,a,b) dbeta(q,a,b) 
joint <- function(x,q) LL_Bern(x,q)*prior_beta(q,1,1)
plot(seq(0,1,0.01),joint(data,seq(0,1,0.01)),type = "l")

#メトロポリスアルゴリズム
#一様乱数との比較で採択or棄却を決める
Metropolis <- function(current) {
  propose <- runif(1,0,1) #提案分布として一様分布を採用
  #propose <- currentを中心とした正規分布とすれば0,1以外の時でも対応
  postOdds <- joint(data, propose)/joint(data, current)
  pmove <- min(postOdds,1)
  if(pmove >= runif(1,0,1)) propose else current
}

nsteps <- 11000 #mcmcサンプリングの数
mcmcsample <- rep(NA,nsteps + 1) #要素が欠損している11000次元ベクトルを生成
mcmcsample[1] <- 0.5

for (i in 1:nsteps){
	mcmcsample[i+1] <- Metropolis(mcmcsample[i])
}

plot(mcmcsample[-(1:1000)],type="l",col="blue")

class(mcmcsample)


#プロット
hist(mcmcsample[-(1:1000)],freq = FALSE ,col="blue")
lines(seq(0,1,0.01),dbeta(seq(0,1,0.01),1+sum(data),1+length(data)-sum(data)),col="red")