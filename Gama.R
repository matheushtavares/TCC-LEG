library(lmomco)
library(MASS)
library(stats)
library(stats4)

#1 Simulando uma amostra da distribuição Gama
# Gerando uma semente
set.seed(3000)
# Definindo o tamanho da amostra
n = 80
# Definindo o par?metro beta
rt = 0.75
# Definindo o par?metro k
sh = 3.0
# Simulando a distribuição
x.gam <- rgamma(n=n,rate=rt,shape=sh)
plot(ecdf(x.gam),main="",ylab="F(x)")
plot(density(x.gam),ylim=c(0,0.25),main="",ylab="f(x)",xlab="x")
curve(dgamma(x,rate=rt,shape=sh),n=n,from=0,to=15,lty=2,add=TRUE)
rug(x.gam)
# Cálculo da média e variância amostrais
mean <- mean(x.gam)
var <- var(x.gam)

#2 Método da Máxima Verossimilhança
# Função de Log-verossimilhança
logl <- function(beta,k,dados,n) {
  nll <- -n*k*log(beta)+n*log(gamma(k))-(k-1)*sum(log(dados))+beta*sum(dados)
  return(nll)
}
#Otimizador
est <- mle(minuslogl=logl,start=list(beta=2,k=1),fixed=list(dados=x.gam, n=length(x.gam)), method="L-BFGS-B", lower=c(0.001,0.001), upper=c(Inf,Inf))
summary(est)

# Parâmetros estimados pelo método
MVrt <- getElement(coef(est),1)
MVsh <- getElement(coef(est),2)
MV.gam <- curve(dgamma(x,rate=MVrt,shape=MVsh),n=n,from=0,to=20,col="firebrick1",add=TRUE)

#3 Método dos Momentos
# Cálculo de beta
est.sh <- function(dados,mean,n) {
  esh <- (mean^2)/((1/n)*sum((dados^2)-(mean^2)))
  return(esh)
}
MMsh = est.sh(x.gam,mean,n)
# Cálculo de k
est.rt <- function(dados,mean,n) {
  ert <- (mean)/((1/n)*sum((dados^2)-(mean^2)))
  return(ert)
}
MMrt = est.rt(x.gam,mean,n)
MM.gam <- curve(dgamma(x,rate=MMrt,shape=MMsh),n=n,from=0,to=20,col="forestgreen",add=TRUE)

#4 Método dos Momentos-L
# Cálculo dos momentos-L amostrais
lmoments <- lmoms(x.gam)
# Estimando os parâmetros a partir dos momentos-L
para.lm <- lmom2par(lmoments, type="gam")
MLMsh <- para.lm[[2]][[1]]
MLMrt <- 1/para.lm[[2]][[2]]
MLM.gam <- curve(dgamma(x,rate=MLMrt,shape=MLMsh),n=n,from=0,to=20,col="darkorchid",add=TRUE)
legend("topright",c("Simulada","Te?rica","MMV","MM","MML"),col=c(1,1,"firebrick1","forestgreen","darkorchid"),lty=c(1,2,1,1,1))