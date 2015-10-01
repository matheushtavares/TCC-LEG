library(lmomco)
library(MASS)
library(stats)
library(stats4)
library(FAdist)

#1 Simulando uma função Log-normal
# Gerando uma semente
set.seed(4000)
# Definindo o tamanho da amostra
n = 20
# Definindo a log-média
m = 1
# Definindo o log-desvio padrão
sd = 0.5
# Simulando a distribuição
x.ln <- rlnorm(n=n,meanlog=m,sdlog=sd)
plot(ecdf(x.ln),main="Distribuição Log-Normal Acumulada",ylab="F(x)")
y.ln <- density(x.ln,n=n)
plot(density(x.ln),ylim=c(0,0.5),main="",ylab="f(x)")
curve(rlnorm(x,meanlog=m,sdlog=sd),n=n,from=0,to=10,lty=2,add=TRUE)
rug(x.ln)
mean <- mean(x.ln)
var <- var(x.ln)

#2 Método da Máxima Verossimilhança
MVest <- fitdistr(x.ln,"lognormal")
# Parâmetro obtidos pelo método
MVm <- MVest[[1]][[1]]
MVsd <- MVest[[1]][[2]]
MV.ln <- curve(dlnorm(x,meanlog=MVm,sdlog=MVsd),n=n,from=0,to=8,col="firebrick1",add=TRUE)

#3 Método dos Momentos
est.m <- function(dados,n) {
  mean <- 2*log((1/n)*sum(dados))-0.5*log((1/n)*sum(dados^2))
  return(mean)
}

est.var <- function(dados,n) {
  var <- log((1/n)*sum(dados^2))-2*log((1/n)*sum(dados))
  return(var)
}
# Obtendo os parâmetros pelo método
MMm = est.m(x.ln,n)
MMsd = sqrt(est.var(x.ln,n))
MM.ln <- curve(dlnorm(x,meanlog=MMm,sdlog=MMsd),n=n,from=0,to=7,col="forestgreen",add=TRUE)

#4 Método dos Momentos-L
# Cálculo dos momentos-L amostrais
lmoments <- lmoms(x.ln)
# Estimando os parâmetros a partir dos momentos-L
para.lm <- lmom2par(lmoments, type="ln3")
MLMsd <- para.lm[[2]][[3]]
MLMm <- para.lm[[2]][[2]]
MLMth <- para.lm[[2]][[1]]
MLM.ln <- curve(dlnorm3(x,scale=MLMm,shape=MLMsd,thres=MLMth),n=n,from=0,to=8,col="darkorange",add=TRUE)
legend("topright",c("Simulada","Teórica","MMV","MM","MML"),col=c(1,1,"firebrick1","forestgreen","darkorchid"),lty=c(1,2,1,1,1))