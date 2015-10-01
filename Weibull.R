library(FAdist)
library(lmomco)
library(MASS)
library(moments)
library(rootSolve)
library(stats)
library(stats4)

#1 Simulando uma amostra da distribuição de Weibull
#Gerando uma semente
set.seed(3000)
# Definindo o tamanho da amostra
n = 50
# Definindo o par?metro beta
sc = 1.0
# Definindo o par?metro k
sh = 5.0
# Definindo o par?metro csi
thres = 0.0
# Simulando a distribuição
x.wei <- rweibull3(n=n,scale=sc,shape=sh,thres=thres)
plot(ecdf(x.wei),main="",ylab="F(x)")
plot(density(x.wei),ylim=c(0,2.5),main="",ylab="f(x)",xlab="x")
curve(dweibull3(x,scale=sc,shape=sh,thres=thres),n=n,from=0,to=1.6,lty=2,add=TRUE)
rug(x.wei)
# Cálculo da média e variância amostrais
mean <- mean(x.wei)
var <- var(x.wei)
sd <- sd(x.wei)

#2 Método da Máxima Verossimilhança
# Função de Log-verossimilhança
logl <- function(lambda,alpha,csi,dados,n) {
  nll <- -n*(log(alpha)-alpha*log(lambda))-(alpha-1)*sum(log(dados-csi))+(1/(lambda^alpha))*sum((dados-csi)^alpha)
  return(nll)
}

#Otimizador
est <- mle(minuslogl=logl,start=list(lambda=2,alpha=1,csi=0),fixed=list(dados=x.wei, n=length(x.wei)), method="BFGS")
summary(est)

# Parâmetros estimados pelo método
MVsc <- getElement(coef(est),1)
MVsh <- getElement(coef(est),2)
MVth <- getElement(coef(est),3)
MV.wei <- curve(dweibull3(x,scale=MVsc,shape=MVsh,thres=MVth),n=n,from=0,to=1.6,col="firebrick1",add=TRUE)

#3 Método dos Momentos
cv <- sd/mean
skew=skewness(x.wei)
# Cálculo dos 3 primeiros momentos amostrais
momentos <- all.moments(x.wei,order.max=3)
p_mom <- momentos[1]
s_mom <- momentos[2]
t_mom <- momentos[3]
min_shape = 0.1
max_shape = 1000
# Equação para se obter k
achak <- function(k) {
  (((g3 - 3*g1*g2 + 2*g1)/(g2 - g1^2)^1.5)/skew)-1
}
achak <- function(k) {
  ((sqrt(gamma(1+(2/k))-(gamma(1+(1/k)))^2)/gamma(1+(1/k)))/cv)-1
}
# Encontrando as raízes da equação
MMsh <- uniroot.all(achak,c(min_shape,max_shape))
# Cálculo de beta
MMsc = (mean/gamma(1+(1/MMsh)))^MMsh
# Cálculo de csi
MMth = p_mom-MMsc*gamma(1+(1/MMsh))
MM.wei <- curve(dweibull3(x,scale=MMsc,shape=MMsh,thres=MMth),n=n,from=0,to=1.6,col="forestgreen",add=TRUE)

#4 Método dos Momentos-L
# Calculo dos momentos-L amostrais
lmoments <- lmoms(x.wei)
# Estimando os parâmetros a partir dos momentos-L
para.lm <- lmom2par(lmoments, type="wei")
MLMsh <- para.lm[[2]][[3]]
MLMsc <- para.lm[[2]][[2]]
MLMth <- para.lm[[2]][[1]]
MLM.wei <- curve(dweibull3(x,scale=MLMsc,shape=MLMsh,thres=-MLMth),n=n,from=0,to=1.6,col="darkorchid",add=TRUE)
legend("topright",c("Simulada","Te?rica","MMV","MM","MML"),col=c(1,1,"firebrick1","forestgreen","darkorchid"),lty=c(1,2,1,1,1))