library(bbmle)
library(evd)
library(evir)
library(FAdist)
library(lmomco)
library(MASS)
library(nsRFA)
library(rootSolve)
library(stats)
library(stats4)

#1 Simulando uma amostra da distribuição de Pareto
# Gerando uma semente
set.seed(3000)
# Definindo o tamanho da amostra
n = 200
# Definindo o parâmetro csi
sh = 1.0
# Definindo o parâmetro mi
loc = 0.0
# Definindo o parâmetro sigma
sc = 10.0
# Simulando a distribuição
x.gpa <- rgpd(n,sc=sc, loc=loc, shape=sh)
plot(ecdf(x.gpa),main="",ylab="F(x)")
plot(density(x.gpa),ylim=c(0,0.1),main="",ylab="f(x)",xlab="x")
curve(dgpd(x,sc=sc,loc=loc,shape=sh),n=n,from=0,to=120,lty=2,add=TRUE)
rug(x.gpa)

#2 Método da Máxima Verossimilhança
# Função de Log-verossimilhança
logl <- function(sigma,csi,dados,n) {
  nll <- +n*log(sigma)+(1-csi)*sum((-1/csi)*log(1-csi*dados/sigma))
}

# Otimizador
est <- mle(minuslogl=logl,start=list(sigma=5,csi=10),fixed=list(dados=x.gpa, n=length(x.gpa)), method="BFGS")
summary(est)

# Parâmetros estimados pelo método
MVsc <- getElement(coef(est),1)
MVsh <- getElement(coef(est),2)
MV.wei <- curve(dgpd(x,scale=MVsc,shape=MVsh),n=n,from=0,to=120,col="firebrick1",add=TRUE)

#3 Método dos Momentos
mean = mean(x.gpa)
sd = sd(x.gpa)
# Equação para se obter csi
MMsh <- 0.5*mean*(mean^2/(sd+1))
# Equação para se obter sigma
MMsc <- 0.5*(mean^2/(sd-1))
MM.gpa <- curve(dgpd(x,scale=MMsc,shape=MMsh,loc=0.0),n=n,from=0,to=120,col="forestgreen",add=TRUE)

#4 Método dos Momentos-L
# Cálculo dos momentos-L amostrais
lmoments <- lmoms(x.gpa)
# Estimando os parâmetros a partir dos momentos-L
para.lm <- lmom2par(lmoments, type="gpa")
MLMsh <- para.lm[[2]][[1]]
MLMsc <- para.lm[[2]][[2]]
MLMth <- para.lm[[2]][[3]]
MLM.gpa <- curve(dgpd(x,scale=MLMsc,shape=-MLMsh,loc=-MLMth),n=n,from=0,to=1000,col="darkorchid",add=TRUE)
legend("topright",c("Simulada","Teórica","MMV","MM","MML"),col=c(1,1,"firebrick1","forestgreen","darkorchid"),lty=c(1,2,1,1,1))