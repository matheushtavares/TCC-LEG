library(FAdist)
library(lmomco)
library(MASS)
library(moments)
library(PearsonDS)
library(rootSolve)
library(stats)
library(stats4)

#1 Simulando uma amostra da distribuição de Pareto
# Gerando uma semente
set.seed(300)
# Definindo o tamanho da amostra
n = 50
# Definindo o parâmetro a
sh = 3.0
# Definindo o parâmetro lambda
thr = 0.0
# Definindo o parâmetro s
sc = 5.0
# Simulando a distribuição
x.pt3 <- rgamma3(n=n,shape=sh, thres=thr, scale=sc)
plot(ecdf(x.pt3),main="",ylab="F(x)")
plot(density(x.pt3),ylim=c(0,0.2),main="",ylab="f(x)",xlab="x")
curve(dgamma3(x,scale=sc,shape=sh,thres=thr),n=n,from=0,to=50,lty=2,add=TRUE)
rug(x.pt3)

#2 Método da Máxima Verossimilhança
# Função de Log-verossimilhança
logl <- function(sigma,mi,csi,dados,n) {
  nll <- +n*log(mi)-((1-csi)/csi)*sum(log(1-(csi/sigma)*(dados-mi)))
}

# Otimizador
est <- mle(minuslogl=logl,start=list(sigma=5,mi=2,csi=10),fixed=list(dados=x.gpa, n=length(x.gpa)), method="BFGS")
summary(est)

# Parâmetros estimados pelo método
MVsc <- getElement(coef(est),1)
MVsh <- getElement(coef(est),2)
MVth <- getElement(coef(est),3)
MV.wei <- curve(dweibull3(x,scale=MVsc,shape=MVsh,thres=MVth),n=n,from=0,to=1.6,col="firebrick1",add=TRUE)

#3 Método dos Momentos
l.pt3 <- log(x.pt3)
mean = mean(l.pt3)
sd = sd(l.pt3)
skew=skewness(l.pt3)
# Cálculo dos 3 primeiros momentos amostrais
momentos <- all.moments(l.pt3,order.max=3)
p_mom <- momentos[1]
s_mom <- momentos[2]
t_mom <- momentos[3]
# Equação para se obter a
achaa <- function(a) {
  ((log(t_mom)-3*log(p_mom))*(2*log(1-a)-log(1-2*a)))-((log(s_mom)-2*log(p_mom))*(3*log(1-a)-log(1-3*a)))
}
# Encontrando as raízes da equação
MMsh <- uniroot.all(achaa,c(0.1,1000))
csi <- MMsh
# Equação para se obter mi
achami <- function(mi) {
  sd-((mi^2)/((1+2*csi)*(1+csi)^2))
}
# Encontrando as raízes da equação
MMloc <- uniroot.all(achami,c(0.1,1000))
mi <- MMloc
# Cálculo de sigma
achasigma <- function(sigma) {
  mean-((sigma+mi)/(1+csi))
}
# Encontrando as raízes da equação
MM.pt3 <- curve(dgamma3(x,scale=MMsc,shape=MMsh,thres=MMth),n=n,from=0,to=50,col="forestgreen",add=TRUE)

#4 Método dos Momentos-L
# Cálculo dos momentos-L amostrais
lmoments <- lmoms(x.pt3)
# Estimando os parâmetros a partir dos momentos-L
para.lm <- lmom2par(lmoments, type="pe3")
MLMsh <- para.lm[[2]][[2]]
MLMsc <- para.lm[[2]][[3]]
MLMth <- para.lm[[2]][[1]]
MLM.pt3 <- curve(dgamma3(x,scale=MLMsc,shape=MLMsh,thres=-MLMth),n=n,from=0,to=50,col="darkorchid",add=TRUE)
legend("topright",c("Simulada","Teórica","MMV","MM","MML"),col=c(1,1,"firebrick1","forestgreen","darkorchid"),lty=c(1,2,1,1,1))