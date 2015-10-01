library(FAdist)
library(fitdistrplus)
library(lmomco)
library(MASS)
library(moments)
library(rootSolve)
library(stats)
library(stats4)

#1 Simulando uma amostra da distribuição de Gumbel
# Gerando uma semente
set.seed(300)
# Definindo o tamanho da amostra
n = 30
# Definindo o parâmetro beta
sc = 1.0
# Definindo o parâmetro mi
loc = 5.0
# Simulando a distribuição
x.gum <- rgumbel(n=n,scale=sc,location=loc)
plot(ecdf(x.gum),main="",ylab="F(x)")
plot(density(x.gum),ylim=c(0,0.4),main="",ylab="f(x)",xlab="x")
curve(dgumbel(x,scale=sc,location=loc),n=n,from=2,to=10,lty=2,add=TRUE)
rug(x.gum)

#2 Método da Máxima Verossimilhança
fit.MV <- fitdist(x.gum, "gumbel", start=list(scale=5, location=5), method="mle")
summary(fit.MV)

# Parâmetros estimados pelo método
MVsc <- fit.MV[[1]][[1]]
MVloc <- fit.MV[[1]][[2]]
MV.gum <- curve(dgumbel(x,scale=MVsc,location=MVloc),n=n,from=2,to=10,col="firebrick1",add=TRUE)

#3 Método dos Momentos
mean = mean(x.gum)
s = sqrt(var(x.gum))
# Obtendo beta
MMsc = s*sqrt(6)/pi
# Cálculo de mi
MMloc = mean - 0.45006*s
MM.gum <- curve(dgumbel(x,scale=MMsc,location=MMloc),n=n,from=2,to=10,col="forestgreen",add=TRUE)

#4 Método dos Momentos-L
# Cálculo dos momentos-L amostrais
lmoments <- lmoms(x.gum)
# Estimando os parâmetros a partir dos momentos-L
para.lm <- lmom2par(lmoments, type="gum")
MLMsc <- para.lm[[2]][[2]]
MLMloc <- para.lm[[2]][[1]]
MLM.gum <- curve(dgumbel(x,scale=MLMsc,location=MLMloc),n=n,from=2,to=10,col="darkorchid",add=TRUE)
legend("topright",c("Simulada","Teórica","MMV","MM","MML"),col=c(1,1,"firebrick1","forestgreen","darkorchid"),lty=c(1,2,1,1,1))