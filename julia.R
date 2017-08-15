library(TeachingDemos)
library(base64)
library(car)
library(MASS)
library(clusterGeneration)
library(mvtnorm)


## Exercício 1: Entrando com uma matriz de dados
## Idade (anos) e preço de carros usados
## Ex.1.2 pag.28, Johnson and Wichern,1992
y <- matrix(c(3,5,5,7,7,7,8,9,10,11,2.30,1.9,1,0.7,0.3,1,1.05,0.45,0.7,0.3),ncol=2,nrow=10)
y
str(y)
#### Algumas estatísticas multivariadas
my<-colMeans(y)
my 
vy<-var(y)
vy
cvy<-cov(y)
cvy
cvy<-round(cov(y),2)
cvy
roy<-round(cor(y),2)
roy
summary(y)

##Medidas de variabilidade multivariada
vartotal <- sum(diag(vy))
vartotal
vargen <- det(vy)
vargen

##Propriedades da matriz de covariâncias
## Cov(Ya)=a´Cov(Y)a   Cov(YA)A´Cov(Y)A
## Defina um vetor a (combinação linear das variáveis)
## Defina uma matriz A 
## Rotação dos eixos: rotação ortogonal AA’=I)
## Ah<-matrix(c(cosTeta, senTeta, -senTeta, cosTeta),2,2):horário
## Aah<-matrix(c(cosTeta, -senTeta, senTeta, cosTeta),2,2):anti-horário

##Transformação de variáveis
## Normalização das variáveis
## Compare as matrizes covar e cor dos dados originais e padronizados
yz<-scale(y,center=TRUE,scale=TRUE)
yz
colMeans(yz)
var(yz) #revise propriedades
cor(yz) #revise propriedades 

###Representação gráfica
plot(y) #pense na covariância como uma forma quadrática
scatterplotMatrix(~y[,1]+y[,2], reg.line=FALSE, smooth=FALSE, 
                  spread=FALSE, span=0.5, diagonal = 'boxplot', data=y)
scatterplotMatrix(~y[,1]+y[,2], reg.line=lm, smooth=FALSE, 
                  spread=FALSE, span=0.5, diagonal = 'boxplot', data=y)
##scatterplotMatrix(~y[,1]+y[,2]+y[,3]+y[,4]+y[,5]+y[,6], reg.line=FALSE, smooth=FALSE, spread=FALSE, span=0.5, diagonal = 'boxplot', data=y)

library(ggplot2)
dt<-as.data.frame(caes)
ggplot(dt, aes(x=dt[,1],y=dt[,2],label=rownames(dt))) + geom_text()

#### Faces de Chernoff – Discuta os resultados
library(TeachingDemos)
faces( as.matrix(y), fill=T) ##unidades amostrais
faces( as.matrix(t(y)), fill=T) ##variáveis: pode NÃO fazer sentido!

####Pesquise como construir o gráfico radar no R

#### Gerando dados da Normal Multivariada (Bivariada)
library(MASS)
##set.seed(1298)
mu<-c(0,0)
sigma<-matrix(c(2,1,1,2),ncol=2)
n<-10
yr<-mvrnorm(n,mu,sigma)
yr
mi<-colMeans(yr)
mi
s<-cov(yr)
s
plot(yr)

##Alternativas para gerar dados da normal multivariada
library(clusterGeneration)
library(mvtnorm)
## Definir “n” e “p”
mu<-rnorm(p)
R<-rcorrmatrix(p,alphad=1) 
yr<-rmvnorm(n,mu,R)
yr


y=caes
#### Distância Euclidiana (ordinária) ENTRE OBSERVAÇÕES
de<-dist(y)
de
de37<-sqrt(t(y[3,]-y[7,])%*%(y[3,]-y[7,])) ##raiz quadrada de dist 
de37
round(de,2)
dim(de)
de <- as.matrix(de)
dim(de)
de37<-de[3,7]
de37


#### Distância Euclidiana Ponderada (dist de Pearson) ENTRE OBSERVAÇÕES
sdy<-cbind(sd(y[,1]),sd(y[,2]))
sdy
ysd<-sweep(y,2,sdy,FUN="/") #alternativa 1 de ponderação
ysd
yzs<-scale(y, center = FALSE, scale = apply(y, 2, sd, na.rm = TRUE)) #alternativa 2 de ponderação
yzs ##ysd=yzs
dep<-dist(ysd)
dep
round(dep,2)
dep<-as.matrix(dep)
dep37<-dep[3,7]
dep37  #compare as distâncias Euclidiana e de Pearson

### Distância de Mahalanobis das observações AO CENTRÓIDE 
dm<-mahalanobis(y,my,cvy)
dm  #dist de Mahalanobis ao quadrado
my<-as.matrix(my)
dm4my<-t(y[4,]-my)%*%solve(cvy)%*%(y[4,]-my) ##dist ao quadrado
dm4my
dm3my<-t(y[3,]-my)%*%solve(cvy)%*%(y[3,]-my)
dm3my

de4my<-t(y[4,]-my)%*%(y[4,]-my)
de4my  #dist Euclidiana da obs 4 ao centróide
de3my<-t(y[3,]-my)%*%(y[3,]-my)
de3my

##Calcule a Distância de Mahalanobis das observações padronizadas ao centróide

### Distância de Mahalanobis ENTRE OBSERVAÇÕES 
###A obs 7 está mais próxima da obs 3 ou da obs 1?
dm37<-t(y[3,]-y[7,])%*%solve(cvy)%*%(y[3,]-y[7,])
dm37
dm17<-t(y[1,]-y[7,])%*%solve(cvy)%*%(y[1,]-y[7,])
dm17

de[3,7]
de[1,7]

dep[3,7]
dep[1,7]

plot(y)

##Calcule a matriz de Distâncias de Mahalanobis ENTRE todas as observações

# Função bivbox do livro de Everitt
######################
bivbox<-function(a, d = 7, mtitle = "Bivariate Boxplot",
                 method = "robust",xlab="X",ylab="Y")
{
  #a is data matrix
  #d is constant(usually 7 para 99% de confiança)
  p <- length(a[1,  ])
  if(method == "robust") {
    param <- biweight(a[, 1:2]); m1 <- param[1]; m2 <- param[2]
    s1 <- param[3]; s2 <- param[4]; r <- param[5]
  }
  else {
    m1 <- mean(a[, 1]); m2 <- mean(a[, 2]); 
    s1 <- sqrt(var(a[, 1])); s2 <- sqrt(var(a[, 2])); r <- cor(a[, 1:2])[1, 2]
  }
  x <- (a[, 1] - m1)/s1; y <- (a[, 2] - m2)/s2
  e <- sqrt((x * x + y * y - 2 * r * x * y)/(1 - r * r))
  e2 <- e * e; em <- median(e); emax <- max(e[e2 < d * em * em])
  r1 <- em * sqrt((1 + r)/2); r2 <- em * sqrt((1 - r)/2); theta <- ((2 * pi)/360) * seq(0, 360, 3)
  xp <- m1 + (r1 * cos(theta) + r2 * sin(theta)) * s1; yp <- m2 + (r1 * cos(theta) - r2 * sin(theta)) * s2
  r1 <- emax * sqrt((1 + r)/2); r2 <- emax * sqrt((1 - r)/2); theta <- ((2 * pi)/360) * seq(0, 360, 3)
  xpp <- m1 + (r1 * cos(theta) + r2 * sin(theta)) * s1; ypp <- m2 + (r1 * cos(theta) - r2 * sin(theta)) * s2
  maxxl <- max(xpp); minxl <- min(xpp); maxyl <- max(ypp); minyl <- min(ypp)
  b1 <- (r * s2)/s1; a1 <- m2 - b1 * m1; y1 <- a1 + b1 * minxl; y2 <- a1 + b1 * maxxl
  b2 <- (r * s1)/s2; a2 <- m1 - b2 * m2; x1 <- a2 + b2 * minyl; x2 <- a2 + b2 * maxyl
  maxx <- max(c(a[, 1], xp, xpp, x1, x2)); minx <- min(c(a[, 1], xp, xpp, x1, x2))
  maxy <- max(c(a[, 2], yp, ypp, y1, y2)); miny <- min(c(a[, 2], yp, ypp, y1, y2))
  plot(a[, 1], a[, 2], xlim = c(minx, maxx), ylim = c(miny, maxy), xlab =xlab, ylab =ylab,
       lwd = 2, pch = 1)
  lines(xp, yp, lwd = 2); lines(xpp, ypp, lty = 2, lwd = 2)
  segments(minxl, y1, maxxl, y2, lty = 3, lwd = 2); segments(x1, minyl, x2, maxyl, lty = 4, lwd = 2)
}
########################

bivbox(y, method ="O")


## Exercício 2: Refaça o roteiro de aula “para os dados dos cães - Manly (2004)”
#Responda as questões propostas na Oficina CEC1
caes<-read.table("F:/MAE0330/2016/caes.csv", header = TRUE, sep=';', dec='.')
caes
