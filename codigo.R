#teste

library(TeachingDemos)
library(base64)
library(car)
library(MASS)
library(clusterGeneration)
library(mvtnorm)


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




endereco = "/home/cicconella/Analise-Multivariada/MAE0330-Caes"

caes = read.table(endereco, header = T, sep = ";")
caes
dim(caes)

lucas = c(6.8,15.5,15.5, 6, 23.9,29,1)
liting = c(14,25.5,23.8,11.5,48.3,41.5)
#caes[2,4] = 7.7
victor =  c(9,16,22,8,35,42)
adriana = c(6.57,18.74,19.12,4.28,30.47,35.87)
manuela = c(10,20,16,8,28,32)
carlos = c(11.19,3.22,6.9,5.36,4.40, 8)
lai = c(13,23,21,8,32,35)
#caes[4,] = caes[3,]
debora = c(12.44,24.79, 25.70,9.9, 37.5, 43.95)
rodrigo = c(9.8, 20.6,19.7,8.2,39.1,36)
vini = c(10.9,24.1,23,9,36,44.1)
alex = c(11,27,22,9,33,40)

obs = alex

caes = rbind(caes, obs)
caes

media = colMeans(caes)
print(media)

var = round(cov(caes),digits = 2)
print(var)

cor = round(cor(caes),2)
print(cor)

vartotal <- sum(diag(var))
vartotal
vargen <- det(var)
vargen

caesEstrela<-scale(caes,center=TRUE,scale=TRUE)
caesEstrela

mediaE = colMeans(caesEstrela)
round(mediaE,2)

varE = round(var(caesEstrela),2)
varE

corE = round(cor(caesEstrela),2)
corE

vartotalE <- sum(diag(varE))
vartotalE
vargenE <- round(det(varE), 2)
vargenE

#require("knitr")

#?kable
#scatterplotMatrix(caes, reg.line=FALSE, smooth=FALSE, 
#                  spread=FALSE, span=0.5, diagonal = 'boxplot', data=y)
#scatterplotMatrix(caes, reg.line=lm, smooth=FALSE, 
#                  spread=FALSE, span=0.5, diagonal = 'boxplot', data=y)


euclidiana = dist(caes)
euclidiana

bivbox(caesEstrela[,c(1,4)], method = "O")
bivbox(caes[,c(1,2)], method = "O")

maha = mahalanobis(caes, media, var)
maha

plot(maha,ylim = c(1,13))
abline(qchisq(.90, df=6), 0)



####Euclidiana

euclidiana = dist(caes)

euclidiana = as.matrix(euclidiana)
euclidiana

#Medias das distancias
print(mean(euclidiana[7,-7]))

modernos = caes[-7,]
ancestral = caes[7,]

modernos
ancestral

media = colMeans(modernos)
centroide = matrix(media, ncol=6)

#Distancia do centroide
x = rbind(media,ancestral)
dist(x)

#Dist do minimo
min(euclidiana[7,-7])
which(euclidiana[7,-7]==min(euclidiana[7,-7]))

#Dist do maximo
max(euclidiana[7,-7])
which(euclidiana[7,-7]==max(euclidiana[7,-7]))


####Mahalanobis

caes

cov = round(var(caes[-7,]),2)
vetor = as.numeric(caes[7,])
maha = mahalanobis(modernos, vetor, cov)
maha = sqrt(maha)

#Media das medias 
mean(maha)

#Distancia dos centroides
maha2 = mahalanobis(centroide, vetor, cov)
sqrt(maha2)










