setwd(dir = "Desktop/Code/R/")
source(file = "Multiplot.R")                         # From 'cookbook for R'
source(file = "central_posterior_interval_update.R") # Mine, hosted on Github
library("ggplot2")

#First plot--------------------------------------------
#Example beta densities.
theta <- seq(from = 0, to= 1, by = 0.001)

betaDens <- function(theta, n, x) {              # Of course, can just use dbeta(x+1, n-x+1) instead
  return((theta)^x*(1-theta)^(n-x))
}

Dat <- data.frame(theta, betaDens(theta, 5, 1), betaDens(theta, 10, 4), betaDens(theta, 30, 18), betaDens(theta, 1000, 500))
names(Dat) <- c("Theta", "Density1", "Density2", "Density3", "Density4")

# n=5; x=1
p1 <- ggplot(data = Dat, aes(x=Theta, y=Density1)) +
  geom_line() +
  labs(title="n=5, x=1", y="Density", x=NULL) +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
# n=10; x=4
p2 <- ggplot(data = Dat, aes(x=Theta, y=Density2)) +
  geom_line() +
  labs(title="n=10, x=4", y=NULL, x=NULL) +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
# n=30; x=18
p3 <- ggplot(data = Dat, aes(x=Theta, y=Density3)) +
  geom_line() +
  labs(title="n=30, x=18", y="Density") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
# n=1000; x=500
p4 <- ggplot(data = Dat, aes(x=Theta, y=Density4)) +
  geom_line() +
  labs(title="n=1000, x=500", y=NULL) +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())

multiplot(p1, p3, p2, p4, cols=2)
#end First plot-----------------------------------------

#Normalize density of n=30; x=18
densityNormal <- Dat$Density3/(sum(Dat$Density3*.001))
#Find HDI
interval.around.mode1D(densityNormal, p = 0.95)

#Second plot--------------------------------------------
x <- 0:10
theta <- seq(from = 0, to = 1, by = 0.001)
alpha <- seq(from = -0.5, to = 10.5, by = 0.01)
beta <- seq(from = -0.5, to = 10.5, by = 0.01)

densAlpha <- c(rep(0, times=50), 2/10-2/100*alpha[51:1051], rep(0, times=50))
densBeta <- c(rep(0, times=50), 2/10-2/100*beta[51:1051], rep(0, times=50))
densTheta <- dbeta(theta, shape1 = 3, shape2 = 4)
densBinom <- dbinom(x = x, size = 10, prob = 0.3)

Dat2 <- cbind.data.frame(alpha, beta, densAlpha, densBeta)
names(Dat2) <- c("Alpha", "Beta", "densAlpha", "densBeta")
p5 <- ggplot(Dat2, aes(x=Alpha, y=densAlpha)) +
  geom_line() +
  labs(title="Flat, Independent Hyperpriors for Alpha and Beta", y="Density") +
  xlim(c(-0.5, 10.5)) +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
p6 <- ggplot(Dat2, aes(x=Beta, y=densBeta)) +
  geom_line() +
  labs(title=" ", y="Density") +
  xlim(c(-0.5, 10.5)) +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())

multiplot(p5, p6)
#end Second plot---------------------------------------


#Third plot--------------------------------------------
qplot(theta, densTheta, geom="line", xlab="Theta", ylab="Density", main="Beta(3,4) Distribution")
#end Third plot----------------------------------------



#Fourth plot--------------------------------------------
Dat3 <- data.frame(x, densBinom)
p7 <- ggplot(data = Dat3, aes(x=x, y=densBinom)) +
  geom_segment(xend=x, yend=0) +
  scale_x_continuous(breaks=seq(0, 10, by=2)) +
  labs(x="X", y="Mass", title="Binomial PMF with parameters (10, 0.3)")
p7
#end Fourth plot----------------------------------------

# Get sequence of alpha, beta values:
a <- seq(from = 0.001, to = 10, length.out=1000)
b <- seq(from = 0.001, to = 10, length.out=1000)

#Hierarchical Model posterior for Phi (derived on slide):
computePhiDensityMatrix <- function(x, y) {
  PhiDensityMatrix <- as.matrix(rep(0, times = length(y)))
  for (i in 1:length(x)) {
    phiDensity <- (2/10-2/100*x[i])*(2/10-2/100*y)*(beta(x[i]+5, y+5))*(beta(x[i]+3, y))*(beta(x[i]+5, y+2))*(beta(x[i]+7, y+4))/(beta(x[i], y))^4
                  # beta(...) is the beta function, not a variable
    PhiDensityMatrix <- cbind(PhiDensityMatrix, phiDensity)
  }
  return(PhiDensityMatrix[1:dim(PhiDensityMatrix)[1],2:dim(PhiDensityMatrix)[2]])
}

# P(phi)
Matrix <- computePhiDensityMatrix(x = a, y = b)  # 10-15 seconds to compute
Matrix <- Matrix/sum(Matrix)

#Fifth plot ---------------------------------------------
contour(x=a, y=b, t(Matrix)) # t(Matrix) because contour() reads "x" as rows and "y" as columns, possibly counterintuitively
title(main = "Contour plot for density of Phi", xlab = "Alpha", ylab = "Beta")
#end Fifth plot -----------------------------------------

# Vector of samples, first define Rejection Sample function
rejectionSample <- function(dataMatrix, n) {
  vectorProb <- rep(NA, times = n)
  vectorPar1 <- rep(NA, times = n)
  vectorPar2 <- rep(NA, times = n)
  count <- 0
  numTries <- 0
  while (count < n) {
    x <- sample(seq(1:dim(dataMatrix)[2]), size = 1)
    y <- sample(seq(1:dim(dataMatrix)[1]), size = 1)
    ran <- runif(1)
    Max <- max(dataMatrix)
    numTries <- numTries + 1
    if (dataMatrix[y,x]/Max >= ran) {
      count <- count + 1
      vectorProb[count] <- dataMatrix[y,x]
      vectorPar1[count] <- x/100
      vectorPar2[count] <- y/100
  } }
  ranAlpha <<- vectorPar1[!is.na(vectorPar1)]
  ranBeta <<- vectorPar2[!is.na(vectorPar2)]
  ranSample <<- vectorProb[!is.na(vectorProb)]
  print(c("Acceptance rate:", n/numTries), quote=F)
}

rejectionSample(Matrix, n=1000) # Takes ~15 seconds to wait for 1000 samples, good enough. ~16% acceptance rate.

#Sixth plot----------------------------------------------
Dat3 <- data.frame(ranAlpha, ranBeta)
p8 <- ggplot(data=Dat3, aes(ranAlpha)) +
  geom_histogram(fill="grey", col="black", bins=15) +
  labs(x="Alpha", y=NULL, title="Sampled Marginal of Alpha") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
p9 <- ggplot(data=Dat3, aes(ranBeta)) +
  geom_histogram(fill="grey", col="black", bins=15) +
  labs(x="Beta", y=NULL, title="Sampled Marginal of Beta") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
multiplot(p8, p9, cols = 2)
#end sixth plot------------------------------------------

#P(theta | phi)
#Theta1
theta1 <- rep(NA, times = length(ranAlpha))
for (i in 1:length(ranAlpha)) {
  theta1[i] <- rbeta(n=1, ranAlpha[i]+5, ranBeta[i]+5)
}

#theta4
theta4 <- rep(NA, times = length(ranAlpha))
for (i in 1:length(ranAlpha)) {
  theta4[i] <- rbeta(n=1, ranAlpha[i]+7, ranBeta[i]+4)
}

#Seventh plot------------------------------------------
Dat4 <- data.frame(theta1, theta4)
p10 <- ggplot(data=Dat3, aes(theta1)) +
  geom_histogram(fill="grey", col="black", bins=15) +
  labs(x="Theta", y=NULL, title="Posterior of Theta 1") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
p11 <- ggplot(data=Dat3, aes(theta4)) +
  geom_histogram(fill="grey", col="black", bins=15) +
  labs(x="Theta", y=NULL, title="Posterior of Theta 4") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
multiplot(p10, p11, cols = 2)
#end seventh plot--------------------------------------

#P(x|x)
xnew <- rep(NA, times = length(theta4))
for (i in 1:length(theta4)) {
  xnew[i] <- rbinom(n=1, size = 11, theta4[i])
}

#Eighth plot------------------------------------------
Dat4 <- data.frame(xnew)
p12 <- ggplot(data=Dat4, aes(xnew)) +
  geom_histogram(fill="grey", col="black", bins=11) +
  labs(x="X=Binom(11, Theta4)", y=NULL, title="Posterior Predictive Distribution of X for Experiment 4") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
p12
#end eighth plot--------------------------------------

#P(theta5 | x)
theta5 <- rep(NA, times = length(ranAlpha))
for (i in 1:length(ranAlpha)) {
  theta5[i] <- rbeta(n=1, ranAlpha[i], ranBeta[i])
}

#P(xnew | x)
xnew <- rep(NA, times = length(theta5))
for (i in 1:length(theta5)) {
  xnew[i] <- rbinom(n=1, size = 10, theta5[i])
}

#Ninth plot----------------------------------------------
Dat4 <- data.frame(theta5, xnew)
p13 <- ggplot(data=Dat4, aes(x=theta5)) +
  geom_histogram(fill="grey", col="black", bins=15) +
  labs(x="Theta", y=NULL, title="Posterior of Theta 5") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())

p14 <- ggplot(data=Dat4, aes(x=xnew)) +
  geom_histogram(fill="grey", col="black", bins=10) +
  labs(x="X=Binom(10, Theta5)", y=NULL, title="Posterior Predictive Distribution of X") +
  theme(axis.text.y=element_blank(), axis.ticks = element_blank())
multiplot(p13, p14, cols=2)
#end ninth plot----------------------------------------------
