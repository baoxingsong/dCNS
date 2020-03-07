# file score could be produced by function implemented in permutationLsqLambdaK
# this script is used to estimate the K and lambda value to calculate E-value and P-value

data1 = read.table("randomScores")
u = unique(data1$V1)
s = u[order(u)]

data = data.frame(x=s, y=s)
for ( i in 1:length(s) ){
	data[i, 2] =  length(which(data1$V1 >= data$x[i])) / nrow(data1)
}
x=data$x
y=data$y
library(minpack.lm)
nonlin_mod=nlsLM(y~(1-exp(-1*k*1000000*exp(-1*l*x))), control=nls.lm.control(maxiter=550), start=list(k=0.3, l=0.2))
plot(x, y, xlab="s", ylab="p-value")
lines(x,predict(nonlin_mod),col="red")
nonlin_mod




4000000
nonlin_mod
Nonlinear regression model
  model: y ~ (1 - exp(-1 * k * 4e+06 * exp(-1 * l * x)))
   data: parent.frame()
       k        l
0.001959 0.344325
 residual sum-of-squares: 0.003286

Number of iterations to convergence: 25
Achieved convergence tolerance: 1.49e-08



1000000
nonlin_mod
Nonlinear regression model
  model: y ~ (1 - exp(-1 * k * 1e+06 * exp(-1 * l * x)))
   data: parent.frame()
       k        l
0.006662 0.382291
 residual sum-of-squares: 0.001207

Number of iterations to convergence: 29
Achieved convergence tolerance: 1.49e-08






## for masking
1000000
nonlin_mod
Nonlinear regression model
  model: y ~ (1 - exp(-1 * k * 1e+06 * exp(-1 * l * x)))
   data: parent.frame()
      k       l
0.01011 0.38877
 residual sum-of-squares: 0.0002461

Number of iterations to convergence: 25
Achieved convergence tolerance: 1.49e-08










int matchingScore = 2;
int mismatchingPenalty = -3;
int openGapPenalty = -3;
int extendGapPenalty = -1;


data1 = read.table("score")
u = unique(data1$V1)
s = u[order(u)]

data = data.frame(x=s, y=s)
for ( i in 1:length(s) ){
	data[i, 2] =  length(which(data1$V1 >= data$x[i])) / nrow(data1)
}
x=data$x
y=data$y
library(minpack.lm)
nonlin_mod=nlsLM(y~(1-exp(-1*k*1000000*exp(-1*l*x))), control=nls.lm.control(maxiter=550), start=list(k=0.3, l=0.2))
plot(x, y, xlab="s", ylab="p-value")
lines(x,predict(nonlin_mod),col="red")

nonlin_mod

Nonlinear regression model
  model: y ~ (1 - exp(-1 * k * 1e+06 * exp(-1 * l * x)))
   data: parent.frame()
        k         l
1.807e-05 5.740e-02
 residual sum-of-squares: 0.001187

Number of iterations to convergence: 32
Achieved convergence tolerance: 1.49e-08



