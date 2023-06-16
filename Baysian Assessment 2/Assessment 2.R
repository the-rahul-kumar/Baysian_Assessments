setwd(getwd())

# Set the Seed to 1234
set.seed(1234)

# Plot the mixture distribution
x <- seq(0, 1,length=1000)
fx <- 1/3*dbeta(x,1,5) + 1/3*dbeta(x,3,5) + 1/3*dbeta(x,10,5)

f_x <- function(x) {
  new_x = 1/3*dbeta(x,1,5) + 1/3*dbeta(x,3,5) + 1/3*dbeta(x,10,5)
}

par(mfrow=c(1,1))
plot(x,fx,type="l",col="blue",lwd = 3,main = "Mixture Distribution Plot",
     xlab = "X Value Sequence", ylab = "Beta Distribution")

par(mfrow=c(1,1))
plot(x,fx,type="l",col="blue",lwd = 3, main = "Accept/Reject Algorithm on Unif(0,1)",
     xlab = "X Value Sequence", ylab = "Beta Distribution")

# K occurs at where the uniform distribution is maximum, this is to ensure
# that we have the tightest distribution and not reject excessive values.
# It can be seen that the highest point in the distribution occurs at x = 0
# therefore K is as follows

K <- 1/3*dbeta(0,1,5) + 1/3*dbeta(0,3,5) + 1/3*dbeta(0,10,5)
K
# Plot K as a straight line
abline(a = K, b = 0,col="darkorange1",lty=2,lwd = 2)

# Simulate 10,000 values in a Uniform Distribution
N <- 10000
X <- runif(N)
U <- runif(N, 0, f_x(0))

points(X,U, pch=16, cex=0.5)

accepted_x <- c()
for(i in 1:N){
  x_test <- runif(1,0,1)
  acceptance <- f_x(x_test) / (K)
  u <- runif(1,0,1)
  if(u < acceptance){
    accepted_x <- c(accepted_x,x_test)
  }
}

ind=(U < f_x(X))
points(X[ind],U[ind], cex=0.5, col="red")
par(mfrow=c(1,1))
hist(accepted_x,probability=T)
lines(x,fx,type="l",col="blue",lwd = 3)
# points(X[ind],U[ind], cex=0.5, col="red")

# The Observed acceptance rate is 
sum(ind==T)/N

# Theoretical Acceptance Rate is 1/K
1/K


# Importance Sampling Algorithm

ind=(U < f_x(X))
mean(X[ind])

w=f_x(X)/dunif(X)
W=w/sum(w) 

par(mfrow=c(1,1))
hist(X[ind],probability=T)
lines(x,fx,type="l",col="blue",lwd = 3)
d=density(X,weights=W,from=0,to=1)
lines(d,col="red",lwd=3,lty=2)

X_IS <- sample(X,size=N,prob=W,replace=T)
lines(density(X_IS),col="orange",lwd=3)
quantile(X_IS,probs=c(0.25,0.50,0.75))

