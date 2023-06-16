setwd(getwd())

library(mvtnorm)

# Read CSV
raw_data <- read.csv(file = 'winequality-red.csv')

# Find NA in data, its clean, however just incase removing them too
sum(is.na(raw_data))
data <- na.omit(raw_data) 

# Add column called good_Wine where if quality is >= 6.5 its 1, otherwise its 0
df <- transform(data,good_wine=ifelse(quality>=6.5,1,0))
# Remove the Quality Column
df <- subset(df, select=-c(quality))

model <- glm(good_wine ~ fixed.acidity+volatile.acidity+citric.acid+
              residual.sugar+chlorides+free.sulfur.dioxide+total.sulfur.dioxide+
              density+pH+sulphates+alcohol,family = binomial(link="logit") , data=df)

# Print summary to estimate useful parameters 
s <- summary(model)
s
# capture summary as a text file to input in report
#capture.output(s, file = "model_summary.txt")


# Probability of total sulfur dioxide
# Range of total.sulfur.dioxide
total.sulfur.dioxide_range <- seq(from=min(df$total.sulfur.dioxide),
                               to=max(df$total.sulfur.dioxide), by=1)

# Grab the intercepts
b0 <- model$coef[1] # Intercept
b1 <- model$coef[2] # fixed.acidity
b2 <- model$coef[3] # volatile.acidity
b3 <- model$coef[4] # citric.acid
b4 <- model$coef[5] # residual.sugar
b5 <- model$coef[6] # chlorides
b6 <- model$coef[7] # free.sulfur.dioxide
b7 <- model$coef[8] # total.sulfur.dioxide
b8 <- model$coef[9] # density
b9 <- model$coef[10] # pH
b10 <- model$coef[11] # sulphates
b11 <- model$coef[12] # alcohol

# Set the other variables to their mean
fixed.acidity_mean <- mean(df$fixed.acidity)
volatile.acidity_mean <- mean(df$volatile.acidity) 
citric.acid_mean <- mean(df$citric.acid)
residual.sugar_mean <- mean(df$residual.sugar)
chlorides_mean <- mean(df$chlorides) 
free.sulfur.dioxide_mean <- mean(df$free.sulfur.dioxide)
density_mean <- mean(df$density)
pH_mean <- mean(df$pH)
sulphates_mean <- mean(df$sulphates) 
alcohol_mean <- mean(df$alcohol)

# Write the equation
total.sulfur.dioxide_logits <- b0 + fixed.acidity_mean*b1 + 
  volatile.acidity_mean*b2 + citric.acid_mean*b3 +
  residual.sugar_mean*b4 + chlorides_mean*b5 + free.sulfur.dioxide_mean*b6 +
  total.sulfur.dioxide_range*b7 + density_mean*b8 + pH_mean*b9 + 
  sulphates_mean*b10 + alcohol_mean*b11

# Calculate the probability
total.sulfur.dioxide_probs <- exp(total.sulfur.dioxide_logits)/
  (1 + exp(total.sulfur.dioxide_logits))

plot(total.sulfur.dioxide_range, total.sulfur.dioxide_probs,
     type="l", 
     lwd=3, 
     lty=2, 
     col="blue", 
     xlab="Total Sulfur Dioxide", ylab="P(Good Wine)", main="Probability of Good Wine")


# Log posterior distribution
lpost.LR <- function(beta,x,y)
{
  eta <- as.numeric(x %*% beta)
  logp <- eta - log(1+exp(eta))
  logq <- log(1-exp(logp))
  logl <- sum(logp[y==1]) + sum(logq[y==0])
  lprior <- sum(dnorm(beta,0,10,log=T))
  return(logl + lprior)
}

# Simulations fixed at 10^4
S <- 10^4

X <- cbind(rep(1,nrow(df)), df$fixed.acidity, df$volatile.acidity, 
         df$citric.acid, df$residual.sugar, df$chlorides, df$free.sulfur.dioxide, 
         df$total.sulfur.dioxide, df$density, df$pH, df$sulphates, df$alcohol)
y <- df$good_wine

beta_mat1 <- matrix(NA,nrow=S,ncol=ncol(X))
beta_mat2 <- matrix(NA,nrow=S,ncol=ncol(X))
beta_mat3 <- matrix(NA,nrow=S,ncol=ncol(X))
beta_mat4 <- matrix(NA,nrow=S,ncol=ncol(X))

# # Initialisation 1 - MLE
beta_mat1[1,] <- as.numeric(coefficients(model))

# # Initialisation 2 - All 0
beta_mat2[1,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# Initialisation 3 - Mean
beta_mat3[1,] <- c(1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0)

# # Initialisation 4 - Unobserved Variable
beta_mat4[1,] <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)

y_new <- c(1)
x_new <- c(1,7.5,0.6, 0.0, 1.70, 0.085, 5, 45, 0.9965, 3.40, 0.63, 12)

# Initialisation Beta 1
Omega_prop <- solve(t(X) %*% X)
k <- ncol(beta_mat1)
acc <- 0
for(iter in 2:S)
{
  # 1. Propose a new set of values
  beta_star <- rmvnorm(1,beta_mat1[iter-1,],0.5*Omega_prop)
  
  # 2. Compute the posterior density on the proposed value and on the old value  
  newpost=lpost.LR(t(beta_star),X,y)
  oldpost=lpost.LR(matrix(beta_mat1[iter-1,],ncol=1),X,y)
  
  # 3. Acceptance step
  if(runif(1,0,1)>exp(newpost-oldpost)){
    beta_mat1[iter,]=beta_mat1[iter-1,]
  } else{
    beta_mat1[iter,]=beta_star
    acc=acc+1
  }
  # 4. Print the stage of the chain
  if(iter%%1000==0){print(c(iter,acc/iter))}
  
  # 5. Prediction 
  p_new <- exp(sum(beta_mat1[iter,] * x_new) ) / (1 + exp(sum(beta_mat1[iter,] * x_new) ))
  y_new[iter] <- rbinom(1,1,prob=p_new)
}

# Initialisation Beta 2
k <- ncol(beta_mat2)
acc <- 0
for(iter in 2:S)
{
  # 1. Propose a new set of values
  beta_star <- rmvnorm(1,beta_mat2[iter-1,],0.5*Omega_prop)
  
  # 2. Compute the posterior density on the proposed value and on the old value  
  newpost=lpost.LR(t(beta_star),X,y)
  oldpost=lpost.LR(matrix(beta_mat2[iter-1,],ncol=1),X,y)
  
  # 3. Acceptance step
  if(runif(1,0,1)>exp(newpost-oldpost)){
    beta_mat2[iter,]=beta_mat2[iter-1,]
  } else{
    beta_mat2[iter,]=beta_star
    acc=acc+1
  }
  # 4. Print the stage of the chain
  if(iter%%1000==0){print(c(iter,acc/iter))}
  
  # 5. Prediction 
  p_new <- exp(sum(beta_mat2[iter,] * x_new) ) / (1 + exp(sum(beta_mat2[iter,] * x_new) ))
  y_new[iter] <- rbinom(1,1,prob=p_new)
}


# Initialisation Beta 3
k <- ncol(beta_mat3)
acc <- 0
for(iter in 2:S)
{
  # 1. Propose a new set of values
  beta_star <- rmvnorm(1,beta_mat3[iter-1,],0.5*Omega_prop)
  
  # 2. Compute the posterior density on the proposed value and on the old value  
  newpost=lpost.LR(t(beta_star),X,y)
  oldpost=lpost.LR(matrix(beta_mat3[iter-1,],ncol=1),X,y)
  
  # 3. Acceptance step
  if(runif(1,0,1)>exp(newpost-oldpost)){
    beta_mat3[iter,]=beta_mat3[iter-1,]
  } else{
    beta_mat3[iter,]=beta_star
    acc=acc+1
  }
  # 4. Print the stage of the chain
  if(iter%%1000==0){print(c(iter,acc/iter))}
  
  # 5. Prediction 
  p_new <- exp(sum(beta_mat3[iter,] * x_new) ) / (1 + exp(sum(beta_mat3[iter,] * x_new) ))
  y_new[iter] <- rbinom(1,1,prob=p_new)
}

# Initialisation Beta 4
k <- ncol(beta_mat4)
acc <- 0
for(iter in 2:S)
{
  # 1. Propose a new set of values
  beta_star <- rmvnorm(1,beta_mat4[iter-1,],0.5*Omega_prop)
  
  # 2. Compute the posterior density on the proposed value and on the old value  
  newpost=lpost.LR(t(beta_star),X,y)
  oldpost=lpost.LR(matrix(beta_mat4[iter-1,],ncol=1),X,y)
  
  # 3. Acceptance step
  if(runif(1,0,1)>exp(newpost-oldpost)){
    beta_mat4[iter,]=beta_mat4[iter-1,]
  } else{
    beta_mat4[iter,]=beta_star
    acc=acc+1
  }
  # 4. Print the stage of the chain
  if(iter%%1000==0){print(c(iter,acc/iter))}
  
  # 5. Prediction 
  p_new <- exp(sum(beta_mat4[iter,] * x_new) ) / (1 + exp(sum(beta_mat4[iter,] * x_new) ))
  y_new[iter] <- rbinom(1,1,prob=p_new)
}
# Plot the Graphs
plot(beta_mat1[,1],type="l", ylab=expression(beta[0]))
lines(beta_mat2[,1],col="red")
lines(beta_mat3[,1],col="blue")
lines(beta_mat4[,1],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[1],col="red",lty=2)

plot(beta_mat1[,2],type="l", ylab=expression(beta[1]))
lines(beta_mat2[,2],col="red")
lines(beta_mat3[,2],col="blue")
lines(beta_mat4[,2],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[2],col="red",lty=2)

plot(beta_mat1[,3],type="l", ylab=expression(beta[2]))
lines(beta_mat2[,3],col="red")
lines(beta_mat3[,3],col="blue")
lines(beta_mat4[,3],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[3],col="red",lty=2)

plot(beta_mat1[,4],type="l", ylab=expression(beta[3]))
lines(beta_mat2[,4],col="red")
lines(beta_mat3[,4],col="blue")
lines(beta_mat4[,4],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[1],col="red",lty=2)
abline(h=model$coefficients[4],col="red",lty=2)

plot(beta_mat1[,5],type="l", ylab=expression(beta[4]))
lines(beta_mat2[,5],col="red")
lines(beta_mat3[,5],col="blue")
lines(beta_mat4[,5],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[1],col="red",lty=2)
abline(h=model$coefficients[5],col="red",lty=2)

plot(beta_mat1[,6],type="l", ylab=expression(beta[5]))
lines(beta_mat2[,6],col="red")
lines(beta_mat3[,6],col="blue")
lines(beta_mat4[,6],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[5],col="red",lty=2)

plot(beta_mat1[,7],type="l", ylab=expression(beta[6]))
lines(beta_mat2[,7],col="red")
lines(beta_mat3[,7],col="blue")
lines(beta_mat4[,7],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[7],col="red",lty=2)

plot(beta_mat1[,8],type="l", ylab=expression(beta[7]))
lines(beta_mat2[,8],col="red")
lines(beta_mat3[,8],col="blue")
lines(beta_mat4[,8],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[8],col="red",lty=2)

plot(beta_mat1[,9],type="l", ylab=expression(beta[8]))
lines(beta_mat2[,9],col="red")
lines(beta_mat3[,9],col="blue")
lines(beta_mat4[,9],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[9],col="red",lty=2)

plot(beta_mat1[,10],type="l", ylab=expression(beta[9]))
lines(beta_mat2[,10],col="red")
lines(beta_mat3[,10],col="blue")
lines(beta_mat4[,10],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[10],col="red",lty=2)

plot(beta_mat1[,11],type="l", ylab=expression(beta[10]))
lines(beta_mat2[,11],col="red")
lines(beta_mat3[,11],col="blue")
lines(beta_mat4[,11],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[11],col="red",lty=2)

plot(beta_mat1[,12],type="l", ylab=expression(beta[11]))
lines(beta_mat2[,12],col="red")
lines(beta_mat3[,12],col="blue")
lines(beta_mat4[,12],col="green")
legend("topright",inset = 0.02, legend=c("Beta_Mat1", "Beta_Mat2", "Beta_Mat3", "Beta_Mat4"),
       col=c("black", "red","blue","green"), lty=1, cex=0.8)
abline(h=model$coefficients[12],col="red",lty=2)