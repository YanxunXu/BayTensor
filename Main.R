library(MultiwayRegression)
library(MASS)
library(mvnfast)
source("Functions.R")
source("BayTensor_MCMC.R")
source('BayTensor_Fast.R')
source('Generate_Data.R')

#### Simulation Setups ##########
start_time<-proc.time()
t_true<- 3; s_true<- 3 # core tensor truth (t_true,t_true,s_true,s_true)
t_and_s <- c(t_true, s_true)

N<-100 # Number of observations
SNR <- 5
dimX<-c(16,12);dimY<-c(10,8)

seed<- 12345
data<-Generate_data(N, dimX, dimY, seed, 0, t_and_s, SNR)
## Pick core_dim ##
core_dim<-c(3,3,2,2) ## Initial guess
K <- c(t_true, t_true, s_true, s_true) + 2

##############################################
# Run BayTensor_Fast Algorithm
##############################################
sim <- BayTensor_Fast(data$Y, data$X, core_dim, 90, K)
# Prediction on test data  
rpe = NULL
Nnew = 1000
for(test in c(1:5)){
  seed<-floor(test*pi*100)
  test_data<-Generate_data(Nnew, dimX, dimY,seed, data$B, t_and_s, SNR)
  
  Y_est <- tensor_tensor(test_data$X, sim$B_est ,c(2:(length(dimX)+1)),c(1:length(dimX)))
  Y_est_vec <- as.vector(Y_est)
  errY <- Y_est_vec - as.vector(test_data$Y)
  
  bbbb<-sum(errY^2)/sum(as.vector(test_data$Y)^2)
  
  rpe = c(rpe, bbbb)
}
cat('Mean RPE: ', mean(rpe), ', SD: ', sd(rpe))




##############################################
# Run BayTensor_MCMC Algorithm
##############################################
b<- 0.15 
I<-100 # Select Core Tensor Dimension
lambda <- 1 # Possion Parameter
n_iter <- 1100
K <- c(t_true+2,t_true+2,s_true+2,s_true+2)
core_dim <- c(3,3,2,2)

sim <- BayTensor_MCMC(data$Y, data$X, core_dim, seed, b, n_iter, K, lambda, I)

rpe = NULL
Nnew = 1000
for(test in c(1:5)){
  seed<-floor(test*pi*100)
  test_data<-Generate_data(Nnew, dimX, dimY,seed, data$B, t_and_s, SNR)
  
  Y_est_vec <- rep(0, prod(dim(test_data$Y)))
  for(k in c(1:1000)){
    Y_est <- tensor_tensor(test_data$X, sim$B_estimate[[k]] ,c(2:(length(dimX)+1)),c(1:length(dimX)))
    Y_est_vec <- Y_est_vec + as.vector(Y_est) + sqrt(sim$Sigma2_estimate[[k]]) * rnorm(prod(dim(test_data$Y)))
  }
  Y_est_vec <- Y_est_vec / 1000
  errY <- Y_est_vec - as.vector(test_data$Y)
  
  bbbb<-sum(errY^2)/sum(as.vector(test_data$Y)^2)
  
  rpe = c(rpe, bbbb)
}
cat('Mean RPE: ', mean(rpe), ', SD: ', sd(rpe))
