library(ncvreg)
library(R.matlab)
library(glmnet)

setwd("Y:/MyDocuments/Xiaoqian/GMC-computation/exp/relaxed-GMC")

results <- readMat("relaxed_rho_test.mat")

nReps <- as.numeric(results$nReps)
numRho <- as.numeric(results$numRho)
p <- as.numeric(results$p)
beta <- as.vector(results$beta)

Err_rLasso <- matrix(NA, nrow = nReps, ncol = numRho) 
Err_Lasso <- matrix(NA, nrow = nReps, ncol = numRho)
Err_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
Err_MCP <- matrix(NA, nrow = nReps, ncol = numRho) 


Err0_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
Err0_MCP <- matrix(NA, nrow = nReps, ncol = numRho)

Pd_rLasso <- matrix(NA, nrow = nReps, ncol = numRho)
Pd_Lasso <- matrix(NA, nrow = nReps, ncol = numRho) 
Pd_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
Pd_MCP <- matrix(NA, nrow = nReps, ncol = numRho) 



Pd0_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
Pd0_MCP <- matrix(NA, nrow = nReps, ncol = numRho)

F1_rLasso <- matrix(NA, nrow = nReps, ncol = numRho) 
F1_Lasso <- matrix(NA, nrow = nReps, ncol = numRho) 
F1_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
F1_MCP <- matrix(NA, nrow = nReps, ncol = numRho)

TP_rLasso <- matrix(NA, nrow = nReps, ncol = numRho)
TP_Lasso <- matrix(NA, nrow = nReps, ncol = numRho) 
TP_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
TP_MCP <- matrix(NA, nrow = nReps, ncol = numRho)

FP_rLasso <- matrix(NA, nrow = nReps, ncol = numRho)
FP_Lasso <- matrix(NA, nrow = nReps, ncol = numRho) 
FP_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
FP_MCP <- matrix(NA, nrow = nReps, ncol = numRho)

FN_rLasso <- matrix(NA, nrow = nReps, ncol = numRho) 
FN_Lasso <- matrix(NA, nrow = nReps, ncol = numRho) 
FN_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
FN_MCP <- matrix(NA, nrow = nReps, ncol = numRho)

lambda_rLasso <- matrix(NA, nrow = nReps, ncol = numRho) 
lambda_Lasso <- matrix(NA, nrow = nReps, ncol = numRho) 
alpha_rLasso <- matrix(NA, nrow = nReps, ncol = numRho) 
lambda_SCAD <- matrix(NA, nrow = nReps, ncol = numRho) 
lambda_MCP <- matrix(NA, nrow = nReps, ncol = numRho)

for (j in 1:numRho) {
  
  
  #X <- results$Xlist[[j]][[1]]
  X <- results$Xlist[j, 1][[1]][[1]]
  
  for (i in 1:nReps) {
    
    
    #y <- results$ylist[[(j-1)*nReps+i]][[1]]
    y <- results$ylist[i, j][[1]][[1]]
    
    ##Lasso
    fit_Lasso <- cv.glmnet(X, y, nfolds = 5, gamma=seq(0, 1, by=0.25), relax = TRUE)
    
    # relaxed lasso
    lambda_rLasso[i, j] <- fit_Lasso$relaxed$lambda.min
    alpha <- fit_Lasso$relaxed$gamma.min
    alpha_rLasso[i, j] <- alpha
    fit <- glmnet(X, y, lambda = lambda_rLasso[i, j], intercept=FALSE)
    beta_rLasso <- fit$beta
    
    
    TP_rLasso[i, j] = length(intersect(which(beta_rLasso!=0),which(beta!=0)));
    FP_rLasso[i, j] = length(intersect(which(beta_rLasso!=0),which(beta==0)));
    FN_rLasso[i, j] = length(intersect(which(beta_rLasso==0),which(beta!=0)));
    F1_rLasso[i, j] = 2*TP_rLasso[i, j] / (2*TP_rLasso[i, j] + FP_rLasso[i,j] + FN_rLasso[i, j]);
    
    
    if(length(which(beta_rLasso!=0))==0){
      beta_hat <- rep(0, p)
      # compute relax estimate
      beta_relax <- alpha*beta_rLasso + (1-alpha)*beta_hat
      Err_rLasso[i, j] <- norm(beta-beta_relax, '2')
      Pd_rLasso[i, j] <- norm(X%*%beta-X%*%beta_relax, '2')
    }else{
      refit <- lm(y~X[, which(beta_rLasso!=0)]-1)
      beta_hat <- rep(0, p)
      beta_hat[which(beta_rLasso!=0)] = as.vector(refit$coefficients)
      # compute relax estimate
      beta_relax <- alpha*beta_rLasso + (1-alpha)*beta_hat
      Err_rLasso[i, j] <- norm(beta-beta_relax, '2')
      Pd_rLasso[i, j] <- norm(X%*%beta-X%*%beta_relax, '2')
    }
    
    # Lasso
    lambda_Lasso[i, j] <- fit_Lasso$lambda.min
    fit <- glmnet(X, y, lambda = lambda_Lasso[i, j], intercept=FALSE)
    beta_Lasso <- fit$beta
    TP_Lasso[i, j] = length(intersect(which(beta_Lasso!=0),which(beta!=0)));
    FP_Lasso[i, j] = length(intersect(which(beta_Lasso!=0),which(beta==0)));
    FN_Lasso[i, j] = length(intersect(which(beta_Lasso==0),which(beta!=0)));
    F1_Lasso[i, j] = 2*TP_Lasso[i, j] / (2*TP_Lasso[i, j] + FP_Lasso[i,j] + FN_Lasso[i, j]);
    
    Err_Lasso[i, j] <- norm(beta-beta_Lasso, '2')
    Pd_Lasso[i, j] <- norm(X%*%beta-X%*%beta_Lasso, '2')
    
    
    
    ##SCAD
    fit_SCAD <- cv.ncvreg(X, y, penalty="SCAD", nfolds = 5)
    
    
    lambda_SCAD[i, j] <- fit_SCAD$lambda.min
    fit <- ncvfit(X, y, penalty="SCAD", lambda = fit_SCAD$lambda.min)
    beta_SCAD <- fit$beta
    
    
    TP_SCAD[i, j] = length(intersect(which(beta_SCAD!=0),which(beta!=0)));
    FP_SCAD[i, j] = length(intersect(which(beta_SCAD!=0),which(beta==0)));
    FN_SCAD[i, j] = length(intersect(which(beta_SCAD==0),which(beta!=0)));
    F1_SCAD[i, j] = 2*TP_SCAD[i, j] / (2*TP_SCAD[i, j] + FP_SCAD[i,j] + FN_SCAD[i, j]);
    
    Err0_SCAD[i, j] <- norm(beta-beta_SCAD, '2')
    Pd0_SCAD[i, j] <- norm(X%*%beta-X%*%beta_SCAD, '2')
    
    if(length(which(beta_SCAD!=0))==0){
      beta_hat <- rep(0, p)
      Err_SCAD[i, j] <- norm(beta-beta_hat, '2')
      Pd_SCAD[i, j] <- norm(X%*%beta-X%*%beta_hat, '2')
    }else{
      refit <- lm(y~X[, which(beta_SCAD!=0)]-1)
      beta_hat <- rep(0, p)
      beta_hat[which(beta_SCAD!=0)] = as.vector(refit$coefficients)
      Err_SCAD[i, j] <- norm(beta-beta_hat, '2')
      Pd_SCAD[i, j] <- norm(X%*%beta-X%*%beta_hat, '2')
    }
    
    ##MCP
    fit_MCP <- cv.ncvreg(X, y, penalty="MCP", nfolds = 5)
    
    
    lambda_MCP[i, j] <- fit_MCP$lambda.min
    fit <- ncvfit(X, y, penalty="MCP", lambda = fit_MCP$lambda.min)
    beta_MCP <- fit$beta
    
    
    TP_MCP[i, j] = length(intersect(which(beta_MCP!=0),which(beta!=0)));
    FP_MCP[i, j] = length(intersect(which(beta_MCP!=0),which(beta==0)));
    FN_MCP[i, j] = length(intersect(which(beta_MCP==0),which(beta!=0)));
    F1_MCP[i, j] = 2*TP_MCP[i, j] / (2*TP_MCP[i, j] + FP_MCP[i,j] + FN_MCP[i, j]);
    
    
    Err0_MCP[i, j] <- norm(beta-beta_MCP, '2')
    Pd0_MCP[i, j] <- norm(X%*%beta-X%*%beta_MCP, '2')
    
    if(length(which(beta_MCP!=0))==0){
      beta_hat <- rep(0, p)
      Err_MCP[i, j] <- norm(beta-beta_hat, '2')
      Pd_MCP[i, j] <- norm(X%*%beta-X%*%beta_hat, '2')
    }else{
      refit <- lm(y~X[, which(beta_MCP!=0)]-1)
      beta_hat <- rep(0, p)
      beta_hat[which(beta_MCP!=0)] = as.vector(refit$coefficients)
      Err_MCP[i, j] <- norm(beta-beta_hat, '2')
      Pd_MCP[i, j] <- norm(X%*%beta-X%*%beta_hat, '2')
    }
    
    save.image("sparse_rho_ncv_np5_snr5_rep50.RData")
  }
  
}