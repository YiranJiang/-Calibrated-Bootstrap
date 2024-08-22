library("glmnet")
NN <- 1000
M <- 1000
which.beta <- 1
this.betas <- rep(NA,NN)
this.beta_hats <- matrix(nrow = NN, ncol = M)
this.beta_hats_bt_joint <- matrix(nrow = NN, ncol = M)
this.beta_hats_bt_marginal <- matrix(nrow = NN, ncol = M)
this.beta_hats_rb_joint <- matrix(nrow = NN, ncol = M)
this.beta_hats_rb_marginal <-  matrix(nrow = NN, ncol = M)

sample.residual <- function(xy,res.ML, pred.ML){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  n <- length(Y)
  p <- ncol(X)
  centered.residuals <- res.ML$residuals - mean(res.ML$residuals)
  this.residuals <- sample(centered.residuals,size = n, replace = TRUE)
  this.Y <- pred.ML + this.residuals
  xy$Y <- this.Y
  xy
}

standard_bt <- function(xy, pen, NN){
  bs_joint <- rep(NA,NN)
  bs_marginal <- rep(NA,NN)
  
  for(i in 1:NN){
    
    idx <- sample.int(n, size=n, replace=TRUE)
    
    XY <- xy[idx,]
    
    
    res.b <- lasso_model(XY, pen)
    res.b$coefficients <- coef(res.b)[,1]
    res.b$coefficients <- res.b$coefficients[2:(p+1)]
    
    bs_marginal[i] <- abs(res.b$coefficients[which.beta] - res.ML$coefficients[which.beta])
    
    beta_diff <- res.b$coefficients - res.ML$coefficients
    bs_joint[i] <- sum(beta_diff*as.numeric(XX %*% beta_diff))/(p)
    
  }
  return(list(bs_marginal = bs_marginal,bs_joint = bs_joint))
  
}

residual_bt <- function(xy, res.ML, pred.ML, pen, NN){
  bs_joint <- rep(NA,NN)
  bs_marginal <- rep(NA,NN)
  for (i in 1:NN){
    XY <- sample.residual(xy, res.ML, pred.ML)
    res.b <- lasso_model(XY,pen)
    res.b$coefficients <- coef(res.b)[,1]
    res.b$coefficients <- res.b$coefficients[2:(p+1)]
    
    bs_marginal[i] <- abs(res.b$coefficients[which.beta] - res.ML$coefficients[which.beta])
    
    beta_diff <- res.b$coefficients - res.ML$coefficients
    bs_joint[i] <- sum(beta_diff*as.numeric(XX %*% beta_diff))/(p)
  }
  
  return(list(bs_marginal = bs_marginal,bs_joint = bs_joint))
}

## Calculate Residual Bootstrap Results

this.nn <- 1
level.seq <- seq(0.05,0.95,by = 0.05)

all.cover.mat <- matrix(0,nrow = 4, ncol = length(level.seq))
all.length.vec <- matrix(0,nrow = 4, ncol = length(level.seq))

all.length.joint.bt <- matrix(0,nrow = 5, ncol = length(level.seq))
all.length.joint.rb <- matrix(0,nrow = 5, ncol = length(level.seq))

all.length.marginal.bt <- matrix(0,nrow = 5, ncol = length(level.seq))
all.length.marginal.rb <- matrix(0,nrow = 5, ncol = length(level.seq))

args <- commandArgs(trailingOnly = TRUE)

save.directory <- "./output-bootstrap"

n <- as.integer(args[1])

this.id <- as.integer(args[2])


cat("Loading ",this.id,"...\n")
for (jj in 1:5){
  set.seed(12345)
  load(paste0("output_n_",n,"/result-",this.id,"-",jj,"-n_",n,"-lasso-dr-out.RDS")) ## Load the data set from the corresponding CB experiment

  
  pred.ML <- (data.matrix(xy[,-ncol(xy)]) %*% res.ML$coefficients)[,1]
  res.ML$residuals <- xy$Y - pred.ML
  
  this.result <- residual_bt(xy, res.ML, pred.ML, pen, 1000)
  
  this.beta_hats_rb_joint[this.nn,] <- this.result$bs_joint
  this.beta_hats_rb_marginal[this.nn,] <- this.result$bs_marginal
  
  
  set.seed(12345)
  
  
  this.result <- standard_bt(xy, pen, 1000)
  
  this.beta_hats_bt_joint[this.nn,] <- this.result$bs_joint
  this.beta_hats_bt_marginal[this.nn,] <- this.result$bs_marginal
  
  
  this.beta.result <- sum((beta - res.ML$coefficients) * as.numeric(XX %*% (beta - res.ML$coefficients)))/(p)
  
  this.beta.result.marginal <- abs(beta[which.beta] - res.ML$coefficients[which.beta])
  
  for(i in 1:length(level.seq)){
    all.length.vec[1,i] <- all.length.vec[1,i] + 2*quantile(this.beta_hats_bt_joint[this.nn,], level.seq[i])
    all.length.vec[2,i] <- all.length.vec[2,i] + 2*quantile(this.beta_hats_bt_marginal[this.nn,], level.seq[i])
    all.length.vec[3,i] <- all.length.vec[3,i] + 2*quantile(this.beta_hats_rb_joint[this.nn,], level.seq[i])
    all.length.vec[4,i] <- all.length.vec[4,i] + 2*quantile(this.beta_hats_rb_marginal[this.nn,] , level.seq[i])
    
    all.length.joint.bt[jj,i] <- 2*quantile(this.beta_hats_bt_joint[this.nn,], level.seq[i])
    all.length.joint.rb[jj,i] <- 2*quantile(this.beta_hats_rb_joint[this.nn,], level.seq[i])
    all.length.marginal.bt[jj,i]  <- 2*quantile(this.beta_hats_bt_marginal[this.nn,], level.seq[i])
    all.length.marginal.rb[jj,i]  <- 2*quantile(this.beta_hats_rb_marginal[this.nn,], level.seq[i])
    
    
    if (this.beta.result <= quantile(this.beta_hats_bt_joint[this.nn,], level.seq[i])){
      all.cover.mat[1,i] <- all.cover.mat[1,i] + 1
    }
    
    if (this.beta.result.marginal <= quantile(this.beta_hats_bt_marginal[this.nn,], level.seq[i])){
      all.cover.mat[2,i] <- all.cover.mat[2,i] + 1
      
    }
    
    if (this.beta.result <= quantile(this.beta_hats_rb_joint[this.nn,], level.seq[i])){
      all.cover.mat[3,i] <- all.cover.mat[3,i] + 1
    }
    
    if (this.beta.result.marginal <= quantile(this.beta_hats_rb_marginal[this.nn,], level.seq[i])){
      all.cover.mat[4,i] <- all.cover.mat[4,i] + 1
      
    }
    
  }
  
  this.nn <- this.nn + 1
  
} 

save.directory <- "./output-bootstrap"

saveRDS(all.cover.mat, paste0(save.directory,"/result-",this.id,"-n_",n,"-lasso.RDS"))
saveRDS(list(all.length.vec = all.length.vec, all.length.joint.bt = all.length.joint.bt, all.length.joint.rb = all.length.joint.rb,
             all.length.marginal.bt = all.length.marginal.bt, all.length.marginal.rb = all.length.marginal.rb),
        paste0(save.directory,"/length-",this.id,"-n_",n,"-lasso.RDS"))


