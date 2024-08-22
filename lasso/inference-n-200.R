library("glmnet")
library("natural")

source("functions.R")
NN <- 5
args <- commandArgs(trailingOnly = TRUE)

save.directory <- "./output"

this.id <- as.integer(args[1])

n <- 200
kappa <- 0.5


p <- n * kappa
beta <- rep(0,n*kappa)
beta[1] <- 3


which.beta <- 1
sigma_sq <- 1

pen.seq <- 10^seq(1, -3, length = 100)
pen <- pen.seq[43] ## cross-validated results
level.seq <- seq(0.05,0.95,by = 0.05)
this.cover.mat <- matrix(0,nrow = 2, ncol = length(level.seq))
this.length.vec <- matrix(0,nrow = 2, ncol = length(level.seq))

set.seed(54321)
xy <- data.simu(beta, n, kappa,sqrt(sigma_sq))
X <- data.matrix(xy[,-ncol(xy)])
XX <- t(X) %*% X
set.seed(12345*this.id)

for (nn in 1:NN){
  Y.mean <- as.numeric(X %*% beta)
  xy[,ncol(xy)]<- Y.mean + sqrt(sigma_sq) * rnorm(n)
  y <- data.matrix(xy[,ncol(xy)])[,1]
  
  if (T){
      this.fit <- nlasso_cv(x = X, y = y)
      this.sigma_sq <- this.fit$sig_df**2;print(this.sigma_sq)
      res.ML <- lasso_model(xy, pen)
      res.ML$coefficients <- coef(res.ML)[,1]
      res.ML$coefficients <- res.ML$coefficients[2:(p+1)]      
      res.ML$sigmasq <- get_sigma(xy, res.ML$beta)
  
  }else{
      res.ML <- lasso_model(xy, pen)
      res.ML$coefficients <- coef(res.ML)[,1]
      res.ML$coefficients <- res.ML$coefficients[2:(p+1)]
      res.ML$sigmasq <- get_sigma(xy, res.ML$coefficients)
      this.sigma_sq <- res.ML$sigmasq*n /(n - sum(res.ML$coefficients != 0))
  
  }

  level.target <- c(0.05,0.5,0.95)

  this.N <- 100
  this.round <- 100

  out.1 <- boot.ra (xy, beta, this.sigma_sq,sigma_sq, res.ML, sir=FALSE, m=as.integer(n), M=c(this.round,this.N),
                      level.target=level.target, pen = pen, sample_pen = 2*pen, upper_scale = 5, lower_scale= 0.75, c = 2000, decay = 1)
  ## sample_pen: the penalty value used when obtaining the resampled beta (larger penalty is used here for numerical effectiveness and stability). 
  
  save.image(file = paste0(save.directory,"/result-",this.id,"-",nn,"-n_200-lasso-dr-out.RDS"))

  F.vec <- round(c(out.1$F_mat),3)
  
  sample.N <- 1000
  
  ## In this RA setting, considering F function take values in a discrete set 0, 0.01, ..., 1.00, a slightly different resampling approach compared to the main manuscript is used here. 
  
  this.weight <- c(rep(sample.N/this.N, this.N - 1),sample.N/(2*this.N), sample.N/(2*this.N))  ## Make the number of samples having F function value 0.00, 0.01, 0.02, ... 1.00 Uniformly distributed
  
  beta.mat <- matrix(NA, nrow = sample.N, ncol = p)
  this.count <- 1
  this.seq <- round(seq(0,1,by = 1/this.N),3)
  
  for (i in 1:length(this.seq)){
    this.target <- which(F.vec == this.seq[i])
    j.seq <- generate_indices(i,length(this.seq))
    idx <- 1
    while(length(this.target) == 0){
      this.target <- which(F.vec == this.seq[j.seq[idx]])
      idx <- idx + 1
    }
    
    this.sample <- sample.effective(this.target, this.weight[i]) - 1
    
    this.col <- floor(this.sample / this.round) + 1
    this.row <- this.sample %% this.round + 1
    
    for (j in 1:length(this.row)){
      beta.mat[this.count,] <- out.1$beta_list[[this.col[j]]][this.row[j],]
      this.count <- this.count + 1
    }
    
  }
  
  cal.vec <- rep(NA, nrow(beta.mat))
  cal.vec.marginal <- rep(NA, nrow(beta.mat))
  
  
  for (i in 1:nrow(beta.mat)){
    beta_diff <- beta.mat[i,] - res.ML$coefficients
    cal.vec[i] <- sum(beta_diff*as.numeric(XX %*% beta_diff))/(p)
    cal.vec.marginal[i] <- abs(beta.mat[i,which.beta] - res.ML$coefficients[which.beta])
    
  }
  
  this.beta.result <- sum((beta - res.ML$coefficients) * as.numeric(XX %*% (beta - res.ML$coefficients)))/(p)
  
  this.beta.result.marginal <- abs(beta[which.beta] - res.ML$coefficients[which.beta])
  
  
  
  for(i in 1:length(level.seq)){
    this.length.vec[1,i] <- this.length.vec[1,i] + 2*quantile(cal.vec, level.seq[i])
    this.length.vec[2,i] <- this.length.vec[2,i] + 2*quantile(cal.vec.marginal, level.seq[i])
    
    if (this.beta.result <= quantile(cal.vec, level.seq[i])){
      this.cover.mat[1,i] <- this.cover.mat[1,i] + 1
    }
    
    if (this.beta.result.marginal <= quantile(cal.vec.marginal, level.seq[i])){
      this.cover.mat[2,i] <- this.cover.mat[2,i] + 1
      
    }
    
  }
  
  
   print(this.cover.mat/nn)

}

saveRDS(this.cover.mat, paste0(save.directory,"/result-",this.id,"-n_200-lasso.RDS"))
saveRDS(this.length.vec, paste0(save.directory,"/length-",this.id,"-n_200-lasso.RDS"))