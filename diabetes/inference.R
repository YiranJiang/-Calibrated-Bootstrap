install_and_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      library(package, character.only = TRUE)
    }
  }
}

packages <- c("lars","glmnet")
install_and_load(packages)

source("functions.R")

data(diabetes)
X <- diabetes$x
X <- scale(X,T,T)
y <- diabetes$y
y <- y - mean(y)

  
n <- nrow(X)
p <- ncol(X)

xy <-  data.frame(cbind(X, y))
names(xy) <- c(paste("X", 1:p, sep="."), "Y")
xnam <- paste0("X.", 1:p)
(fmla <- as.formula(paste("Y ~ 0 + ", paste(xnam, collapse= "+"))))
attr(xy, "model") <- fmla


set.seed(12345)
pen <- cv_lambda(xy)
beta.matrix <- matrix(NA,nrow = p, ncol = 1000)
F.matrix <- matrix(NA,nrow = p, ncol = 1000)

if(T){ 
  max.x <- sqrt(p/n)*5.0


  this.lm <- lm(attr(xy, "model"),data = xy)
  this.lm <- summary(this.lm)
  this.sigma_sq <- this.lm$sigma**2
  sigma_sq <- this.sigma_sq
  
  res.ML <- lasso_model(xy, pen)
  res.ML$coefficients <- coef(res.ML)[,1]
  res.ML$coefficients <- res.ML$coefficients[2:(p+1)]
  
  res.ML$coefficients

  level.target <- c(0.05,0.5,0.95)


  this.N <- 100
  this.round <- 100
  
  
  out.1 <- boot.ra (xy, beta, this.sigma_sq,sigma_sq, res.ML, sir=FALSE, m=as.integer(n), M=c(this.round,this.N),
                      level.target=level.target, pen = pen, sample_pen = 2*pen, upper_scale = 2, lower_scale= 10/n, c = 2000, decay = 1)
  ## sample_pen: the penalty value used when obtaining the resampled beta (larger penalty is used here for numerical effectiveness and stability). 
  
  
  set.seed(123)
  
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
  

  cal.mat.marginal <- matrix(NA, nrow = nrow(beta.mat), ncol = p)
  
  
  for (i in 1:nrow(beta.mat)){
    beta_diff <- abs(beta.mat[i,] - res.ML$coefficients)
    cal.mat.marginal[i,] <- beta_diff
    
  }

  alpha <- 0.05
  cb.full.ci.matrix <- matrix(NA, nrow = p, ncol = 2)
  
  for (j in 1:p){
    cb.full.ci.matrix[j,] <- c(res.ML$coefficients[j] - quantile(cal.mat.marginal[,j], 1-alpha),
                               res.ML$coefficients[j] + quantile(cal.mat.marginal[,j], 1-alpha))
  }
  

    

}

saveRDS(cb.full.ci.matrix, "cb-diabetes-result.RDS")





