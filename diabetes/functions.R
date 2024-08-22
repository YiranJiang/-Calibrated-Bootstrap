library(glmnet)

data.simu <- function(beta, n=5000,kappa=0.5,sigma = 3) {
  p <- as.integer(floor(kappa*n))
  X <- matrix(rnorm(n*p), ncol=p, nrow=n)
  X <- scale(X)
  col_means <- colMeans(X)
  col_sds <- apply(X, 2, sd)

  Y <- X%*% beta + rnorm(n,sd = sigma)
  data <-  data.frame(cbind(X, Y))
  names(data) <- c(paste("X", 1:p, sep="."), "Y")
  xnam <- paste0("X.", 1:p)
  (fmla <- as.formula(paste("Y ~ 0 + ", paste(xnam, collapse= "+"))))
  attr(data, "model") <- fmla
  data
}

my.lm <- function(formula, data, ...){
  res <- lm(formula, data, ...)
  res$beta <- res$coefficients
  res$sd <- sqrt(sum(res$residuals^2)/length(res$residuals))
  res
}

lasso_model <- function(xy, pen){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  lasso_model <- glmnet(X, Y, intercept = F, alpha=1, lambda = pen, standardize = F)  # alpha=1 specifies Lasso
  return(lasso_model)
}

cv_lambda <- function(xy){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  lambda_seq <- 10^seq(1, -3, length = 100)
  cv_result <- cv.glmnet(X, Y, alpha = 1, intercept = F, lambda = lambda_seq)
  optimal_lambda <- cv_result$lambda.min
  cat("optimal Lambda: ",optimal_lambda, "\n")
  return(optimal_lambda)
}


max.x <- 0.5

get_sigma <- function(xy, beta){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  return(sum((X%*% beta - Y)**2)/nrow(X))
}


loss <- function(beta, sigmasq, xy, pen=0){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  n <- length(Y) 
  z <- Y - as.numeric(X%*%beta)

  (t(z) %*% z)[1,1]/(2) + n*pen*sum(abs(beta))
}



plot.level <- function(xy, max.x, level.target, ...){
  
  #  res <- lm(attr(xy, "model"), data=xy)
  #  beta <- res$coefficients
  #  plot(c(0, max.x), c(0,1), type="n", xlab="t", ylab="Plausibility", ...)
  plot(c(0, max.x), c(0,1), type="n", xlab="t", ylab="Probability", ...)
  abline(h=level.target, col="gray")
  #  axis(4, at = level.target, labels = level.target)
}

SIR.sample <- function(xy, res.b, res.ML, size){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  n <- length(Y) #print(head(X)) print(beta)
  
  z <- Y - as.numeric(X%*%res.ML$coefficients)
  s <- sqrt(res.ML$sigmasq)
  lnDensity.ML <- dnorm(z,sd=s,log=TRUE)
  
  z <- Y - as.numeric(X%*%res.b$coefficients)
  s <- sqrt(res.b$sigmasq)
  lnDensity.b <- dnorm(z,sd=s,log=TRUE)
  
  weights <- exp(lnDensity.b - lnDensity.ML)
  weights <- weights/sum(weights)
  idx <- sample(1:n, size=size, replace=TRUE, prob = weights)
  #  cat("\tlength(unique(idx)):", length(unique(idx)),"\n")
  #  xy[idx,]
  structure(xy[idx,], weights=weights)
}

# q = 0: Standard Bootstrap
my.sample <- function(x, size, replace=TRUE, prob=NULL, q=0){
  y <- runif(size)<q
  z <- sample(x, size=size-sum(y), replace = TRUE)
  c(x[y], z)
}

sample.effective <- function(array,n){
  if (length(array) >= n){
    return(sample(array,n,replace = FALSE))
  }else{
    
    return(c(rep(array, floor(n / length(array))), sample(array, n %% length(array),replace = FALSE)))
  }
}

generate_indices <- function(i, n) {
  # Initialize the list of valid indices
  valid_indices <- c()
  
  # Maximum offset to consider; can be adjusted based on n
  max_offset <- n - 1
  
  # Generate indices dynamically
  for (offset in 1:max_offset) {
    if (i + offset <= n) valid_indices <- c(valid_indices, i + offset)
    if (i - offset >= 1) valid_indices <- c(valid_indices, i - offset)
  }
  
  # Return sorted list of unique valid indices
  return(valid_indices)
}

SA.qqplot <- function(level.m, xy, model, res.ML, XX, add=TRUE, M=1000, ...) {
  n <- nrow(xy)
  p <- ncol(xy)-1
  y <- numeric(0)
  x <- qchisq(1-level.m$pl, df=p)/p
  for(k in 1:nrow(level.m)){
    cat(k, "x:", x[k],"\n")
    log.m <- level.m$weighted[k]
    m <- exp(log.m)
    m <- as.integer(floor(m) + (runif(1) <= m - floor(m)))
    S <- numeric(0)
    for(i in 1:M){
      idx <- sample.int(n, size=m, replace=TRUE)
      XY <- xy[idx,]
      res.b <- lm(model, data=XY)
      
      b <- res.b$coefficients - res.ML$coefficients
      a <- sum(b*as.numeric(XX %*% b))/p
      S[i] <- a
    }
    y[k] <- quantile(S, probs = 1-level.m$pl[k])
    cat(k, "y:", y[k],"\n")
  }
  if(add){
    points(x, y, ...)
  } else {
    plot(x, y, xlab="Theoretical Quantile", ylab="Empirical Quantile")
    abline(a=0, b=1, col=2)
  }
}

SA.qqplot2 <- function(level.m, xy, model, res.ML, XX, Counts, beta.list, add=TRUE, M=1000, ...) {
  n <- nrow(xy)
  p <- ncol(xy)-1
  y <- numeric(0)
  
  prob <- 1/Counts
  prob <- rep(1/length(Counts), length(Counts)) #prob/sum(prob)
  
  S <- numeric(0)
  my.Counts <- Counts - Counts
  for(i in 1:M){
    idx <- sample(1:length(prob), size=1, prob=prob)
    my.Counts[idx] <- my.Counts[idx] + 1
    beta.matrix <- beta.list[[idx]]
    b <- beta.matrix[sample.int(nrow(beta.matrix), size=1),]
    b <- b - res.ML$coefficients
    a <- sum(b*as.numeric(XX %*% b))/p
    S[i] <- a
  }
  
  x <- qchisq(((1:M)-0.5)/M, df=p)/p
  y <- sort(S)
  if(add){
    points(x, y, ...)
  } else {
    plot(x, y, xlab="Theoretical Quantile", ylab="Empirical Quantile")
    abline(a=0, b=1, col=2)
  }
  print(rbind(Counts, my.Counts))
  
  if(FALSE){
    prob[1:2] <- prob[length(prob)-c(0,1)] <- 0
    prob <- prob/sum(prob)
    
    S <- numeric(0)
    my.Counts <- Counts - Counts
    for(i in 1:M){
      idx <- sample(1:length(prob), size=1, prob=prob)
      my.Counts[idx] <- my.Counts[idx] + 1
      beta.matrix <- beta.list[[idx]]
      b <- beta.matrix[sample.int(nrow(beta.matrix), size=1),]
      b <- b - res.ML$coefficients
      a <- sum(b*as.numeric(XX %*% b))/p
      S[i] <- a
    }
    
    x <- qchisq(((1:M)-0.5)/M, df=p)/p
    y <- sort(S)
    plot(x, y, ..., xlab="Theoretical Quantile", ylab="Empirical Quantile")
    abline(a=0, b=1, col=2)
    print(rbind(Counts, my.Counts))
  }
  
}

distributional.fixed <- function(m, xy, beta, sigma_sq,true_sigma_sq, model, res.ML, loss.ML, XX, sir, x.matrix,pen=0, add=TRUE, M=c(1000,100), ...) {
  x <- xy[,-ncol(xy), drop=FALSE]
  XX <- t(X)%*%X
  y <- xy[,ncol(xy)]
  n <- nrow(xy)
  p <- ncol(xy)-1
  F.list <- list()

  losses <- c()
  results <- c()
  for(b in 1:M[1]){
    print(b)
    
    if (length(m) == 1){
      idx <- sample.int(n, size=m, replace=TRUE)
    }else{
      idx <- sample.int(n, size=sample(as.numeric(m), size=1), replace=TRUE)
    }
    print(length(unique(idx)))
    XY <- xy[idx,]

    res.b <- lasso_model(XY, pen)
    res.b$coefficients <- coef(res.b)[,1]
    res.b$coefficients <- res.b$coefficients[2:(p+1)]

    loss.B <- loss(res.b$coefficients, sigma_sq, xy, pen=pen)

    losses <- c(losses, loss.B)
    
    S.y.theta.b <- -(loss.B-loss.ML)
    
    results <- c(results,S.y.theta.b/p)
    
    S <- numeric(0)
    for(db in 1:M[2]){
      while(TRUE){

        xy.simu <- if(sir) {
          SIR.sample(xy, res.b, res.ML, size=n)
        } else  {
          Y.mean <- as.numeric(x.matrix %*% res.b$coefficients)
          cbind(x, Y=Y.mean + sqrt(sigma_sq)*rnorm(n))
        }
        
        res.db <- lasso_model(xy.simu, pen)
        res.db$coefficients <- coef(res.db)[,1]
        res.db$coefficients <- res.db$coefficients[2:(p+1)]
        
        res.db$sigmasq <- get_sigma(xy.simu, res.db$coefficients)

        loss.DB.ML <- loss(res.db$coefficients, sigma_sq, xy.simu, pen=pen)
        loss.B <- loss(res.b$coefficients, sigma_sq, xy.simu, pen=pen)
        
        S[db] <- -(loss.B-loss.DB.ML)
        
        if(!is.na(S[db])) break
        cat("\n\033[0;31mloss.DB.ML =", loss.DB.ML, "\n")
        cat("\033[0;31mloss.DB    =", loss.DB.ML, "\n")
        cat("\033[0;31mS[",db, "] = NA, q=", q, "\nres.db$coefficients:\n")
        print(res.db$coefficients[is.na(res.db$coefficients)])
        cat("\033[0m")
        
        if(any(is.na(w)))
          cat("\033[0;31mw[",db, "] = ", w[is.na(w)], "\033[0m\n")
      }
    }
    
    F.S.y.theta.b <- sum(S <= S.y.theta.b)/length(S)
    
    F.list <- append(F.list, list(list(value=F.S.y.theta.b, beta=res.b$coefficients,
                                       S = S.y.theta.b, check2 = S.y.theta.b/p,
                                       sigma = sqrt(sigma_sq)
    )))
  }
  
  
  NN <- M[1]
  true.SS <- rep(NA, NN)
  true.SS2 <- rep(NA, NN)
  true.beta <- matrix(NA,nrow = NN, ncol = p)
  

  true.SS <- rep(NA, NN)
  for (i in 1:NN){
    new_xy <- xy
    new_xy[,ncol(new_xy)] <- (x.matrix %*% beta)[,1] + sqrt(true_sigma_sq) * rnorm(n)
    this.b <- lasso_model(new_xy, pen)
    this.b$coefficients <- coef(this.b)[,1]
    this.b <- this.b$coefficients[2:(p+1)]
    loss.B <- loss(beta, sigma_sq, xy, pen=pen)
    loss.ML <- loss(this.b, sigma_sq, xy, pen=pen)
    
    true.SS[i] <- (loss.B-loss.ML)/p
    
    true.SS2[i] <- true.SS[i]/get_sigma(new_xy, this.b)
    true.beta[i,] <- this.b
  }
    
  true.SS <- sort(true.SS)
  true.SS2 <- sort(true.SS2)

  x <- unlist(lapply(F.list, function(x) x$S))
  y <- unlist(lapply(F.list, function(x) x$value))
  plot(x,y, ylim=c(0,1), main="(a)", xlab=expression(T[y, theta]), ylab=expression(F[theta](T[y,theta])))
  
  z <- lapply(F.list, function(x) c(x$value, x$check2))

  z <- unlist(z)
  dim(z) <- c(2, length(z)/2)
  y <- numeric(0)
  K <- numeric(0)
  param.mat <- matrix(NA,nrow = M[1], ncol = p+1)
  zz <- numeric(0)
  for(i in 1:M[1]){
    U <- runif(1)
    k <- which.min(abs(U-z[1,]))
    y[i] <- z[2,k]
    zz[i] <- z[1,k]
    K[i] <- k
    param.mat[i,1:p] <- F.list[[k]]$beta
    param.mat[i,p+1] <- F.list[[k]]$sigma
  }
  x <- qunif(((1:ncol(z))-0.5)/ncol(z))
  plot(true.SS, sort(y), main="(b)",
       xlab="Theoretical Quantile", ylab="Empirical Quantile")
  abline(a=0, b=1, col="gray")
  
  list(F.list=F.list, K=K, param.mat = param.mat, true.beta = true.beta)
}





boot.ra <- function(xy, beta, sigma_sq, true_sigma_sq, res.ML, M =c(1000L, 200L), sir = TRUE, m, level.target=0.05,
                      n = nrow(xy),
                      p = ncol(xy)-1,
                      pen= 0, #p,
                      sample_pen = 0,
                      level.fun = NULL, # function() runif(1)
                      distributional.plot=FALSE,
                      upper_scale = 2,
                      lower_scale = 0.5,
                      c = 2000,
                      decay = 1
){
  m_0 <- m
  log.m <- log(m)
  model <- attr(xy, "model")
  if(missing(m)) m <- n 
  
  x <- xy[,-ncol(xy)]
  y <- xy[, ncol(xy)]
  x.matrix <- data.matrix(x)
  XX <- t(x.matrix)%*%x.matrix
  
  max.x <- M[1]


  do.SA <- T
  if(do.SA ){
    SA.i <-  0
    log.m <- log(m)
  }
  
  

  weight <- 1.0/min(50,M[1])
  level.m <- data.frame(pl=level.target, m=rep(m, length(level.target)),
                     log.m = rep(log.m, length(level.target)),
                     weighted = rep(log.m, length(level.target)),
                     m0 = rep(m, length(level.target)))
  sa.a  <- lapply(as.list(1:length(level.target)), function(x) numeric(0))
  sa.b  <- lapply(as.list(1:length(level.target)), function(x) numeric(0))
  sa.pl <- lapply(as.list(1:length(level.target)), function(x) numeric(0))
  

  loss.ML <- loss(beta = res.ML$coefficients, sigma_sq, xy, pen=pen)

  if(distributional.plot){
    out <- distributional.fixed(m, xy, beta, sigma_sq, true_sigma_sq, model, res.ML, loss.ML,
                          XX=t(x.matrix)%*%x.matrix,
                          sir, x.matrix, pen = pen, add=TRUE, M=M, col=4)
    return(out)
  }
  pl <- numeric(0)
  
  Counts <- rep(0, length(level.target))
  names(Counts) <- level.target
  N0 <- as.integer(M[1]/2)
  F.levels <- c((level.target[-1] + level.target[-length(level.target)])/2, 1.0)
  beta.list <- lapply(as.list(1:length(level.target)), function(x)
    matrix(nrow=0, ncol=p))
  
  loss_mat <- matrix(NA, nrow = M[1], ncol = nrow(level.m))
  result_mat <- matrix(NA, nrow = M[1], ncol = nrow(level.m))
  F_mat <- matrix(NA, nrow = M[1], ncol = nrow(level.m))
  beta_list <- list()
  
  for (k in 1:nrow(level.m)){
    beta_list[[k]] <- matrix(NA, nrow = M[1], ncol = p)
  }
  
  for(b in 1:M[1]){
    
    SA.i <-  SA.i + 1
    
    for(k in 1:nrow(level.m)){
      
      {
        log.m <- level.m$log.m[k]
        m0 <- level.m$m0[k]
        m <- exp(log.m)
        m <- as.integer(floor(m)+(runif(1)<= m-floor(m)))
        
        level.target <- level.m$pl[k]
        SA.s <- sqrt(M[2]*level.target*(1-level.target))
      }
      
      if(FALSE && !is.null(level.fun)) {
        level.target <- level.fun()
        if(b > M[1]/4) {
          if(nrow(level.m)==1){
            level.m <- data.frame(pl=0.5, m=m, log.m = log.m)
          }
          log.m <- level.m$log.m[k <- which.min(abs(level.m$pl - level.target))]
          m <- exp(log.m)
          m <- as.integer(floor(m)+(runif(1)<= m-floor(m)))
        }
      }
      
      while(TRUE){
        idx <- sample.int(n, size=m, replace=TRUE)

        XY <- xy[idx,]
        
        res.b <- lasso_model(XY, sample_pen)

        
        res.b$coefficients <- coef(res.b)[,1]
        res.b$coefficients <- res.b$coefficients[2:(p+1)]
        beta_list[[k]][b,] <- res.b$coefficients
        loss.B <- loss(res.b$coefficients, sigma_sq, xy, pen=pen)

        
        S.y.theta.b <- -(loss.B-loss.ML)

        
        if(!is.na(S.y.theta.b)) break;
        stop("Whoops")
      }
      
      if(is.na(S.y.theta.b)) stop("S.y.theta.b is NA")
      
      
      S <- numeric(0)
      for(db in 1:M[2]){
        while(TRUE){
          xy.simu <- if(sir) {
            SIR.sample(xy, res.b, res.ML, size=n)
          } else  {
            Y.mean <- as.numeric(x.matrix %*% res.b$coefficients)
            cbind(x, Y=Y.mean + sqrt(sigma_sq)*rnorm(n))
          }
          
          res.db <- lasso_model(xy.simu, pen)
          res.db$coefficients <- coef(res.db)[,1]
          res.db$coefficients <- res.db$coefficients[2:(p+1)]
          loss.DB.ML <- loss(res.db$coefficients, sigma_sq, xy.simu, pen=pen)
          loss.B <- loss(res.b$coefficients, sigma_sq, xy.simu, pen=pen)
          
          S[db] <- -(loss.B-loss.DB.ML)
          if(!is.na(S[db])) break
          cat("\n\033[0;31mloss.DB.ML =", loss.DB.ML, "\n")
          cat("\033[0;31mloss.DB    =", loss.DB.ML, "\n")
          cat("\033[0;31mS[",db, "] = NA, q=", q, "\nres.db$coefficients:\n")
          print(res.db$coefficients[is.na(res.db$coefficients)])
          cat("\033[0m")
          
          if(any(is.na(w)))
            cat("\033[0;31mw[",db, "] = ", w[is.na(w)], "\033[0m\n")
        }
      }
      
      sa.pl[[k]][b] <- F.S.y.theta.b <- sum(S <= S.y.theta.b)/length(S)
      
      
      
      if(do.SA){
        # Use Stochastic Approximation to adjust m
        # so that there are about level.target*100% of the F.S.y.theta.b smaller or equal
        # to level.target

        Z <- (( F.S.y.theta.b <= level.target) - level.target)/SA.s
        sa.b[[k]][SA.i] <- as.numeric(F.S.y.theta.b <= level.target)
        F_mat[SA.i,k] <- F.S.y.theta.b
          
        if (SA.i %% 100 == 0){
          c <- c*decay
        }
        
        g <- c/(SA.i)
        

        m <- m0 + g*Z

        
        if(m> upper_scale*n) m <- upper_scale*n else if(m< n*lower_scale) m <- n*lower_scale
        m0 <- m

        m <- as.integer(floor(m) + (runif(1) <= m - floor(m)))
        
        log.m <- log(m)
        
        sa.a[[k]][SA.i] <- m/n
        level.m$m[k] <- m
        level.m$m0[k] <- m0
        level.m$log.m[k] <- log.m
        level.m$weighted[k] <- (1-weight)*level.m$weighted[k] + weight * log.m
        

        SA.b <- sa.b[[k]]
        SA.a <- sa.a[[k]]
        pl <- sa.pl[[k]]
        
        
      }
      
    }
  }
  
  
  list(level.m=level.m, SA.a =SA.a, SA.b = SA.b, loss_mat = loss_mat, result_mat = result_mat,
       F_mat = F_mat, beta_list = beta_list)
}



