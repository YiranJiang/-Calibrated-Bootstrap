interval.fixed <- function(m, xy,  M=c(1000,100), alpha, ...) {
  model <- attr(xy, "model")
  x <- xy[,-ncol(xy)]
  y <- xy[, ncol(xy)]
  x.matrix <- data.matrix(x)
  pen=0
  n <- nrow(xy)
  p <- ncol(xy)-1
  
  XX <- t(x.matrix)%*%x.matrix
  
  res.ML <- lm(model, data=xy)
  res.ML$sigmasq <- 1
  
  loss.ML <- loss(beta = res.ML$coefficients, res.ML$sigmasq, xy, pen=pen)
  
  SS <- list()
  SS_loss <- list()
  
    for(k in seq(1 + length(m))){
      if (k >= 2){
        cat("m:", m[k-1],"\n")
      }
      
    S_loss <- numeric(0)
    S <- numeric(0)
    for(i in 1:M[1]){
      if(i %%100 == 0){
        print(i)
      }
      if(k==1)
      {
        idx <- sample.int(n, size=n, replace=TRUE)
      } else {
        idx <- sample.int(n, size=m[k-1], replace=TRUE)
      }
      XY <- xy[idx,]
      res.b <- lm(model, data=XY)
      
      b <- res.b$coefficients - res.ML$coefficients
      res.b$sigmasq <- 1
      a <- sum(b*as.numeric(XX %*% b))/p
      
      S_loss[i] <-  loss(res.b$coefficients, res.b$sigmasq, xy, pen=pen)
      S[i] <- a
    }
    SS_loss <- append(SS_loss, list(S_loss))
    SS <- append(SS, list(S))
  }
  
  this.order_bt <- order(SS[1][[1]])
  
  nbins <- 10
  delta <- 1.0/nbins
  level.target <- seq(delta/2.0, 1.0, by = delta)
  
  count <- 0
  cat(" ","Truth ", "n-B"," CB", "\n")
  for (alpha in level.target){
    count <- count + 1
    
    this.index_bt <- (this.order_bt[floor(length(this.order_bt)*(1-alpha))])
    
    this.order_cb <- order(SS_loss[count + 1][[1]])
    this.index_cb <- (this.order_cb[floor(length(this.order_cb)*(1-alpha))])
    
    cat(alpha, (qchisq(1-alpha,df = p)/p), SS[1][[1]][this.index_bt], SS[count+1][[1]][this.index_cb],"\n")
  }
  
  
  
  return(list(SS = SS, SS_loss = SS_loss))
}

loss <- function(beta, sigmasq, xy, pen=0){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  n <- length(Y) 
  z <- Y - as.numeric(X%*%beta)
  s <- sqrt(sigmasq)
  sum(-dnorm(z,sd=s,log=TRUE)) - 0.5*pen*log(sigmasq) 
}


load("results2.RData")

this.out <- interval.fixed(out.1$level.m$m, xy, model,  M=c(1000, NULL))

save(this.out, file = "Joint_interval.RData")
  








