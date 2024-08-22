data.simu <- function(n=5000, kappa=0.5) {
  rnorm(n)
  p <- as.integer(floor(kappa*n))
  X <- matrix(rnorm(n*p), ncol=p, nrow=n)
  beta <- rep(0, p)
  Y <- X%*% beta + rnorm(n)
  data <-  data.frame(cbind(X, Y))
  names(data) <- c(paste("X", 1:p, sep="."), "Y")
  xnam <- paste0("X.", 1:p)
  (fmla <- as.formula(paste("Y ~ 0 + ", paste(xnam, collapse= "+"))))
  print(fmla)
  attr(data, "model") <- fmla
  data
}

my.lm <- function(formula, data, ...){
  res <- lm(formula, data, ...)
  res$beta <- res$coefficients
  res$sd <- sqrt(sum(res$residuals^2)/length(res$residuals))
  res
}

max.x <- 0.5

loss <- function(beta, sigmasq, xy, pen=0){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  n <- length(Y) 
  z <- Y - as.numeric(X%*%beta)
  s <- sqrt(sigmasq)
  sum(-dnorm(z,sd=s,log=TRUE)) - 0.5*pen*log(sigmasq) 
}

plot.level <- function(xy, max.x, level.target, ...){

  plot(c(0, max.x), c(0,1), type="n", xlab="t", ylab="Probability", ...)
  abline(h=level.target, col="gray")
}

SIR.sample <- function(xy, res.b, res.ML, size){
  X <- data.matrix(xy[,-ncol(xy)])
  Y <- xy$Y
  n <- length(Y) 

  z <- Y - as.numeric(X%*%res.ML$coefficients)
  s <- sqrt(res.ML$sigmasq)
  lnDensity.ML <- dnorm(z,sd=s,log=TRUE)

  z <- Y - as.numeric(X%*%res.b$coefficients)
  s <- sqrt(res.b$sigmasq)
  lnDensity.b <- dnorm(z,sd=s,log=TRUE)

  weights <- exp(lnDensity.b - lnDensity.ML)
  weights <- weights/sum(weights)
  idx <- sample(1:n, size=size, replace=TRUE, prob = weights)

  structure(xy[idx,], weights=weights)
}

# q = 0: Standard Bootstrap
my.sample <- function(x, size, replace=TRUE, prob=NULL, q=0){
  y <- runif(size)<q
  z <- sample(x, size=size-sum(y), replace = TRUE)
  c(x[y], z)
}



distributional.fixed <- function(m, xy, model, res.ML, loss.ML, XX, sir, x.matrix, add=TRUE, M=c(1000,100), ...) {
  x <- xy[,-ncol(xy), drop=FALSE]
pen=0
  n <- nrow(xy)
  p <- ncol(xy)-1
  F.list <- list()
  for(b in 1:M[1]){
    idx <- sample.int(n, size=sample(as.numeric(m), size=1), replace=TRUE)
    XY <- xy[idx,]
    res.b <- lm(model, data=XY)
    res.b$sigmasq <- 1 

    loss.B <- loss(res.b$coefficients, res.b$sigmasq, xy, pen=pen)
    S.y.theta.b <- -(loss.B-loss.ML)

    S <- numeric(0)
    for(db in 1:M[2]){
      while(TRUE){
        cat("\r", db)
    #simulated data: XY.db
        xy.simu <- if(sir) {
          SIR.sample(xy, res.b, res.ML, size=n)
        } else  {
          Y.mean <- as.numeric(x.matrix %*% res.b$coefficients)
          cbind(x, Y=Y.mean + sqrt(res.b$sigmasq)*rnorm(n))
        }

        res.db <- lm(model, data=xy.simu)
      	res.db$sigmasq <- 1 
        loss.DB.ML <- loss(res.db$coefficients, res.db$sigmasq, xy.simu, pen=pen)
        loss.B <- loss(res.b$coefficients, res.b$sigmasq, xy.simu, pen=pen)
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

    beta <- res.b$coefficients - res.ML$coefficients
    a <- sum(beta*as.numeric(XX %*% beta))/p
    F.list <- append(F.list, list(list(value=F.S.y.theta.b, beta=res.b$coefficients,
				 S = S.y.theta.b, check=a, beta = res.b$coefficients)))
  }

  x <- unlist(lapply(F.list, function(x) x$S))
  y <- unlist(lapply(F.list, function(x) x$value))
  plot(x,y, ylim=c(0,1), main="(a)", xlab=expression(T[y, theta]), ylab=expression(F[theta](T[y,theta])))

  z <- lapply(F.list, function(x) c(x$value, x$check))
  print(z)
  print(unlist(z))
  z <- unlist(z)
  dim(z) <- c(2, length(z)/2)
  y <- numeric(0)
  K <- numeric(0)
  for(i in 1:M[1]){
    U <- runif(1)
    k <- which.min(abs(U-z[1,]))
    y[i] <- z[2,k]
    K[i] <- k
  }
  x <- qchisq(((1:ncol(z))-0.5)/ncol(z), df=p)/p
  plot(x, sort(y), main="(b)",
       xlab="Theoretical Quantile", ylab="Empirical Quantile")
  abline(a=0, b=1, col="gray")
  abline(v=qchisq(c(0.05,0.95), df=30)/30, lty=2)


  list(F.list=F.list, K=K)
}




boot.ra <- function(xy, M =c(1000L, 200L), sir = TRUE, m, level.target=0.05,
		      n = nrow(xy),
		      p = ncol(xy)-1,
		      pen= 0, #p,
		      level.fun = NULL, 
		      distributional.plot=FALSE
){
  m_0 <- m
  pen <- 0
  model <- attr(xy, "model")
  if(missing(m)) m <- n 

  x <- xy[,-ncol(xy)]
  y <- xy[, ncol(xy)]
  x.matrix <- data.matrix(x)

  max.x <- M[1]
  if(!distributional.plot)
    plot.level(xy, main="(a)",
	    max.x=max.x, level.target=level.target)

  do.SA <- TRUE
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

  res.ML <- lm(model, data=xy)
  res.ML$sigmasq <- 1 

  loss.ML <- loss(beta = res.ML$coefficients, res.ML$sigmasq, xy, pen=pen)

  if(distributional.plot){
    out <- distributional.fixed(m, xy, model, res.ML, loss.ML,
			  XX=t(x.matrix)%*%x.matrix,
			  sir, x.matrix, add=TRUE, M=M, col=4)
    return(out)
  }
  pl <- numeric(0)

  Counts <- rep(0, length(level.target))
  names(Counts) <- level.target
  N0 <- as.integer(M[1]/2)
  F.levels <- c((level.target[-1] + level.target[-length(level.target)])/2, 1.0)
  beta.list <- lapply(as.list(1:length(level.target)), function(x) 
		      matrix(nrow=0, ncol=p))

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
    

    while(TRUE){
      idx <- sample.int(n, size=m, replace=TRUE)
      while (length(unique(idx)) <= p){
        idx <- sample.int(n, size=m, replace=TRUE)
      }
      cat("length(idx):", length(idx), " length(unique(idx)):", length(unique(idx)),"\n")
      XY <- xy[idx,]


      res.b <- lm(model, data=XY)
      res.b$sigmasq <- 1 

      loss.B <- loss(res.b$coefficients, res.b$sigmasq, xy, pen=pen)
      S.y.theta.b <- -(loss.B-loss.ML)
      
      if(!is.na(S.y.theta.b)) break;
      stop("S.y.theta.b is NA")
    }
    
    cat("[",b,"] S.y.theta.b: ", S.y.theta.b,"\n")
    if(is.na(S.y.theta.b)) stop("S.y.theta.b is NA")


    S <- numeric(0)
    for(db in 1:M[2]){
      while(TRUE){

        xy.simu <- if(sir) {
          SIR.sample(xy, res.b, res.ML, size=n)
        } else  {
          Y.mean <- as.numeric(x.matrix %*% res.b$coefficients)
          cbind(x, Y=Y.mean + sqrt(res.b$sigmasq)*rnorm(n))
        }

        res.db <- lm(model, data=xy.simu)
      	res.db$sigmasq <- 1 
        loss.DB.ML <- loss(res.db$coefficients, res.db$sigmasq, xy.simu, pen=pen)
        loss.B <- loss(res.b$coefficients, res.b$sigmasq, xy.simu, pen=pen)
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

    if(b >= N0){
      which.bin <- sum(F.S.y.theta.b <= F.levels)
      if(which.bin==k){
        Counts[which.bin] <- Counts[which.bin] + 1
        beta.list[[k]] <- rbind(beta.list[[k]], res.b$coefficients)
      }
    }


    if(do.SA){ 
          # Use Stochastic Approximation to adjust m
          # so that there are about level.target*100% of the F.S.y.theta.b smaller or equal
          # to level.target

          Z <- (( F.S.y.theta.b <= level.target) - level.target)/SA.s
          sa.b[[k]][SA.i] <- as.numeric(F.S.y.theta.b <= level.target)
       
          g <- 2000/(SA.i)

          m <- m0 + g*Z
          m0 <- m
          if(m0>2*n) m0 <- 2*n else if(m0<max(p+1,n/2)) m0 <- max(p+1,n/2)

	  m <- as.integer(floor(m) + (runif(1) <= m - floor(m)))

	  if(m>2*n) m <- 2*n else if(m< max(p+1,n/2)) m <- max(p+1,n/2)
	  log.m <- log(m)
	  
          sa.a[[k]][SA.i] <- m/n
          level.m$m[k] <- m
          level.m$m0[k] <- m0
	  level.m$log.m[k] <- log.m
	  level.m$weighted[k] <- (1-weight)*level.m$weighted[k] + weight * log.m
        
      cat("log.m:", log.m," m:", m,"\033[0m\n")

      SA.b <- sa.b[[k]]
      SA.a <- sa.a[[k]]
      pl <- sa.pl[[k]]
      
      
      if(b%%100 == 0 && k==nrow(level.m)){
        save(sa.b, sa.a, sa.pl, file = "results1.RData")
      }

      if(b==M[1]){
        # lines(seq(0, max.x, length=M[1])[1:length(SA.a)]/2, SA.a, col=k)
        lines(seq(0, max.x, length=M[1])[1:length(SA.a)]/2, SA.a-1, col=k, lty=2)
        lines(seq(0, max.x, length=M[1])[1:length(SA.b)], cumsum(SA.b)/(1:length(SA.b)), col=k)
      }
    }

  }
  }



  XX <- t(x.matrix)%*%x.matrix
  s2 <- 1
  SS <- list()

  for(k in c(1, 2)){ 
    ## k = 1: n-bootstrap
    ## k = 2: m-out-of-n-bootstrap
    
    S <- numeric(0)
    for(i in 1:M[1]){
      if(k==1)
      {
        idx <- sample.int(n, size=n, replace=TRUE)
      } else {
        if(nrow(level.m)==1) {
          log.m <- level.m$weighted
	} else {
          log.m <- sample(level.m$weighted, size=1)
	}
	m <- exp(log.m)
	m <- as.integer(floor(m) + (runif(1) <= m - floor(m)))
	cat("i:", i, "m:", m,"\n")
        idx <- sample.int(n, size=m, replace=TRUE)
      }
      XY <- xy[idx,]
      res.b <- lm(model, data=XY)

      b <- res.b$coefficients - res.ML$coefficients
      a <- sum(b*as.numeric(XX %*% b))/p
      S[i] <- a/s2
    }
    
    SS <- append(SS, list(S))
  }

  for(k in c(1, 2)){
    
    S <- SS[[k]]
    if(k==1){
      x <- qchisq(((1:length(S))-0.5)/length(S), df=p)/p
      y <- sort(S)
      plot(x, y, ylim=range(x,y, SS[[2]]),
	  xlab="Theoretical Quantile", ylab="Empirical Quantile",
          main="(b)")
    } else {
      points(qchisq(((1:length(S))-0.5)/length(S), df=p)/p, sort(S), pch=2)

    }
    abline(a=0, b=1, col=2)
    abline(v = qchisq(level.m$pl, df=p)/p, col="gray")
    SS <- append(SS, list(S))
  }


  list(level.m=level.m, SS=SS)
}

create.pdf <- T

# Resampling Approximation: find m for prespecified confidence levels
if(T){
  if(create.pdf) pdf(file="lr_joint_ra.pdf", width=12, height=7)
    
  par(mfcol=c(1,2), ask=FALSE)
  
  kappa <- 0.3; n <- 500
  
  max.x <- sqrt(kappa)*5.0
  
  xy <- data.simu(n, kappa)
  
  nbins <- 10
  delta <- 1.0/nbins
  level.target <- seq(delta/2.0, 1.0, by = delta)
  
  out.1 <- boot.ra (xy, sir=FALSE, m=as.integer(n*1.2), M=c(2000,10),
  	      level.target=level.target)
  save(out.1, xy, file = "results2.RData")
  
  
  if(create.pdf) dev.off()
}


## Distributional Resampling
if(T){
  if(create.pdf) pdf(file="lr_joint_dr.pdf", width=12, height=7)
  par(mfcol=c(1,2), ask=FALSE)

  nbins <- 10
  delta <- 1.0/nbins
  level.target <- seq(delta/2.0, 1.0, by = delta)

  out.2 <- boot.ra (xy, sir=FALSE, m=out.1$level.m$m, 
                      M=c(1000,500), ## For illustrative purpose, a large number of repetitions is used here. In practice, this is not required
		      level.target=level.target, distributional.plot=TRUE)
  save(out.2, file = "results3.RData")

  if(create.pdf) dev.off()
}


