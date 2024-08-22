simulation.n <- 100

bootstrap.sd <- function(vec, B = 10000, this.seed = 123){
  set.seed(this.seed)
  this.means <- rep(NA, B)
  for (b in 1:B){
    new.vec <- vec[sample(seq(length(vec)),length(vec),replace = T)]
    this.means[b] <- mean(new.vec)
  }
  return(sd(this.means))
}

level.seq <- seq(0.05,0.95, by = 0.05)

all.mat.joint.bt <- matrix(nrow = 0, ncol = length(level.seq))
all.mat.marginal.bt <- matrix(nrow = 0, ncol = length(level.seq))
all.mat.joint.rb <- matrix(nrow = 0, ncol = length(level.seq))
all.mat.marginal.rb <- matrix(nrow = 0, ncol = length(level.seq))


for (i in 1:100){
  this.out <- readRDS(file = paste0("output-bootstrap/length-",i,"-n_",simulation.n,"-lasso.RDS"))
  
  all.mat.joint.bt <- rbind(all.mat.joint.bt,this.out$all.length.joint.bt/2)
  all.mat.marginal.bt <- rbind(all.mat.marginal.bt,this.out$all.length.marginal.bt)
  all.mat.joint.rb <- rbind(all.mat.joint.rb,this.out$all.length.joint.rb/2)
  all.mat.marginal.rb <- rbind(all.mat.marginal.rb,this.out$all.length.marginal.rb)
  
}


cat(round(mean(all.mat.joint.bt[,19]),3), round(bootstrap.sd(all.mat.joint.bt[,19]),3), "\n",
    round(mean(all.mat.joint.bt[,17]),3), round(bootstrap.sd(all.mat.joint.bt[,17]),3), "\n",
    round(mean(all.mat.joint.bt[,15]),3), round(bootstrap.sd(all.mat.joint.bt[,15]),3))

cat(round(mean(all.mat.marginal.bt[,19]),3), round(bootstrap.sd(all.mat.marginal.bt[,19]),3), "\n",
    round(mean(all.mat.marginal.bt[,17]),3), round(bootstrap.sd(all.mat.marginal.bt[,17]),3), "\n",
    round(mean(all.mat.marginal.bt[,15]),3), round(bootstrap.sd(all.mat.marginal.bt[,15]),3))

cat(round(mean(all.mat.joint.rb[,19]),3), round(bootstrap.sd(all.mat.joint.rb[,19]),3), "\n",
    round(mean(all.mat.joint.rb[,17]),3), round(bootstrap.sd(all.mat.joint.rb[,17]),3), "\n",
    round(mean(all.mat.joint.rb[,15]),3), round(bootstrap.sd(all.mat.joint.rb[,15]),3))

cat(round(mean(all.mat.marginal.rb[,19]),3), round(bootstrap.sd(all.mat.marginal.rb[,19]),3), "\n",
    round(mean(all.mat.marginal.rb[,17]),3), round(bootstrap.sd(all.mat.marginal.rb[,17]),3), "\n",
    round(mean(all.mat.marginal.rb[,15]),3), round(bootstrap.sd(all.mat.marginal.rb[,15]),3))


all.cover.mat <- matrix(0,nrow = 4, ncol = length(level.seq))
all.length.vec <- matrix(0,nrow = 4, ncol = length(level.seq))

for (i in 1:100){
  all.cover.mat <- all.cover.mat + readRDS(file = paste0("output-bootstrap/result-",i,"-n_,",simulation.n,"-lasso.RDS"))
}

all.cover.mat[,c(19,17,15)] / 500
