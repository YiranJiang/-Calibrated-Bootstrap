save.directory <- "./output"
simulation.n <- 100

N <- 100
level.seq <- seq(0.05,0.95,by = 0.05)
this.cover.mat <- matrix(0,nrow = 2, ncol = length(level.seq))
this.sum <- 0
this.id.vec <- c()
for (this.id in 1:100){
  
  if (file.exists(paste0(save.directory,"/result-",this.id,"-n_,",simulation.n,"-lasso.RDS"))){
    this.cover.mat <- this.cover.mat + readRDS(paste0(save.directory,"/result-",this.id,"-n_",simulation.n,"-lasso.RDS"))
    this.sum <- this.sum + 1
    this.id.vec <- c(this.id.vec, this.id)
  }

}

print(this.cover.mat/(this.sum*5)) ## Each node runs 5 repetitions
