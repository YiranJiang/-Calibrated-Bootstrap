install_and_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      library(package, character.only = TRUE)
    }
  }
}

packages <- c("selectiveInference", "lars","tidyr","ggplot2","dplyr","glmnet")
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

res.ML <- lasso_model(xy, pen)

res.ML$coefficients <- coef(res.ML)[,1]
res.ML$coefficients <- res.ML$coefficients[2:(p+1)]
res.ML$coefficients

alpha <- 0.05


## Calculating Exact POSI Intervals
this.posi <- fixedLassoInf(X,y,coef(res.ML,s = pen,extract = TRUE)[-1],pen*n,alpha = alpha)
this.names <- colnames(X)
this.names[c(5,6,7,8,9,10)] <- c("S1","S2","S3","S4","S5","S6")
this.predictor.names <- this.names[this.posi$vars]



## Calculating Standard Bootstrap Intervals

sb.beta.matrix <- matrix(nrow = 1000, ncol = p)

for(i in 1:1000){
  idx <- sample.int(n, size=n, replace=TRUE)
  
  XY <- xy[idx,]

  res.b <- lasso_model(XY, pen)
  res.b$coefficients <- coef(res.b)[,1]
  res.b$coefficients <- res.b$coefficients[2:(p+1)]

  sb.beta.matrix[i,] <- res.b$coefficients
}


sb.full.ci.matrix <- matrix(nrow = p, ncol = 2)
for (i in 1:p){
  sb.full.ci.matrix[i,] <- quantile(sb.beta.matrix[,i], c(alpha/2, 1-alpha/2))
}


cb.full.ci.matrix <- readRDS("cb-diabetes-result.RDS")

ci_matrix1 <- cb.full.ci.matrix[this.posi$vars,]
ci_matrix2 <- sb.full.ci.matrix[this.posi$vars,]
ci_matrix3 <- this.posi$ci

methods <- c("Calibrated Bootstrap","Standard Bootstrap", "Exact POSI")
ci_data1 <- data.frame(method = methods[1], lower = ci_matrix1[, 1], upper = ci_matrix1[, 2], predictor = this.predictor.names)
ci_data2 <- data.frame(method = methods[2], lower = ci_matrix2[, 1], upper = ci_matrix2[, 2], predictor = this.predictor.names)
ci_data3 <- data.frame(method = methods[3], lower = ci_matrix3[, 1], upper = ci_matrix3[, 2], predictor = this.predictor.names)

combined_data <- rbind(ci_data1, ci_data2, ci_data3)
combined_data$predictor <- factor(combined_data$predictor, levels = this.predictor.names)
combined_data$method <- factor(combined_data$method, levels = methods)
combined_data <- combined_data %>% arrange(method)


## Obtaining the visualized results

pdf(file = "diabetes_plot.pdf")

this.plot <- ggplot(combined_data, aes(x = predictor, ymin = lower, ymax = upper, color = method)) +
  geom_errorbar(position = position_dodge(width = 0.4), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  coord_flip() +
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.85, 0.55),  # Adjust these values as needed
        legend.background = element_rect(color = "black", fill = "white"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))


this.plot
dev.off()

save.image(file = "diabetes.RData")









