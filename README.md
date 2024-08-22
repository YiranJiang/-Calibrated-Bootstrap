
Finite Sample Valid Inference via Calibrated Bootstrap (CB)
===============================================================
While widely used as a general method for uncertainty quantification, the bootstrap method encounters difficulties that raise concerns about its validity in practical applications. This paper introduces a new resampling-based method, termed *calibrated bootstrap*, designed to generate finite sample-valid parametric inference from a sample of size $n$. The central idea is to calibrate an $m\text{-out-of-}n$ resampling scheme, where the calibration parameter $m$ is determined against inferential pivotal quantities derived from the cumulative distribution functions of loss functions in parameter estimation. The method comprises two algorithms. The first, named *resampling approximation* (RA), employs a *stochastic approximation* algorithm to find the value of the calibration parameter $m=m_\alpha$ for a given $\alpha$ in a manner that ensures the resulting $m\text{-out-of-}n$ bootstrapped $1-\alpha$ confidence set is valid. The second algorithm, termed *distributional resampling* (DR),  is developed to further select samples of bootstrapped estimates from the RA step when constructing $1-\alpha$ confidence sets for a range of $\alpha$ values is of interest. The proposed method is illustrated and compared to existing methods using linear regression with and without $L_1$ penalty, within the context of a high-dimensional setting and a real-world data application. The paper concludes with remarks on a few open problems worthy of consideration.

## Related Publication

Yiran Jiang, Chuanhai Liu and Heping Zhang [Finite Sample Valid Inference via Calibrated Bootstrap.] *Major revision at Journal of the Royal Statistical Society Series B*.

<br>

### Reproduce Experimental Results in the Paper:
 
 The paper contains four main experiments:

1. Joint inference in the high-dimensional linear regression example
2. Marginal inference in the high-dimensional linear regression example
3. Inference in the high-dimensional linear regression with Lasso
4. Inference in the diabetes study

<br>

#### Joint inference in the high-dimensional linear regression example:

Required package:  

NA

```{R}
cd lr_example
Rscript ./lr_joint
Rscript ./lr_joint_summary.R
```

<br>


#### Marginal inference in the high-dimensional linear regression example

Required package:  

NA

```{R}
cd lr_example
Rscript ./lr_marginal.R
```

<br>


#### Inference in the high-dimensional linear regression with Lasso

Required package:  

glmnet


Run Simulation Experiments:

```{R}
cd lasso 
Rscript inference-n-100.R $id
Rscript inference-n-200.R $id
Rscript inference-n-500.R $id

```

 - id: experiment ID, can take arbitrary positive integer values

Alternative -- SLURM Job Script (modify the file if required):

```{sh}
sbatch myjob-n-100.sh
sbatch myjob-n-200.sh
sbatch myjob-n-500.sh
```

Summary (modify the directory in the file if required):

```{R}
Rscript summary-cb.R
```


Run Simulation Experiments with Standard Bootstrap and Residual Bootstrap:

```{R}
cd lasso 
Rscript standard-bootstrap-comparison.R $n $id

```
 - n: sample size, can take values 100, 200 or 500
 - id: experiment ID, can take arbitrary positive integer values

Alternative -- SLURM Job Script (modify the file if required):

```{sh}
sbatch myjob-bootstrap.sh
```

Summary (modify the directory in the  file if required):

```{R}
Rscript summary-bootstrap.R 
```

<br>


#### Inference in the diabetes study

Required Packages:

selectiveInference, glmnet, lars, tidyr, ggplot2, dplyr

```{R}
cd diabetes
Rscript inference.R
Rscript summary.R
```


## References

Friedman, Jerome, Hastie, Trevor, and Tibshirani, Rob. "Regularization Paths for Generalized Linear Models via Coordinate Descent." *Journal of Statistical Software*, 33.1 (2010): 1-22.

Efron, Bradley, Hastie, Trevor, Johnstone, Iain, and Tibshirani, Robert. "Least Angle Regression." *The Annals of Statistics*, 32.2 (2004): 407-499.

Lee, Jason D., Sun, Dennis L., Sun, Yuekai, and Taylor, Jonathan E. "Exact Post-Selection Inference, with Application to the Lasso." *The Annals of Statistics*, 44.3 (2016): 907-927.
