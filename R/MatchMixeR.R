## numerically robust matrix inverse
rsolve <- function(mat, d.prop=1e-4, dmin=1e-6, dmax=1e6){
  o <- svd(mat)
  ## singular value threshold
  thresh <- min(max(sum(o$d) * d.prop, dmin), dmax)
  inv.mat <- o$v %*% diag(1/(pmax(o$d, thresh))) %*% t(o$u)
  return(inv.mat)
}

# Per-gene OLS model
#'@param X gene expression matrix
#'@param Y gene expression matrix
## this is the robust version of OLS
OLS <- function(X, Y, min.var=1e-3){
  n <- ncol(X)
  ## A fast row-wise univariate OLS regression based on linear algebra
  Xbar <- rowMeans(X, na.rm=TRUE); Ybar <- rowMeans(Y, na.rm = TRUE)
  X.c <- sweep(X, 1, Xbar); Y.c <- sweep(Y, 1, Ybar)
  CovXY <- rowSums(X.c * Y.c, na.rm = TRUE)/(n-1)
  negCor <- which(CovXY<0)
  VarX <- rowSums(X.c^2, na.rm = TRUE)/(n-1)
  VarY <- rowSums(Y.c^2, na.rm = TRUE)/(n-1)
  ## Pearson corr. coefs. It is OK to have N/As
  cc <- CovXY / sqrt(VarX * VarY)
  ## robust variances
  varXmin <- max(quantile(VarX, 0.05)/10, min.var)
  lowXvar <- which(VarX<varXmin)
  VarX[lowXvar] <- varXmin
  if (length(lowXvar) > nrow(X)/5) warning(paste0("There are many (n=", length(lowXvar), ") genes with very low sample variance on platform X. Please consider running a nonspecific filtering before applying MatchMixerR."))
  ## CovXY[lowXvar] <- varXmin #so that beta1=1 when VarX is too small
  ## regression coefs.
  beta1 <- CovXY / VarX
  beta0 <- Ybar - beta1 * Xbar
  betamat <- cbind("Intercept"=beta0, "Slope"=beta1)
  ## predictions
  Yhat <- sweep(X, 1, beta1, FUN="*") + beta0 %o% rep(1, ncol(X))
  ## RSS vector
  RSS <- rowMeans((Y - Yhat)^2, na.rm=TRUE)
  SS1 <- rowMeans(Yhat^2, na.rm=TRUE)
  SST <- rowMeans(Y^2, na.rm=TRUE)
  var.prop <- SS1/SST; var.prop[SS1==0] <- 0
  ## trouble maker genes
  problematic.genes <- list(lowXvar=lowXvar, negCor=negCor)
  return(list(betamat=betamat, corr=cc, covxy=CovXY, varx=VarX, Yhat=Yhat, RSS=RSS, var.explained.prop=var.prop, Xbar=Xbar, Ybar=Ybar, problematic.genes=problematic.genes))
}

## The main engine of cross-platform normalization procedure: FLMER is based on
## the proposed moment-based method
#'@param X gene expression matrix
#'@param Y gene expression matrix
#'@export
flmer <- function(Xmat, Ymat, min.var=1e-3){
  Xmat <- as.matrix(Xmat); Ymat <- as.matrix(Ymat)
  ## 11/12/2018. Use new notations; p==ngenes; n==sample size
  p <- nrow(Xmat); n <- ncol(Xmat); N <- p*n
  ## calculate some useful items
  Xibar <- rowMeans(Xmat,na.rm=TRUE); Xbar <- mean(Xibar)
  Xmat.c <- Xmat - Xibar
  Yibar <- rowMeans(Ymat,na.rm=TRUE); Ybar <- mean(Yibar)
  Ymat.c <- Ymat - Yibar
  covXY <- rowSums(Xmat.c * Ymat.c, na.rm=TRUE)
  VarX <- rowSums(Xmat.c * Xmat.c, na.rm=TRUE)
  ## robust variances
  varXmin <- max(quantile(VarX, 0.05)/10, min.var)
  lowXvar <- which(VarX<varXmin)
  VarX[lowXvar] <- varXmin
  beta.yixi <- covXY / VarX
  ## The overall regression
  Xc <- Xmat - mean(Xmat, na.rm=TRUE); Yc <- Ymat - mean(Ymat, na.rm=TRUE)
  VarX.overall <- max(sum(Xc^2, na.rm=TRUE), varXmin)
  beta.yx <- sum(Xc*Yc, na.rm=TRUE) / VarX.overall
  ## In this case, we assume  that there is no collinearity in Z
  qprime <- 2*p
  var.epsilon <- sum((Ymat.c - Xmat.c * beta.yixi)^2, na.rm=TRUE) / (N - qprime)
  ## just to ensure that var.epsilon != 0
  var.epsilon <- max(var.epsilon, min.var)
  ## Preparing for calculating S, \| \mathbf{R}_{Z|\mathbf{X}} \|^2,
  ## and \|Z' \mathbf{R}_{Z|\mathbf{X}} \|^2 terms.
  XiYibar <- covXY/n + Xibar*Yibar
  ## Xi2bar = \|Xi\|^2 / n
  Xi2bar <- VarX/n + Xibar^2
  ## Xs is Xmat standardized by the global mean/sd
  Xs <- (Xmat - Xbar) / sqrt(VarX.overall)
  Xsbar <- rowMeans(Xs, na.rm=TRUE)
  ## XsX is the inner producd between Xs and X. It is denoted as
  ## \varsigma in the notes.
  XsX <- rowSums(Xs*Xmat, na.rm=TRUE)
  XsNorm2 <- sum(Xs^2, na.rm=TRUE)
  ## the normalized version of S
  S <- n^2 * (sum((Yibar - Ybar -(Xibar - Xbar)*beta.yx)^2) + sum((XiYibar -Ybar*Xibar -(Xi2bar - Xbar*Xibar)*beta.yx)^2)) / var.epsilon
  ## Calculate \| \mathbf{R}_{Z|\mathbf{X}} \|^2
  Projxz.norm2 <- n + XsNorm2 * ( n^2*sum(Xsbar^2) + sum(XsX^2)) + n*sum(Xibar^2)/p
  Rzx.norm2 <- sum(Xmat^2, na.rm=TRUE) + N -Projxz.norm2
  ## General terms; without n^4 yet.
  S1 <- sum(Xibar^2); S2 <- sum(Xsbar^2); S3 <- sum(XsX^2)/n^2
  ZPZ.norm2 <- p^2/N^2 + S2^2 + 2*( p*S1/N^2 + S2*S3) + S1^2/N^2 + 2*(sum(Xibar*XsX)^2)/(N*n^2) + S3^2
  ## subtractions; to multiply by *n^4
  Minus.terms <- sum((1/N + Xsbar^2)^2) + 2*sum((Xibar/N + Xsbar*XsX/n)^2) + sum((Xibar^2/N + XsX^2/n^2)^2)
  ## added terms
  XiNorm2 <- rowSums(Xmat^2, na.rm=TRUE)
  Add.terms <- sum((1/n -1/N -Xsbar^2)^2) + 2*sum((Xibar/n -Xibar/N -Xsbar*XsX/n)^2) + sum((XiNorm2/n^2 -Xibar^2/N -XsX^2/n^2)^2)
  ZRzx.norm2 <- n^4*(ZPZ.norm2 - Minus.terms + Add.terms)
  ## Estimate lambda based on Moment matching
  lambdahat <- max(0, (S - Rzx.norm2) / ZRzx.norm2)
  ## Inference of the fixed effects. First, compute some small
  ## "building blocks".
  Axx <- XiNorm2 - (lambdahat/(1+n*lambdahat)) * n^2 * Xibar^2
  A1x <- n*Xibar/(1+n*lambdahat)
  A1y <- n*Yibar/(1+n*lambdahat)
  Axy <- XiYibar*n - (lambdahat/(1+n*lambdahat)) * n^2 * Xibar* Yibar
  W11 <- n/(1+n*lambdahat) - lambdahat*A1x^2/(1+lambdahat*Axx)
  W1x <- A1x/(1+lambdahat*Axx)
  W1y <- A1y -lambdahat*A1x*Axy / (1+lambdahat*Axx)
  Wxx <- Axx/(1+lambdahat*Axx)
  Wxy <- Axy / (1+lambdahat*Axx)
  ## Now the actual estimation
  covBeta <- var.epsilon * rsolve(matrix(c(sum(W11), sum(W1x), sum(W1x), sum(Wxx)), 2))
  betahat <- drop(covBeta %*% c(sum(W1y), sum(Wxy)) / var.epsilon)
  names(betahat) <- c("Intercept", "Slope")
  dimnames(covBeta) <- list(c("Intercept", "Slope"), c("Intercept", "Slope"))
  ## Mixed effects terms via EBLUP
  gamma0hat <- lambdahat * (W1y - betahat[1]*W11 - betahat[2]*W1x)
  gamma1hat <- lambdahat * (Wxy - betahat[1]*W1x - betahat[2]*Wxx)
  ## Individual betas
  betamat <- t(rbind(gamma0hat, gamma1hat) + betahat)
  colnames(betamat) <- c("Intercept", "Slope")
  ## t-statistics for the fixed effects
  t.fixed <- betahat/sqrt(pmax(diag(covBeta), min.var))
  ## the predicted Yhat
  Yhat <- Xmat * betamat[, "Slope"] + betamat[, "Intercept"]
  return(list(betahat=betahat, betamat=betamat, Yhat=Yhat,
              lambdahat=lambdahat, var.epsilon=var.epsilon,
              covBeta=covBeta, t.fixed=t.fixed))
}

## Implement the covariance transformed FLMER as follows.
#'Match-MixeR
#'
#'fit mixed effect regression model
#'@param X gene expression matrix
#'@param Y gene expression matrix
#'@export
MM <- function(Xmat, Ymat){
  rr1 <- OLS(Xmat, Ymat)
  cov1 <- cov(rr1$betamat)
  s0 <- sqrt(cov1[1,1]); s1 <- sqrt(cov1[2,2]); rho <- cov1[1,2]/s0/s1
  a12 <- rho/sqrt(1-rho^2); a22 <- s1/(s0 * sqrt(1-rho^2))
  A <- matrix(c(1, a12,
                0, a22), 2, byrow=TRUE)
  ## the X transformation
  Xtilde <- a12 + a22*Xmat
  ## apply flmer() to the covariance transformed data
  rr3 <- flmer(Xtilde, Ymat)
  ## the reverse transformation. Note that each *row* of "betamat" is
  ## beta_i; so we have to transpose the reverse matrix multiplication.
  betamat <- rr3$betamat %*% t(A); colnames(betamat) <- colnames(rr3$betamat)
  ## other misc. items
  betahat <- drop(A %*% rr3$betahat); names(betahat) <- names(rr3$betahat)
  lambdahat <- rr3$lambdahat
  covGamma <-  rr3$lambdahat * rr3$var.epsilon * (A %*% t(A))
  Yhat <- rr3$Yhat
  var.epsilon <- rr3$var.epsilon
  covBeta <- A %*% rr3$covBeta %*% t(A)
  t.fixed <- betahat/sqrt(diag(covBeta))
  return(list(betahat=betahat, betamat=betamat, Yhat=Yhat,
              lambdahat=lambdahat, var.epsilon=var.epsilon,
              var.explained.prop=rr1$var.explained.prop,
              covGamma=covGamma, covBeta=covBeta, t.fixed=t.fixed,
              varx=rr1$VarX, problematic.genes=rr1$problematic.genes))
}


