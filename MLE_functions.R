# BiTEN Project: MLE estimation functions
# Frank Marrs
# 09/13/16
#
# This file contains functions supporting MLE estimation of A and B matrices on trade data based on inputs
# 



require(Matrix)   # ?
require(countrycode)   # countrycode
require(fields)   # image.plot
require(glmnet)  # cv.glmnet
require(MASS)  # ginv
require(abind)  # abind
require(foreach)
require(doMC)
require(foreign)  # read.dta
require(arcdiagram)  # arcplot



##########################
### Fitting functions  ###
##########################

# These functions perform various linear model fits based on inputs


# HOFF multiplicative bilinear model: 
# Fit via MLE block coordinate descent 
# D is matrix of lagged Y (covariates)
# Y is oservations
# X is covariates
# for multiple time periods
# Returns A, B, beta, betaOLS, Yhat, 2LL, 2LLinit
fit_MLE_array_bilinear <- function(Y, D, X, verbose=T, printout=10, tol=1e-6, init="I", sigma_init=1, use_cov=T, imax=1e6)
{
  # Get sizes and check
  if(dim(Y)[3] != dim(D)[3]){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  if(length(dim(Y)) < 3){ stop("Y is not an array")}
  
  # Set NAs to zero to start
  i_na <- is.na(Y)   # indices of NAs
  Y[i_na] <- 0
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  m <- dim(D)[1]
  n <- dim(D)[2]
  tmax <- dim(Y)[3]
  
  # m <- rankA 
  # k <- rankB
  # nper <- length(Y)/(2*(S*k + L*m) + dim(X)[4])    # Warning for data size
  
  # Initialize
  if(strtrim(init,3) == "ran") {
    if(verbose == T){
      cat("Randomly initializing A and B \n")
    }
    A <- matrix(rnorm(S*m, 0, sigma_init), S, m)
    B <- matrix(rnorm(L*n, 0, sigma_init), L, n)
    
  } else if (strtrim(init, 1) == "I") {
    if(verbose == T){
      cat("Initializing A,B with identity \n\n")
    }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
    
    A <- diag(S)
    B <- diag(L)
    
  } else if (strtrim(init, 3) == "reg") {
    if(verbose == T){
      cat("Initializing A,B with regression \n\n")
    }
    if(S>L){ first_fit <- "A"} else {first_fit <- "B"}
    # out <- initialize_mats(D, Y, 1, n, first=first_fit)   # initialize with n-rank regression
    
    A <- diag(S)   # update step will update this one first, so it's really immaterial
    # B <- tcrossprod(out$W,out$Z)
    B <- init_bilinear(Y,D)
    
  } else { stop("Invalid initialization type")}
  
  Y[i_na] <- ( amprod( amprod(D, A, 1), t(B), 2) )[i_na]   # initialize NAs 
  
  # Initialization
  Ainit <- A
  Binit <- B
  if(use_cov){
    betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
  } else { betaOLS <- NA}
  
  
  # Find optimal values
  change <- 100
  count <- 0
  # SSE <- 1
  while(change > tol & count < imax){
    
    # Update for beta
    Ystar <- Y - amprod( amprod(D, A, 1), t(B), 2) 
    if(use_cov){
      beta <- matrix(lm(c(Ystar) ~ -1 + t(mat(X, 4)) )$coef, nrow=1)
      Xbeta <- drop(amprod(X, beta, 4))
    } else {
      beta <- NA
      Xbeta <- 0
    }
    
    
    if(count == 0) {   # save initial LL and beta
      beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( (Ystar - Xbeta )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Update A and B
    Ytilde <- Y - Xbeta
    result <- update_MLE_bilinear(D, Ytilde, A, B)
    Anew <- result$A
    Bnew <- result$B
    changeAB <- max(abs(c(c(A - Anew), c(B-Bnew))))  # doesn't make much difference stopping A,B vs using LL
    A <- Anew
    B <- Bnew
    
    Yhat <- Y - Ytilde + amprod(amprod(D, A, 1), t(B), 2 )  # estimate
    Y[i_na] <- Yhat[i_na]  # update NAs
    
    LLnew <- -length(Y)*log( sum( (Ytilde - Yhat )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    # SSEnew <- LLnew <- mean( (Y-Yhat)^2,na.rm=TRUE )
    # change <- (SSEnew - SSE)/SSE
    # SSE <- SSEnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T) {
    cat("\n************************************************ \n")
    
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  # Prediction
  Yhat <- Xbeta + amprod( amprod(D, A, 1), t(B), 2) 
  
  return(list(A=A, B=B, Yhat=Yhat, beta= beta, betaOLS = betaOLS, LL2 = LL, LL2_init=LLinit, iter=count))
}


# Fit via MLE block coordinate descent 
# asymmetric A and B
# D is matrix of lagged Y (covariates)
# Y is oservations
# X is covariates
# for multiple time periods
# Returns A, B, beta, betaOLS, Yhat, 2LL, 2LLinit
fit_MLE_array <- function(Y, D, X, rankA=1, rankB=1, verbose=T, printout=10, tol=1e-9, init="random", sigma_init=1, use_cov=T)
{
  # Get sizes and check
  if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  if(length(dim(Y)) < 3){ stop("Y is not an array")}
  
  i_na <- is.na(Y)  # save NAs
  Y[i_na] <- 0  # set to zero
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  tmax <- dim(Y)[3]
  
  k <- rankA 
  m <- rankB
  nper <- length(Y)/(2*(S*k + L*m) + dim(X)[4])    # Warning for data size
  
  # Initialize
  if(strtrim(init,3) == "ran") {
    if(verbose == T){
      cat("Randomly initializing U,V,W,Z \n")
    }
    U <- matrix(rnorm(S*k, 0, sigma_init), S, k)
    V <- matrix(rnorm(S*k, 0, sigma_init), S, k)
    W <- matrix(rnorm(L*m, 0, sigma_init), L, m)
    Z <- matrix(rnorm(L*m, 0, sigma_init), L, m)
    
    A <- tcrossprod(U,V)
    BT <- tcrossprod(Z,W)
    
  } else if (strtrim(init, 3) == "reg") {
    if(verbose == T){
      cat("Initializing UVWZ with regression \n\n")
    }
    out <- initialize_mats(D, Y, k, m, first="B")
    
    U <- out$U
    V <- out$V
    W <- out$W
    Z <- out$Z
    
    A <- tcrossprod(U,V)
    BT <- tcrossprod(Z,W)
    
  } else if (strtrim(init, 1) == "I") {
    if(verbose == T){
      cat("Initializing A, B as identity \n\n")
    }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
    
    A <- diag(S)
    BT <- diag(L)
    
  } else { stop("Invalid initialization type")}
  
  
  
  Ainit <- A
  BTinit <- BT
  
  # Preliminary matrices to avoid recalculating
  DDT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x])))    # t(D) %*% D
  DTD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x])))    # D %*% t(D)
  
  if(use_cov){
    betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
  } else { betaOLS <- NA}
  
  # Find optimal values
  change <- 100
  count <- 0
  # if(nper < 5){ warning(paste0("Only ", round(nper, 3), " data points per coefficient"))}
  
  
  while(change > tol){
    
    # Update for beta
    Ystar <- Y - amprod(D, A, 1) - amprod(D, BT, 2)
    if(use_cov){
      beta <- matrix(lm(c(Ystar) ~ -1 + t(mat(X, 4)) )$coef, nrow=1)
      Xbeta <- drop(amprod(X, beta, 4))
    } else {
      beta <- NA
      Xbeta <- 0
    }
    
    if(count == 0) {   # save initial LL and beta
      beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( (Ystar - Xbeta )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Update A and B
    Ytilde <- Y - Xbeta
    DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], Ytilde[,,x])))    
    DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Ytilde[,,x])))    
    result <- update_MLE_asymmetric(D, U, V, W, Z, DDT, DTD, DYT, DTY)
    U <- result$U
    V <- result$V
    W <- result$W
    Z <- result$Z
    Anew <- tcrossprod(U,V)
    BTnew <- tcrossprod(Z,W)
    changeAB <- max(abs(c(c(A - Anew), c(BT-BTnew))))  # doesn't make much difference stopping A,B vs using LL
    A <- Anew
    BT <- BTnew
    
    Yhat <- Y - Ytilde + amprod(D, A, 1) + amprod(D, BT, 2)  # estimate
    Y[i_na] <- Yhat[i_na]  # update NAs
    
    # LLnew <- -sum( (Ytilde - (amprod(D, A, 1) + amprod(D, BT, 2)) )^2)
    LLnew <- -length(Y)*log( sum( (Y - Yhat )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T) {
    cat("\n************************************************ \n")
    
    # cat("True 2LL*tau^2: \t", LLtrue, "\n")
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  # Post-process
  A <- tcrossprod(U,V)
  B <- tcrossprod(W,Z)
  
  Yhat <- Xbeta + amprod(D, A, 1) + amprod(D, t(B), 2)
  
  return(list(A=A, B=B, Yhat = Yhat, beta= beta, betaOLS = betaOLS, LLt2 = LL, LLt2_init=LLinit))
}


# BLIN Additive model:
# Fit via MLE block coordinate descent 
# asymmetric A and B
# D is matrix of lagged Y (covariates)
# Y is oservations
# X is covariates
# for multiple time periods
# Returns A, B, beta, betaOLS, Yhat, 2LL, 2LLinit
fit_MLE_array_additive <- function(Y, D, X, type="blin", verbose=T, printout=10, tol=1e-6, init="random", sigma_init=1, use_cov=T, maxit=1e4)
{
  # Get sizes and check
  if(dim(Y)[3] != dim(D)[3]){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  if(length(dim(Y)) < 3){ stop("Y is not an array")}
  
  i_na <- is.na(Y)
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  m <- dim(D)[1]
  n <- dim(D)[2]
  tmax <- dim(Y)[3]
  
  
  # Initialize
  if(strtrim(init,3) == "ran") {
    if(verbose == T){
      cat("Randomly initializing A,B \n")
    }
    A <- matrix(rnorm(S*m, 0, sigma_init), S, m)
    B <- matrix(rnorm(L*n, 0, sigma_init), L, n)
    
  } else if (strtrim(init, 3) == "reg") {
    if(verbose == T){
      cat("Initializing A, B with regression \n\n")
    }
    if(S>L){ first_fit <- "A"} else {first_fit <- "B"}
    out <- initialize_mats(D, Y, 1, n, first=first_fit)   # initialize with n-rank regression
    
    U <- out$U
    V <- out$V
    W <- out$W
    Z <- out$Z
    
    A <- tcrossprod(U,V)
    B <- tcrossprod(W,Z)
    
  } else if (strtrim(init, 1) == "I") {
    if(verbose == T){
      cat("Initializing A, B as identity \n\n")
    }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
    
    A <- diag(S)
    B <- diag(L)
    
  } else { stop("Invalid initialization type")}
  
  
  Ainit <- A
  Binit <- B
  
  # Preliminary matrices to avoid recalculating
  if(type == "sadd"){
    Jl <- matrix(1,L,L)
    Js <- matrix(1,S,S)
    DDT <- L*Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, D[,,x] )))    # LDJD^T
    DTD <- S*Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% D[,,x])))    # SD^T J
  } else if (tolower(strtrim(type, 1)) == "b") {
    Jl <- diag(L)
    Js <- diag(S)
    DDT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x])))    # D %*% t(D)
    DTD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x])))    # t(D) %*% D
  } else { stop( "Invalid type " ) }
  
  
  Y[i_na] <- (amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2) )[i_na]   # initialize NAs
  
  
  if(use_cov){
    betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
  } else { betaOLS <- NA}  
  
  # Find optimal values
  change <- 100
  count <- 0
  # if(nper < 5){ warning(paste0("Only ", round(nper, 3), " data points per coefficient"))}
  
  
  while(change > tol & count < maxit){
    
    # Update for beta
    Ystar <- Y - amprod(amprod(D, A, 1), Jl, 2) - amprod(amprod(D, Js, 1), t(B), 2)
    if(use_cov){
      beta <- matrix(lm(c(Ystar) ~ -1 + t(mat(X, 4)) )$coef, nrow=1)
      Xbeta <- drop(amprod(X, beta, 4))
    } else {
      beta <- NA
      Xbeta <- 0
    }
    
    if(count == 0) {   # save initial LL and beta
      beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( (Ystar - Xbeta )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Update A and B
    Ytilde <- Y - Xbeta
    if(type == "sadd"){
      Jl <- matrix(1,L,L)
      Js <- matrix(1,S,S)
      DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, Ytilde[,,x])))    # D J Y^T
      DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% Ytilde[,,x])))    # D^T J Y
    } else if (tolower(strtrim(type, 1)) == "b") {
      Jl <- diag(L)
      Js <- diag(S)
      DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], Ytilde[,,x])))    # D %*% t(Y)
      DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Ytilde[,,x])))    # t(D) %*% Y
    } else { stop( "Invalid type " ) }
    result <- update_MLE_additive(A, B, D, DDT, DTD, DYT, DTY, type)
    Anew <- result$A
    Bnew <- result$B
    changeAB <- max(abs(c(c(A - Anew), c(B-Bnew))))  # doesn't make much difference stopping A,B vs using LL
    A <- Anew
    B <- Bnew
    
    Yhat <- Y - Ytilde + amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2) 
    Y[i_na] <- Yhat[i_na]   # update NAs
    
    LLnew <- -length(Y)*log( sum( (Y - Yhat)^2) ) - length(Y)    #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T) {
    cat("\n************************************************ \n")
    
    # cat("True 2LL*tau^2: \t", LLtrue, "\n")
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  Yhat <- Xbeta + amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2)
  
  return(list(A=A, B=B, Yhat=Yhat, beta= beta, betaOLS = betaOLS, LLt2 = LL, LLt2_init=LLinit))
}


# BLIN Additive model, multipartite version with possibly different lags:
# Fit via MLE block coordinate descent 
# asymmetric A and B
# D is matrix of lagged Y (covariates)
# Y is oservations
# X is covariates
# for multiple time periods
# Returns A, B, beta, betaOLS, Yhat, 2LL, 2LLinit
fit_MLE_array_additive_multi <- function(Y, D, X, type="blin", verbose=T, 
                                         printout=10, tol=1e-6, init="random", sigma_init=1, use_cov=T)
{
  # Get sizes and check
  if(dim(Y)[3] != dim(D)[3]){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  if(length(dim(Y)) < 3){ stop("Y is not an array")}
  
  i_na <- is.na(Y)
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  m <- dim(D)[1]
  n <- dim(D)[2]
  tmax <- dim(Y)[3]
  
  
  # Initialize
  if(strtrim(init,3) == "ran") {
    if(verbose == T){
      cat("Randomly initializing A,B \n")
    }
    A <- matrix(rnorm(S*m, 0, sigma_init), S, m)
    B <- matrix(rnorm(L*n, 0, sigma_init), L, n)
    
  } else if (strtrim(init, 3) == "reg") {
    if(verbose == T){
      cat("Initializing A, B with regression \n\n")
    }
    if(S>L){ first_fit <- "A"} else {first_fit <- "B"}
    out <- initialize_mats(D, Y, 1, n, first=first_fit)   # initialize with n-rank regression
    
    U <- out$U
    V <- out$V
    W <- out$W
    Z <- out$Z
    
    A <- tcrossprod(U,V)
    B <- tcrossprod(W,Z)
    
  } else if (strtrim(init, 1) == "I") {
    if(verbose == T){
      cat("Initializing A, B as identity \n\n")
    }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
    
    A <- diag(S)
    B <- diag(L)
    
  } else { stop("Invalid initialization type")}
  
  
  Ainit <- A
  Binit <- B
  
  # Preliminary matrices to avoid recalculating
  if(type == "sadd"){
    Jl <- matrix(1,L,L)
    Js <- matrix(1,S,S)
    DDT <- L*Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, D[,,x] )))    # LDJD^T
    DTD <- S*Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% D[,,x])))    # SD^T J
  } else if (tolower(strtrim(type, 1)) == "b") {
    Jl <- diag(L)
    Js <- diag(S)
    DDT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x])))    # D %*% t(D)
    DTD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x])))    # t(D) %*% D
  } else { stop( "Invalid type " ) }
  
  
  Y[i_na] <- (amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2) )[i_na]   # initialize NAs
  
  
  if(use_cov){
    betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
  } else { betaOLS <- NA}  
  
  # Find optimal values
  change <- 100
  count <- 0
  # if(nper < 5){ warning(paste0("Only ", round(nper, 3), " data points per coefficient"))}
  
  
  while(change > tol){
    
    # Update for beta
    Ystar <- Y - amprod(amprod(D, A, 1), Jl, 2) - amprod(amprod(D, Js, 1), t(B), 2)
    if(use_cov){
      beta <- matrix(lm(c(Ystar) ~ -1 + t(mat(X, 4)) )$coef, nrow=1)
      Xbeta <- drop(amprod(X, beta, 4))
    } else {
      beta <- NA
      Xbeta <- 0
    }
    
    if(count == 0) {   # save initial LL and beta
      beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( (Ystar - Xbeta )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Update A and B
    Ytilde <- Y - Xbeta
    if(type == "sadd"){
      Jl <- matrix(1,L,L)
      Js <- matrix(1,S,S)
      DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, Ytilde[,,x])))    # D J Y^T
      DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% Ytilde[,,x])))    # D^T J Y
    } else if (type == "biten") {
      Jl <- diag(L)
      Js <- diag(S)
      DYT <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], Ytilde[,,x])))    # D %*% t(Y)
      DTY <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Ytilde[,,x])))    # t(D) %*% Y
    } else { stop( "Invalid type " ) }
    result <- update_MLE_additive(A, B, D, DDT, DTD, DYT, DTY, type)
    Anew <- result$A
    Bnew <- result$B
    changeAB <- max(abs(c(c(A - Anew), c(B-Bnew))))  # doesn't make much difference stopping A,B vs using LL
    A <- Anew
    B <- Bnew
    
    Yhat <- Y - Ytilde + amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2) 
    Y[i_na] <- Yhat[i_na]   # update NAs
    
    LLnew <- -length(Y)*log( sum( (Y - Yhat)^2) ) - length(Y)    #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T) {
    cat("\n************************************************ \n")
    
    # cat("True 2LL*tau^2: \t", LLtrue, "\n")
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  Yhat <- Xbeta + amprod(amprod(D, A, 1), Jl, 2) + amprod(amprod(D, Js, 1), t(B), 2)
  
  return(list(A=A, B=B, Yhat=Yhat, beta= beta, betaOLS = betaOLS, LLt2 = LL, LLt2_init=LLinit))
}





# Additive model solution via regression for both BiTEN and smoothed additive
# Perform regression of Y onto D for a given type
# D is a 3-mode array of covariates
# Y is a 3-mode array of responses of the same size as D
# X is a 4-mode array of covariates
# type is "sadd" or "biten", type of regression to perform
# use_cov is a boolean flag to use covariates/X or not
# test_rank is a flag to test the rank of X (can be slow based on size of D, Y)
# returns A, B, beta, Yhat (array)
additive_regression <- function(Y, D, X, type="biten", use_cov=T, test_rank=F, penalty=NA, whichlambda="min")
{
  
  # Check sizes
  if(length(dim(Y)) != 3 ){ stop("Y is not a 3-mode array") }
  if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  
  
  # Find sizes, assuming dim(D) = dim(Y)
  S <- nrow(D[,,1])
  L <- ncol(D[,,1])
  tmax <- dim(D)[3]
  
  
  # Build X matrix
  Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov)  # build design matrix
  
  
  if(test_rank){   # see if X is full rank and warn if not
    if(rankMatrix(Xreg)[1] < ncol(Xreg)){ warning("X is not full rank") }
  }
  
  
  # Perform regression and pull out coefficients
  # fit <- lm(c(Y) ~ -1 + Xreg)
  if(is.numeric(penalty)){
    yfit <- c(Y)
    keep <- !is.na(c(yfit))   # glmnet cannot handle NAs
    cv.fit <- cv.glmnet(x=Xreg[keep,], y=yfit[keep], alpha=penalty, family='gaussian', intercept=F)
    # cv.fit <- cv.glmnet(Xreg, y=as.factor(c(Y)), alpha=penalty, family='binomial', intercept=F)
    fit <- cv.fit$glmnet.fit
    if(whichlambda == "1sd"){
      col <- cv.fit$lambda == cv.fit$lambda.1se  # .min??
    } else if (whichlambda == "min"){
      col <- cv.fit$lambda == cv.fit$lambda.min  # .min??
    } else {
      warning("Invalid choice of lambda. Defauling to 1se choice (maximum regularizaiton)")
      col <- cv.fit$lambda == cv.fit$lambda.1se  # .min??
    }
    coefs <- fit$beta[, col]   # also can use <- coef(cv.fit, s="lambda.1se") # or "lambda.min"
    
    Yhat <- array((predict(fit, newx=Xreg, type="response")[, col]), c(S,L,tmax))    # need own predict method
    
  } else {
    if(S < 30 | strtrim(type, 3) == "bit"){   # fit.lm fast for small dimensions
      remove <- which(is.na(c(Y)))  # NAs to remove
      if(length(remove) > 0){
        x <- Xreg[-remove, ]
        y <- c(Y)[-remove]
      } else {
        x = Xreg  ;  y = c(Y)
      }
      fit <- lm.fit(x=x, y=y)   # faster than lm()
      coefs <- fit$coef
      Yhat <- array(NA, c(S,L,tmax))
      Yhat[!is.na(Y)] <- fit$fitted.values
    } else {
      fit <- solve_lm(x=Xreg, y=c(Y))  # handles NAs internally
      coefs <- fit$coef
      Yhat <- array(fit$pred, c(S,L,tmax))
    }
  }
  A <- matrix(coefs[1:S^2], nrow=S)
  B <- matrix(coefs[S^2 + 1:L^2], nrow=L)
  
  if(use_cov){
    p <- dim(X)[4]
    beta <- matrix(coefs[S^2 + L^2 + 1:p], ncol=p)
  } else { beta <- NA }
  
  return(list(A=A, B=B, beta=beta, Yhat=Yhat))
}


# Sparse BLIN multipartite estimation
# Perform regression of Y onto D for a given type
# Y is an m-mode array of responses of the same size as D
# D is an m-length list of m-mode lagged arrays
# X is a m+1-mode array of covariates (not yet implemented)
# use_cov is a boolean flag to use covariates/X or not
# returns m-length list of coefficients
additive_regression_multi <- function(Y, D, X, use_cov=FALSE, test_rank=FALSE, penalty=1,   whichlambda="min")
{
  type="biten"
  penalty <- as.numeric(penalty)
  
  # dimensions
  ms <- dim(Y)
  if(as.logical(use_cov)){
    stop("Not yet coded for covariates")
  }
  
  # Build X matrix
  Xreg <- build_design_additive_multi2(D)
  if(test_rank){   # see if X is full rank and warn if not
    if(rankMatrix(Xreg)[1] < ncol(Xreg)){ warning("X is not full rank") }
  }
  
  
  # Perform regression and pull out coefficients
    yfit <- c(Y)
    keep <- !is.na(c(yfit))   # glmnet cannot handle NAs
    cv.fit <- cv.glmnet(x=Xreg[keep,], y=yfit[keep], alpha=penalty, family='gaussian', intercept=F)
    # cv.fit <- cv.glmnet(Xreg, y=as.factor(c(Y)), alpha=penalty, family='binomial', intercept=F)
    fit <- cv.fit$glmnet.fit
    if(whichlambda == "1sd"){
      col <- cv.fit$lambda == cv.fit$lambda.1se  # .min??
    } else if (whichlambda == "min"){
      col <- cv.fit$lambda == cv.fit$lambda.min  # .min??
    } else {
      warning("Invalid choice of lambda. Defauling to min choice (maximum regularizaiton)")
      col <- cv.fit$lambda == cv.fit$lambda.min  # .min??
    }
    coefs <- fit$beta[, col]   # also can use <- coef(cv.fit, s="lambda.1se") # or "lambda.min"
    
  B <- vector("list", length(ms) - 1)
  mtot <- 0   # running total length
  for(i in 1:length(B)){
    B[[i]] <- matrix(coefs[1:(ms[i]^2) + mtot], ms[i], ms[i] )
    mtot <- mtot + ms[i]^2
  }  
  
  if(use_cov){
    p <- dim(X)[length(dim(X))]
    beta <- matrix(coefs[mtot + 1:p], ncol=p)
  } else { beta <- NA }
  
  return(list(B=B, beta=beta))
}



# Iterative version of additive code. Seems to reproduce regression results (faster for SADD with large S,L), but 
#  there are discrepancies still with BiTEN FULL regression when using Covariates. 
fit_additive_iterative <- function(Y, D, X, type="sadd", use_cov=T,  verbose=T, printout=10, tol=1e-6)
{
  
  # Check sizes
  if(length(dim(Y)) != 3 ){ stop("Y is not a 3-mode array") }
  if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  
  
  # Find sizes, assuming dim(D) = dim(Y)
  S <- nrow(D[,,1])
  L <- ncol(D[,,1])
  tmax <- dim(D)[3]
  i_na <- is.na(Y)   # indices of NAs
  
  
  # Initialize A, B, and beta
  A <- diag(S)  ;  B <- diag(L)
  if(strtrim(type, 3) == "sad"){
    JD <- amprod(D, matrix(1, S, S), 1)
    DJ <- amprod(D, matrix(1, L, L), 2)
    Y[i_na] <- ( amprod(DJ, A, 1) + amprod(JD, t(B), 2) )[i_na]   # initialize NAs 
    DregA <- matrix(0, tmax*S*L, S^2)
    DregB <- matrix(0, tmax*S*L, L^2)
    for(t in 1:tmax){  # vec A then vec B
      DregA[ 1:(S*L) + (t - 1)*S*L, ] <- kronecker(t(DJ[,,t]), diag(S))   # A columns
      DregB[ 1:(S*L) + (t - 1)*S*L, ] <- kronecker(diag(L), JD[,,t])    # B columns
    }
  } else if (strtrim(type, 3) == "bit") {
    Y[i_na] <- ( amprod(D, A, 1) + amprod(D, t(B), 2) )[i_na]   # initialize NAs
    DregA <- matrix(0, tmax*S*L, S^2)
    DregB <- matrix(0, tmax*S*L, L^2)
    for(t in 1:tmax){  # vec A then vec B
      DregA[ 1:(S*L) + (t - 1)*S*L, ] <- kronecker(t(D[,,t]), diag(S))   # A columns
      DregB[ 1:(S*L) + (t - 1)*S*L, ] <- kronecker(diag(L), D[,,t])    # B columns
    } 
  } else (stop("Invalid fit type; use sadd or biten"))
  if(use_cov){
    betaOLS <- solve_lm(t(mat(X, 4)), c(Y))$coef  
  } else { betaOLS <- NA}
  
  
  change <- 100
  count <- 0
  # SSE <- 1
  while(change > tol){
    
    # Update for beta
    if(strtrim(type, 3) == "sad"){
      DB <- amprod(JD, t(B), 2)
      Yab <-  amprod(DJ, A, 1) + DB
    } else {
      DB <- amprod(D, t(B), 2)
      Yab <-  amprod(D, A, 1) + DB 
    }
    Ystar <- Y - Yab
    if(use_cov){
      beta <- matrix(solve_lm(t(mat(X, 4)), c(Ystar))$coef, nrow=1)
      Xbeta <- drop(amprod(X, beta, 4))
    } else {
      beta <- NA
      Xbeta <- 0
    }
    
    if(count == 0) {   # save initial LL and beta
      beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( (Ystar - Xbeta )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Update A and B
    Yxb <- Y - Xbeta - DB
    A <- matrix(solve_lm(DregA, c(Yxb))$coef, S, S)   # A update
    if(strtrim(type, 3) == "sad"){
      AD <- amprod(DJ, A, 1)
    } else {
      AD <- amprod(D, A, 1)
    }
    Yxa <- Y - Xbeta - AD
    B <- matrix(solve_lm(DregB, c(Yxa))$coef, L, L)   # B update
    
    if(strtrim(type, 3) == "sad"){
      DB <- amprod(JD, t(B), 2)
      Yab <-  amprod(DJ, A, 1) + DB
    } else {
      DB <- amprod(D, t(B), 2)
      Yab <-  amprod(D, A, 1) + DB 
    }
    Yhat <- Xbeta + Yab
    Y[i_na] <- Yhat[i_na]  # update NAs
    
    LLnew <- -length(Y)*log( sum( (Y - Yhat )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    # SSEnew <- LLnew <- mean( (Y-Yhat)^2,na.rm=TRUE )
    # change <- (SSEnew - SSE)/SSE
    # SSE <- SSEnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T) {
    cat("\n************************************************ \n")
    
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  
  return(list(A=A, B=B, Yhat=Yhat, beta= beta, betaOLS = betaOLS, LL2 = LL, LL2_init=LLinit))
}


# Fit BLIN model using iterative updates of A and B
fit_MLE_blin_iterative <- function(Y, D, X, verbose=T, printout=10, tol=1e-10, init="random", sigma_init=1, use_cov=F, maxit=1e9)
{
  
  if(use_cov){ stop("Not including covariates yet")}
  
  # Get sizes and check
  if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
  if(length(dim(Y)) < 3){ stop("Y is not an array")}
  
  i_na <- is.na(Y)  # save NAs
  Y[i_na] <- 0  # set to zero
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  tmax <- dim(Y)[3]
  if(use_cov){
    p <- dim(X)[4]
  }
  
  
  # matrices to regress upon
  # Xa <- matrix(0, tmax*S*L, S^2)
  # Xb <- matrix(0, tmax*S*L, L^2)
  # for(t in 1:tmax){
  #   Xa[1:(S*L) + (t-1)*S*L, ] <- kronecker(t(D[,,t]), diag(S))
  #   Xb[1:(S*L) + (t-1)*S*L, ] <- kronecker(diag(L), D[,,t])
  # }
  # XXa <- crossprod(Xa)
  # XXb <- crossprod(Xb)
  XXa <- solve( Reduce("+", lapply(1:tmax, function(t) tcrossprod(D[,,t],D[,,t])) ) )
  XXb <- solve( Reduce("+", lapply(1:tmax, function(t) crossprod(D[,,t],D[,,t])) ) )
  
  
  
  # Initialize
  if(strtrim(init,3) == "ran") {
    if(verbose == T){
      cat("Randomly initializing A,B, beta \n")
    }
    A <- matrix(rnorm(S^2, 0, sigma_init), S, S)
    BT <- matrix(rnorm(L^2, 0, sigma_init), L, L)
    if(use_cov){
      beta <- rnorm(p+1, 0, sigma_init)
    }
    
  } else if (strtrim(init, 3) == "reg") {
    if(verbose == T){
      cat("Initializing A,B,beta with regression \n\n")
    }
    
    avec <- solve(XXa, crossprod(Xa, Y))
    bvec <- solve(XXb, crossprod(Xb, Y))
    
    A <- matrix(avec, S, S)
    BT <- t(matrix(bvec, L, L))
    
    
  } else if (strtrim(init, 1) == "I") {
    if(verbose == T){
      cat("Initializing A, B as identity \n\n")
    }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
    
    A <- diag(S)
    BT <- diag(L)
    
  } else { stop("Invalid initialization type")}
  
  
  Ainit <- A
  BTinit <- BT
  
  
  # if(use_cov){
  #   betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
  # } else { betaOLS <- NA}
  
  # Find optimal values
  change <- 100
  count <- 0
  
  while(change > tol & count < maxit){
    
    if(count == 0) {   # save initial LL and beta
      # beta_init <- beta
      LLinit <- LL <- -length(Y)*log( sum( ( Y - amprod(D, A, 1) - amprod(D, BT, 2) )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    }
    
    # Ytilde <- Y - amprod(D, A, 1)
    # bvec <- solve(XXb, crossprod(Xb, Ytilde))
    # BTnew <- t(matrix(bvec, L, L))
    # 
    # Ytilde <- Y - amprod(D, BTnew, 2)
    # avec <- solve(XXa, crossprod(Xa, Ytilde))
    # Anew <- matrix(avec, S, S)
    
    Ytilde <- Y - amprod(D, A, 1)
    XYb <- Reduce("+", lapply(1:tmax, function(t) crossprod(D[,,t], Ytilde[,,t])))
    BTnew <- t( XXb %*% XYb )
    
    Ytilde <- Y - amprod(D, BTnew, 2)
    XYa <- Reduce("+", lapply(1:tmax, function(t) tcrossprod(Ytilde[,,t], D[,,t])))
    Anew <- XYa %*% XXa
    
    
    changeAB <- max(abs(c(c(A - Anew), c(BT-BTnew))))  # doesn't make much difference stopping A,B vs using LL
    A <- Anew
    BT <- BTnew
    
    Yhat <- amprod(D, A, 1) + amprod(D, BT, 2)  # new estimate
    Y[i_na] <- Yhat[i_na]  # update NAs
    
    # LLnew <- -sum( (Ytilde - (amprod(D, A, 1) + amprod(D, BT, 2)) )^2)
    LLnew <- -length(Y)*log( sum( (Y - Yhat )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
    change <- abs(LLnew - LL) 
    LL <- LLnew
    
    count <- count + 1
    if(count%%printout == 0 & verbose == T){
      cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
    }
  }
  
  if(verbose == T){
    cat("\n************************************************ \n")
    
    # cat("True 2LL*tau^2: \t", LLtrue, "\n")
    cat("Initial 2LL: \t", LLinit, "\n")
    cat("Final 2LL: \t", LL, "\n")
    
    cat("\n************************************************ \n \n")
    
    # cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
    # cat("Est. beta coefficients are: \t\t", beta, "\n")
  }
  
  return(list(A=A, B=t(BT), Yhat = Yhat, LLt2 = LL, LLt2_init=LLinit))   # beta= beta, betaOLS = betaOLS, 
}


# Fit binary datasets using MLE model
#  Handles sparse representation of BiTEN model, or standard representation of BiTEN and SADD models
#
fit_MLE_binary <- function(Y, D, X, use_cov=T, type="biten", penalty=1, sparsedata=F, writeXreg=T, readXreg=T, nameXreg="Xreg.txt", seed=NA, loss="auc")
{
  if(is.numeric(seed)){ set.seed(seed)}
  
  if(!sparsedata){
    # Check sizes
    if(length(dim(Y)) != 3 ){ stop("Y is not a 3-mode array") }
    if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
    if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
    
    # Find sizes, assuming dim(D) = dim(Y)
    S <- nrow(D[,,1])
    L <- ncol(D[,,1])
    tmax <- dim(D)[3]
    
    # Build X matrix
    Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata)  # build design matrix
    
    if(use_cov){
      p <- dim(X)[4]
    }
    
  } else if (sparsedata){
    
    # Check if i,j,t in column names
    if( !("i" %in% names(Y)) | !("j" %in% names(Y)) | !("t" %in% names(Y))){stop("Y must have column names i,j, and t")}
    
    # Find sizes
    S <- max(D$i)
    L <- max(D$j)
    tmax <- max(D$t)
    
    # Sort Y
    rows <- Y$i + (Y$j-1)*S + (Y$t-1)*S*L   # unfolded indices
    indices <- order(rows)
    Y <- Y$ratification_year[indices]  # 1's and zeros, in order
    
    # Build X matrix (full)
    if(readXreg & nameXreg %in% list.files()){
      cat("\n Reading in Xreg \n")
      load(nameXreg)  # should be already sorted/subsetted
    } else {
      if(readXreg){cat("\n nameXreg file not found. Generating Xreg \n")}
      Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata)  # build design matrix
      cat("Design built, dimensions ", dim(Xreg), "\n")
      rows <- rows[indices]  # sorted rows to keep
      Xreg <- Xreg[rows, ]  # subset
      Xreg <- Xreg[,-(S^2 + L^2 + 1)]  # remove intercept
    }

    if(use_cov){
      p <- ncol(X) - 3 - 1   # -3 for i,j,t, -1 for intercept removed
    }
    
  } else { stop("sparsedata must be true/false")}
  
  
  Xreg[which(is.na(Xreg), arr.ind=T)] <- 0   # set NAs to zero
  if(writeXreg){  # write out
    save(Xreg, file=nameXreg)
  }
  
  # Bookkeep zero columns
  remove <- which(colSums(Xreg) == 0)   
  keep <- (1:ncol(Xreg))[-remove]
  
  # Perform regression and pull out coefficients
  if(is.numeric(penalty)){  # blend of lasso and ridge based on penalty value
    cv.fit <- cv.glmnet(Xreg[,keep], y=as.factor(c(Y)), alpha=penalty, family='binomial', intercept=T, type.measure=loss, standardize=F, lambda=10^seq(-6,6,length.out = 100))
    # cv.fit <- cv.glmnet(Xreg, y=as.factor(c(Y)), alpha=penalty, family='binomial', intercept=F)
    fit <- cv.fit$glmnet.fit
    col <- cv.fit$lambda == cv.fit$lambda.min  # .1se
    coefs <- coef(cv.fit)   # usees lambda.min by default
    Yhat <- predict(cv.fit, newx=Xreg[,keep], type="response")   # uses lambda.min by default
    
    fitout <- cv.fit
    
    
  } else {
    # fit <- glm(c(Y) ~ Xreg[,keep] - 1, family="binomial")   # non-lasso'd
    # fit <- glm.fit(x=Xreg, y=c(Y), family=binomial())   # non-lasso'd
    # coefs <- fit$coef
    # Yhat <- fit$fitted.values
    # Use biglm??
    # nlambda=1, 
    fit <- glmnet(x=Xreg[,keep], y=as.factor(c(Y)), family="binomial", lambda=10^seq(1,-10, length.out = 50), intercept=T, standardize = F, maxit=1e6, lambda.min.ratio=1e-12)
    coefs <- coef(fit, s=0)  # include intercept
    Yhat <- predict(fit, newx=Xreg[,keep], type = "response", s=0)
    fitout <- fit
  }
  
  oldcoefs <- coefs   
  coefs <- rep(NA, ncol(Xreg))
  coefs[keep] <- oldcoefs[-1]   # add in NAs for un-estimable quantities, -1 for intercept
  coefs <- c(oldcoefs[1], coefs)  # intercept
  
  A <- matrix(coefs[1:S^2 + 1], S, S)   # + 1 for intercept
  B <- matrix(coefs[S^2 + 1 + 1:L^2], L, L)
  if(use_cov){
    beta <- coefs[c(1,1+S^2 + L^2 + 1:p)]
    beta <- matrix(beta, nrow=length(beta))
  } else { beta <- NA }
  
  if(!sparsedata){
    Yhat <- array(round(Yhat), dim(D))   # reshape 
  }
  
  return(list(A=A, B=B, beta=beta, Yhat=Yhat, fit=fitout, Yreg=Y, Xreg=Xreg, Xkeep=keep, Xremove=remove))
}


#### OLD FUNCTION ####
# P/R cross-validation with automatic lambda selection, standardized X
cv_pr <- function(Y, D, X, outdir, use_cov=T, type="biten", penalty=1, sparsedata=F, writeXreg=T, readXreg=T, 
                   nameXreg="Xreg.txt", seed=NA, ncv=10, verbose=F, maxit=1e5, ncores=1, thresh=1e-7,
                  cutA=0, cutB=0, global_cov=F, remove_dups=F, keep_one_dup=F)
{
  if(! is.numeric(penalty)){ stop("Penalized regression method requires numeric penalty value (alpha") }
  
  dir.create(outdir, showWarnings = F)
  
  if(!sparsedata){
    
    stop("Not implemented for non-sparse data")
    # # Check sizes
    # if(length(dim(Y)) != 3 ){ stop("Y is not a 3-mode array") }
    # if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
    # if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
    # 
    # # Find sizes, assuming dim(D) = dim(Y)
    # S <- nrow(D[,,1])
    # L <- ncol(D[,,1])
    # tmax <- dim(D)[3]
    # 
    # # Build X matrix
    # Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata)  # build design matrix
    # 
    # if(use_cov){
    #   p <- dim(X)[4]
    # }
    
  } else if (sparsedata){
    
    # Check if i,j,t in column names
    if( !("i" %in% names(Y)) | !("j" %in% names(Y)) | !("t" %in% names(Y))){stop("Y must have column names i,j, and t")}
    
    # Add global covariate if desired and not already present
    if(global_cov & !("global_ind" %in% names(X))){
      X <- add_global(X, sort(unique(Y$treaty)), readdir=getwd())
    }
    
    # Find sizes
    S <- max(D$i)
    L <- max(D$j)
    tmax <- max(D$t)
    
    # Sort Y
    rows <- Y$i + (Y$j-1)*S + (Y$t-1)*S*L   # unfolded indices
    indices <- order(rows)
    Yreg <- Y$ratification_year[indices]  # 1's and zeros, in order
    
    # Build X matrix (full)
    if(readXreg & nameXreg %in% list.files()){
      cat("\n Reading in Xreg \n")
      load(nameXreg)  # should be already sorted/subsetted
    } else {
      if(readXreg){cat("\n nameXreg file not found. Generating Xreg \n")}
      Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata, S=S, L=L, tmax=tmax)  # build design matrix
      cat("Design built, dimensions ", dim(Xreg), "\n")
      rows <- rows[indices]  # sorted rows to keep
      Xreg <- Xreg[rows, ]  # subset
      Xreg <- Xreg[,-(S^2 + L^2 + 1)]  # remove intercept
    }
    
    if(use_cov){
      p <- ncol(X) - 3 - 1   # -3 for i,j,t, -1 for intercept removed
    }
    
  } else { stop("sparsedata must be true/false")}
  
  
  Xreg[which(is.na(Xreg), arr.ind=T)] <- 0   # set NAs to zero
  if(writeXreg){  # write out
    save(Xreg, file=nameXreg)
  }
  
  # Bookkeep zero columns
  remove <- which(colSums(Xreg[,1:(S^2)]) <= cutA)    # remove columns with \le cutA 1's in them
  remove <- c(remove, S^2 + which(colSums(Xreg[,S^2 + 1:(L^2)]) <= cutB)) 
  keep <- (1:ncol(Xreg))[-remove]
  
  # Remove duplicate columns
  if(remove_dups & !keep_one_dup){
    nonid <- checkcols(Xreg[,keep])  # columns that are duplicates
    nonid <- sort(unique(c(nonid)))   # vectorize
    nonid <- keep[nonid]    # in terms of original indices
    remove <- c(remove, nonid)   # modify remove and keep columns
    keep <- (1:ncol(Xreg))[-remove]
    
    nonid <- list(keep=NA, remove=nonid)
    
  } else if (!remove_dups) { 
    nonid <- NA
  } else if (remove_dups & keep_one_dup) {
    nonid <- checkcols(Xreg[,keep])  # columns that are duplicates
    # nonid <- sort(unique(c(nonid)))   # vectorize
    nonid <- cbind(keep[nonid[,1]], keep[nonid[,2]])    # in terms of original indices
    
    dups <- sort(unique(c(nonid)))   # columns that are duplicated at some point, in terms of original indices
    a <- as.vector(t(Xreg[,dups]) %*% rnorm(nrow(Xreg)))   # wp1 only same if same number
    ldups <- lapply(sort(unique(a)), function(z) dups[which(a == z)])  # list of all duplicate subsets
    
    remove1 <- unlist(sapply(ldups, function(z) z[-1]))
    remove <- c(remove, remove1)   # modify remove and keep columns
    keep <- (1:ncol(Xreg))[-remove]
    
    nonid_remove <- remove1
    nonid_keep <- sapply(ldups, function(z) z[1])
    nonid <- list(keep=nonid_keep, remove=nonid_remove)
    
    nonid_remove <- sapply(ldups, function(z) z[-1])
    save(ldups, nonid_keep, nonid_remove, keep, remove, file=file.path(outdir, "nonid_data.RData"))
  } else { stop("Invalid input for dupliate column logicals")}
  
  
  
  #### Perform full fit and extract coefficients for lambda
  if(is.numeric(seed)){ set.seed(seed)}  # set seed
  fitfull <- fit <- glmnet(Xreg[,keep], y=as.factor(c(Yreg)), alpha=penalty, family='binomial', intercept=T, lambda.min.ratio=1e-12, standardize=T, maxit=maxit, nlambda=100, thresh=thresh)
  save(fit, file=file.path(outdir, paste0("fullfit.RData")))
  workedfull <- length(fitfull$lambda) 
  if(verbose){
    cat("full fit done, nlambda", workedfull, "\n")
  }
  lambdas <- fitfull$lambda   # save selected lambda sequence
  ####
  
  
  
  
  
  #### Cross-validate
  lambdas <- sort(lambdas, decreasing = T)  # decreasing sequence
  ones_check <- rep(0, ncv)
  if(is.numeric(seed)){ set.seed(seed)}  # set seed
  count <- 0
  while(any(ones_check == 0) & count < 1e4){   # if there are any partitions without 1's, and make loop finite
    count <- count + 1
    cat("cv partition count", count, "\n")
    rowscramble <- sample(1:nrow(Xreg), size = nrow(Xreg), replace = F)
    cvs <- rowscramble %% ncv + 1
    ones_check <- sapply(1:ncv, function(z) sum(Yreg[cvs == z]))   # number of 1's in Yreg for each cv partition
  }
  cvms <- matrix(NA, length(lambdas), ncv)
  assign("last.warning", NULL, envir = baseenv())   # clear warnings
  
    ####
    save(S,L,tmax,cutA,cutB,thresh,nameXreg,keep,remove,seed,ncv,maxit,cvs, file=file.path(outdir, paste0("prefit.RData")))
    ####
  
  if(ncores==1){
    fitlist <- vector("list", ncv)
    for(i in 1:ncv){
      test <- which(cvs == i)
      train <- (1:nrow(Xreg))[-test]
      
      fitlist[[i]] <- fit <- glmnet(Xreg[train,keep], y=as.factor(c(Yreg[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-12, standardize=T, thresh = thresh)
      worked <- length(fit$lambda)
      
      save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
      
      if(verbose){
        cat("\n************************************************\n")
        cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
        print(warnings())
        cat("\n************************************************")
        cat("\n")
        assign("last.warning", NULL, envir = baseenv())
      }
    }
    
  } else if (ncores > 1 ){
    writeLines(c(""), file.path(outdir, "log.txt"))
    
    registerDoMC(cores=ncores)
    mcoptions <- list(preschedule=FALSE, set.seed=T)
    fitlist <- foreach(i=1:ncv, .options.multicore=mcoptions, .packages=c("glmnet") ) %dopar% {  # .combine =cbind
      test <- which(cvs == i)
      train <- (1:nrow(Xreg))[-test]
      
      fit <- glmnet(Xreg[train,keep], y=as.factor(c(Yreg[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-12, standardize=T, thresh = thresh)
      save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
      worked <- length(fit$lambda)
      sink(file.path(outdir, "log.txt"), append=TRUE)   # write out to log file
      cat("\n************************************************\n")
      cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
      print(warnings())
      cat("\n************************************************")
      cat("\n")
      fit
    }
    
  } else {stop("ncores must be numeric >= 1")}
  
  for(i in 1:ncv){
    test <- which(cvs == i)
    train <- (1:nrow(Xreg))[-test]
    
    fit <- fitlist[[i]]
    #glmnet(Xreg[train,keep], y=as.factor(c(Y[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit)
    worked <- length(fit$lambda)
    
    if(worked > 1){   # if fit worked
      Yhats <- predict(fit, newx=Xreg[test,keep], type="response")
      
      pr_temp <- rep(NA, ncol(Yhats))
      for(j in 1:ncol(Yhats)){
        # pr_temp[j] <- pr_curve(Yhats[,j], Y[test], n=2000)$auc
        # pr_temp[j] <- roc_curve(Yhats[,j], Y[test], n=min(length(Y[test]), 2500))$auc
        # pr_temp[j] <- simple_roc_auc(Yhats[,j], Y[test])
        pr_temp[j] <- simple_pr_auc(Yhats[, j], Yreg[test])
        # pr_temp[j] <- auc(roc(Yhats[,j], Y[test]))
      }
      
      cvms[match(fit$lambda, lambdas), i] <- pr_temp   # save lambda values
      # save(fit, pr_temp, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
    }
    if(verbose){
      cat("processed cv", i, "of", ncv, ", nlambda", worked, "\n")
    }
  }
  
  # Calculate minimum lambda value
  mean_cvms <- apply(cvms, 1, mean, na.rm=T)
  imin <- which(mean_cvms == max(mean_cvms, na.rm=T))
  lambda_min <- lambdas[imin]
  ####
  
  
  #### Calculate results from full fit
  fit <- fitfull
  if(workedfull > 1){
    Yhats <- predict(fitfull, newx=Xreg[,keep], type="response")
    cvm_full <- rep(NA, ncol(Yhats))
    for(j in 1:ncol(Yhats)){
      cvm_full[j] <- simple_pr_auc(Yhats[,j], Yreg) # pr_curve(Yhats[,j], Y, n=min(length(Y), 2500))$auc
    }
    Yhat <- predict(fit, newx=Xreg[,keep], type="response", s=lambda_min)
    coefs <- coef(fit, s=lambda_min)
    
    # inflate coefs with NAs for non-estimable quantities
    oldcoefs <- coefs   
    coefs <- rep(NA, ncol(Xreg))
    coefs[keep] <- oldcoefs[-1]   # add in NAs for un-estimable quantities, -1 for intercept
    coefs <- c(oldcoefs[1], coefs)  # intercept
    
    A <- matrix(coefs[1:S^2 + 1], S, S)   # + 1 for intercept
    B <- matrix(coefs[S^2 + 1 + 1:L^2], L, L)
    if(use_cov){
      beta <- coefs[c(1,1+S^2 + L^2 + 1:p)]
      beta <- matrix(beta, nrow=length(beta))
    } else { beta <- NA }
  } else {
    fit <- A <- B <- Yhat <- NA
  }
  ####
  
  
  return(list(A=A, B=B, beta=beta, Yhat=Yhat, fit=fitfull, Yreg=Yreg, Xreg=Xreg, Xkeep=keep, Xremove=remove, cvms=cvms, cvm_full=cvm_full, lambda_min=lambda_min, imin=imin, nonid=nonid))
}





# cutA=cutB=starting column sum to remove, with while-loop
# internal processing of duplicate columns and removed, etc.
cv_pr_envtreat <- function(Y, D, X, outdir, use_cov=T, type="biten", penalty=1, sparsedata=F, writeXreg=T, readXreg=T, 
                           nameXreg="Xreg.txt", seed=NA, ncv=10, verbose=F, maxit=1e5, ncores=1, thresh=1e-7,
                           cutA=0, cutB=0, cAB_upper=25, global_cov=F, remove_dups=F, keep_one_dup=F, cowsfile="COW_country_codes.csv", withnet=T, penbeta=T, cvtype="random")
{
  
  if(! is.numeric(penalty)){ stop("Penalized regression method requires numeric penalty value (alpha") }
  if(withnet){ use_cov <- T}
  dir.create(outdir, showWarnings = F)
  cows <- read.csv(cowsfile, header=T, stringsAsFactors = F)
  
  
  if(!sparsedata){
    
    stop("Not implemented for non-sparse data")
    
  } else if (sparsedata){
    
    # Check if i,j,t in column names
    if( !("i" %in% names(Y)) | !("j" %in% names(Y)) | !("t" %in% names(Y))){stop("Y must have column names i,j, and t")}
    
    # Add global covariate if desired and not already present
    if(global_cov & !("global_ind" %in% names(X))){
      X <- add_global(X, sort(unique(Y$treaty)), readdir=getwd())
    }
    
    # Find sizes
    S <- max(c(max(D$i), max(Y$i), max(X$i)))
    L <- max(c(max(D$j), max(Y$j), max(X$j)))
    tmax <- max(c(max(D$t), max(Y$t), max(X$t)))
    countries <- sort(unique(Y$cowcode))
    treaties <- sort(unique(Y$treaty))
    
    # Sort Y
    rows <- Y$i + (Y$j-1)*S + (Y$t-1)*S*L   # unfolded indices
    indices <- order(rows)
    Yreg <- Y$ratification_year[indices]  # 1's and zeros, in order
    
  } else { stop("sparsedata must be true/false")}
  
  if(withnet){
    cutA <- cutB <- min(cutA,cutB)
    checkfit <- F
    
    while(!checkfit & (cutA <= cAB_upper)){
      # Build X matrix (full)
      if(readXreg & nameXreg %in% list.files()){
        cat("\n Reading in Xreg \n")
        load(nameXreg)  # should be already sorted/subsetted
      } else {
        if(readXreg){cat("\n nameXreg file not found. Generating Xreg \n")}
        Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata, S=S, L=L, tmax=tmax)  # build design matrix
        cat("Design built, dimensions ", dim(Xreg), "\n")
        rows <- rows[indices]  # sorted rows to keep
        Xreg <- Xreg[rows, ]  # subset
        Xreg <- Xreg[,-(S^2 + L^2 + 1)]  # remove intercept
      }
      
      if(use_cov){
        p <- ncol(X) - 3 - 1   # -3 for i,j,t, -1 for intercept removed
      }
      
      Xreg[which(is.na(Xreg), arr.ind=T)] <- 0   # set NAs to zero
      if(writeXreg){  # write out
        save(Xreg, file=nameXreg)
      }
      
      # Bookkeep zero columns
      remove <- which(colSums(Xreg[,1:(S^2)]) <= cutA)    # remove columns with \le cutA 1's in them
      remove <- c(remove, S^2 + which(colSums(Xreg[,S^2 + 1:(L^2)]) <= cutB)) 
      remove <- union(remove, S^2 + L^2 + which(colSums(Xreg[,(S^2 + L^2 +1):ncol(Xreg)]) == 0))
      keep <- (1:ncol(Xreg))[-remove]
      
      # Remove duplicate columns
      if(remove_dups & !keep_one_dup){
        stop("Not sure this works for keep_one_dup = F")
        # nonid <- checkcols(Xreg[,keep])  # columns that are duplicates
        # nonid <- sort(unique(c(nonid)))   # vectorize
        # nonid <- nonid_remove <- keep[nonid]    # in terms of original indices
        # remove <- c(remove, nonid)   # modify remove and keep columns
        # keep <- (1:ncol(Xreg))[-remove]
        # 
        # # nonid <- list(keep=NA, remove=nonid)
        # ldups <- nonid_keep <-  nonid_remove <- NA
        
      } else if (!remove_dups) { 
        ldups <- nonid_keep <-  nonid_remove <- NA
        nonid <- NULL
        
      } else if (remove_dups & keep_one_dup) {
        nonid <- checkcols(Xreg[,keep])  # columns that are duplicates
        
        # nonid <- sort(unique(c(nonid)))   # vectorize
        if(!is.null(nonid)){
          nonid <- cbind(keep[nonid[,1]], keep[nonid[,2]])    # in terms of original indices
          
          dups <- sort(unique(c(nonid)))   # columns that are duplicated at some point, in terms of original indices
          a <- as.vector(t(Xreg[,dups]) %*% rnorm(nrow(Xreg)))   # wp1 only same if same number
          ldups <- lapply(sort(unique(a)), function(z) dups[which(a == z)])  # list of all duplicate subsets
          
          remove1 <- unlist(sapply(ldups, function(z) z[-1]))
          remove <- c(remove, remove1)   # modify remove and keep columns
          keep <- (1:ncol(Xreg))[-remove]
          
          nonid_remove <- remove1
          nonid_keep <- sapply(ldups, function(z) z[1])
          nonid <- list(keep=nonid_keep, remove=nonid_remove)
          
          nonid_remove <- sapply(ldups, function(z) z[-1])
        } else {
          ldups <- nonid_keep <-  nonid_remove <- NA
        }
        
        # save(ldups, nonid_keep, nonid_remove, file=file.path(outdir, "nonid_data.RData"))
      } else { stop("Invalid input for dupliate column logicals")}
      
      if(remove_dups & !is.null(nonid)){
        nout <- cbind(unlist(ldups), rep(1:length(ldups), times=sapply(ldups, length)), 0)
        nout <- cbind(nout, nout[,1] <= S^2, nout[,1] > S^2)
        
        ia <- matrix(1:S^2, S,S)
        ib <- matrix(1:L^2 + S^2, L,L)
        cts <- t(sapply(nout[,1], function(z) rbind(which(ia == z, arr.ind=T), which(ib == z, arr.ind=T))))
        nout <- cbind(nout, cts, NA, NA)  
        
        nout[nout[,5]==1, 6:7] <- c(treaties[nout[nout[,5]==1, 6]], treaties[nout[nout[,5]==1, 7]])
        nout[nout[,4]==1, 8:9] <- c(countries[nout[nout[,4]==1, 6]], countries[nout[nout[,4]==1, 7]])
        nout[nout[,5]==0, 7:6] <- NA
        if(sum(nout[,4]) >= 1){  # A indicator
          nout <- cbind(nout, cows$StateAbb[match(nout[,8], cows$CCode)], cows$StateAbb[match(nout[,9], cows$CCode)])
        } else {nout <- cbind(nout, NA, NA)}
        colnames(nout) <- c("column", "nonid_group", "keep", "Aind", "Bind", "treaty1", "treaty2", "cow1", "cow2", "country1", "country2")
        nout <- nout[order(as.numeric(nout[,1])),]
        write.table(nout, file.path(outdir, "nonid_columns.txt"), row.names=F, quote=F)
      } else {
        nout <- NA
      }
      
      if(length(remove) > 0){
        ia <- matrix(1:S^2, S,S)
        ib <- matrix(1:L^2 + S^2, L,L)
        
        ra <- t(sapply(remove[remove <= S^2], function(z) which(ia == z, arr.ind=T)))
        ra <- cbind(ra[,2], ra[,1])   # transpose!
        rb <- t(sapply(remove[remove >  S^2 & remove <=  S^2 + L^2], function(z) which(ib == z, arr.ind=T)))
        
        rb <- cbind(rb, treaties[rb[,1]], treaties[rb[,2]], NA, NA, NA, NA)
        ra <- cbind(ra, NA, NA, countries[ra[,1]], countries[ra[,2]])
        ra <- cbind(ra, cows$StateNme[match(ra[,5], cows$CCode)], cows$StateNme[match(ra[,6], cows$CCode)])
        
        
        allrem <- rbind(ra,rb)
        allrem <- cbind(allrem, 1*(remove[remove <= S^2 + L^2] %in% nonid_remove))
        colnames(allrem) <- c("i", "j", "treaty1", "treaty2", "cow1", "cow2", "state1", "state2", "nonid")
        
        rbeta <- remove[remove >  S^2 + L^2] - S^2 - L^2
        if(!is.na(rbeta[1])){
          allrem <- cbind(allrem, NA, NA)
          colnames(allrem) <- c("i", "j", "treaty1", "treaty2", "cow1", "cow2", "state1", "state2", "nonid", "ibeta", "coef")
          rc <- cbind(NA, NA, NA, NA, NA, NA, NA, NA, rbeta, (names(X)[!(names(X) %in% c("i","j","t","intercept"))])[rbeta])
        }
        
        write.table(nout, file.path(outdir, "allremoved.txt"), row.names=F, quote=F)
      }
      
      # Full Removal
      remove <- union(remove, unlist(nonid_remove))
      keep <- (1:ncol(Xreg))[-remove]      
      
      # Penalty factor
      pf <- rep(1, length(keep))
      if(!penbeta & use_cov){   # add some zeros to not penalize beta coefficients if desired
        pf <- rep(1, length(keep))
        beta_cols <- intersect((S^2 + L^2 + 1):(ncol(Xreg)), keep)
        pf[keep %in% beta_cols] <- 0   # NO penalty on betas when fed to lasso procedure
      }
      
      #### Save inputs 
      save(ldups, nonid_keep, nonid_remove, keep, remove, Xreg, Yreg, 
           Y, D, X, outdir, use_cov, type, penalty, sparsedata, writeXreg, readXreg, 
           nameXreg, seed, ncv, verbose, maxit, ncores, thresh,
           cutA, cutB, global_cov, remove_dups, keep_one_dup, countries, treaties, 
           nout, allrem, withnet, cowsfile, pf,
           file=file.path(outdir, "prefit2.RData"))
      ####
      
      cat("************* Fit for cutA=cutB=",cutA,"*******************\n\n")
      if(strtrim(cvtype,4) == "year"){ yearsin <- Y$year[indices] } else {yearsin <- NA}
      results <- cv_pr_lasso(Yreg, Xreg[,keep], outdir, penalty, seed=seed, ncv=ncv, verbose=verbose, maxit=maxit, ncores=ncores, thresh=thresh, pf=pf, cvtype=cvtype, years=yearsin)  
      
      cat("Done with while-loop, chose cutA=cutB=", cutA,"\n")
      
      checkfit <- results$nl_full > 25 & all(results$nl_cv > .9*results$nl_full)
      if(!checkfit){
        cutA <- cutB <- cutA + 1
      }
    } # end while loop
    
    
  } else {   # else for withnet
    Xreg <- X[,!(names(X) %in% c("i","j","t","intercept"))]
    Xreg[is.na(Xreg)] <- 0
    rows <- X$i + (X$j-1)*S + (X$t-1)*S*L   # unfolded indices
    indices <- order(rows)
    Xreg <- as.matrix(Xreg[indices, ])  # 1's and zeros, in order
    
    keep <- 1:ncol(Xreg)
    p <- ncol(Xreg)
    remove <- NA
    ldups <- nonid_keep<- nonid_remove <- nout <- allrem <- NA
    
    #### Save inputs 
    save(ldups, nonid_keep, nonid_remove, keep, remove, Xreg, Yreg, 
         Y, D, X, outdir, use_cov, type, penalty, sparsedata, writeXreg, readXreg, 
         nameXreg, seed, ncv, verbose, maxit, ncores, thresh,
         cutA, cutB, global_cov, remove_dups, keep_one_dup, countries, treaties, 
         nout, allrem, withnet, cowsfile,
         file=file.path(outdir, "prefit2.RData"))
    ####
    
    if(penbeta){
      if(strtrim(cvtype,4) == "year"){ yearsin <- Y$year[indices] } else {yearsin <- NA}
      
      results <- cv_pr_lasso(Yreg, Xreg[,keep], outdir, penalty, seed=seed, ncv=ncv, verbose=verbose, maxit=maxit, ncores=ncores, thresh=thresh, cvtype=cvtype, years=yearsin)  # implict even penalization (pf default input)
    } else {
      r <- glm.fit(x=cbind(1,Xreg[,keep]), y=Yreg, family =binomial())
      results$coefs <- r$coefficients
      results$Yhat <- r$fitted.values
      results$fit <- r
      results$cvm_full <- results$lambda_min <- results$imin <- results$nl_full <- results$nl_cv <- NA
      results$prfull <- simple_pr_auc(results$Yhat, Yreg)
    }
  }

  
  
  #### Form coefficient results
  allcoefs <- rep(NA, ncol(Xreg))
  allcoefs[keep] <- results$coefs[-1]
  allcoefs <- c(results$coefs[1], allcoefs)  # add in intercept
  
  
  if(withnet){
    A <- t(matrix(allcoefs[1:(S^2) + 1], S,S))   # transposed A for interpretability
    colnames(A) <- rownames(A) <- countries
    B <- matrix(allcoefs[1:(L^2) + 1 + S^2], L,L)   # A as in coding, not t(A)
    colnames(B) <- rownames(B) <- treaties
  } else { A <- B <- NA}
  
  if(use_cov){
    beta <- matrix(c(allcoefs[1], tail(allcoefs, p)), ncol=1)
    rownames(beta) <- c("intercept", names(X)[!(names(X) %in% c("i","j","t","intercept"))]) 
  } else {beta <- NA}
  
  Yhat <- results$Yhat
  fit <- results$fit
  cvms <- results$cvms
  cvm_full <- results$cvm_full
  lambda_min <- results$lambda_min
  imin <- results$imin
  nl_full <- results$nl_full
  nl_cv <- results$nl_cv
  prfull <- cvm_full[imin]
  if(!withnet & !penbeta){ prfull <- results$prfull}
  save(A,B,beta, results, Yhat, fit, cvms, cvm_full, lambda_min, imin, nl_full, nl_cv, prfull, file=file.path(outdir, "unprocessed_results.RData"))
  ####
  
  
  return(list(A=A, B=B, beta=beta, Yhat=results$Yhat, fit=results$fit, cvms=results$cvms, cvm_full=results$cvm_full, lambda_min=results$lambda_min, imin=results$imin, nl_full=results$nl_full, nl_cv=results$nl_cv, prfull=prfull, cutA=cutA))
}




# cutA=cutB=starting column sum to remove, with while-loop
# internal processing of duplicate columns and removed, etc.
cv_pr_generic <- function(Y, D, X, outdir, use_cov=T, type="biten", penalty=1, sparsedata=F, writeXreg=T, readXreg=T, 
                           nameXreg="Xreg.txt", seed=NA, ncv=10, verbose=F, maxit=1e5, ncores=1, thresh=1e-7,
                           cutA=0, cutB=0, global_cov=F, remove_dups=T, keep_one_dup=T, withnet=T, penbeta=T, cvtype="random",
                          response="response", anames=NULL, bnames=NULL)
{
  
  if(! is.numeric(penalty)){ stop("Penalized regression method requires numeric penalty value (alpha") }
  if(withnet){ use_cov <- T}
  dir.create(outdir, showWarnings = F)
  # cows <- read.csv(cowsfile, header=T, stringsAsFactors = F)
  
  
  
  if(!sparsedata){
    
    stop("Not implemented for non-sparse data")
    
  } else if (sparsedata){
    
    # Check if i,j,t in column names
    if( !("i" %in% names(Y)) | !("j" %in% names(Y)) | !("t" %in% names(Y))){stop("Y must have column names i,j, and t")}
    
    
    # Find sizes
    S <- max(c(max(D$i), max(Y$i), max(X$i)))
    L <- max(c(max(D$j), max(Y$j), max(X$j)))
    tmax <- max(c(max(D$t), max(Y$t), max(X$t)))
    cAB_upper=min(c(25, S, L))
    
    
    # Sort Y
    rows <- Y$i + (Y$j-1)*S + (Y$t-1)*S*L   # unfolded indices
    indices <- order(rows)
    Yreg <- Y[indices, response]  # 1's and zeros, in order
    
    
    if(is.null(anames)){
      countries <- anames <- paste0("a", 1:S)
    } else {
      countries <- anames
    }
    
    if(is.null(bnames)){
      treaties <- bnames <- paste0("b", 1:L)
    } else {
      treaties <- bnames
    }
    
  } else { stop("sparsedata must be true/false")}
  
  if(withnet){
    cutA <- cutB <- min(cutA,cutB)
    checkfit <- F
    
    while(!checkfit & (cutA <= cAB_upper)){
      # Build X matrix (full)
      if(readXreg & nameXreg %in% list.files()){
        cat("\n Reading in Xreg \n")
        load(nameXreg)  # should be already sorted/subsetted
      } else {
        if(readXreg){cat("\n nameXreg file not found. Generating Xreg \n")}
        Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata, S=S, L=L, tmax=tmax, response=response)  # build design matrix
        cat("Design built, dimensions ", dim(Xreg), "\n")
        rows <- rows[indices]  # sorted rows to keep
        Xreg <- Xreg[rows, ]  # subset
        Xreg <- Xreg[,-(S^2 + L^2 + 1)]  # remove intercept
        if(verbose){
          cat("At start: ncol(Xreg) = ", ncol(Xreg), "\n")
        }
      }
      
      if(use_cov){
        p <- ncol(X) - 3 - 1   # -3 for i,j,t, -1 for intercept removed
      }
      
      Xreg[which(is.na(Xreg), arr.ind=T)] <- 0   # set NAs to zero
      if(writeXreg){  # write out
        save(Xreg, file=nameXreg)
      }
      
      # Bookkeep zero columns
      remove <- which(colSums(Xreg[,1:(S^2)]) <= cutA)    # remove columns with \le cutA 1's in them
      remove <- c(remove, S^2 + which(colSums(Xreg[,S^2 + 1:(L^2)]) <= cutB)) 
      remove <- union(remove, S^2 + L^2 + which(colSums(Xreg[,(S^2 + L^2 +1):ncol(Xreg)]) == 0))
      keep <- setdiff(1:ncol(Xreg), remove)
      
      # Remove duplicate columns
      if(remove_dups & !keep_one_dup){
        stop("Not sure this works for keep_one_dup = F")
        # nonid <- checkcols(Xreg[,keep])  # columns that are duplicates
        # nonid <- sort(unique(c(nonid)))   # vectorize
        # nonid <- nonid_remove <- keep[nonid]    # in terms of original indices
        # remove <- c(remove, nonid)   # modify remove and keep columns
        # keep <- (1:ncol(Xreg))[-remove]
        # 
        # # nonid <- list(keep=NA, remove=nonid)
        # ldups <- nonid_keep <-  nonid_remove <- NA
        
      } else if (!remove_dups) { 
        ldups <- nonid_keep <-  nonid_remove <- NA
        nonid <- NULL
        
      } else if (remove_dups & keep_one_dup) {
        nonid <- checkcols(Xreg[,keep])  # columns that are duplicates
        
        # nonid <- sort(unique(c(nonid)))   # vectorize
        if(!is.null(nonid)){
          nonid <- cbind(keep[nonid[,1]], keep[nonid[,2]])    # in terms of original indices
          
          dups <- sort(unique(c(nonid)))   # columns that are duplicated at some point, in terms of original indices
          a <- as.vector(t(Xreg[,dups]) %*% rnorm(nrow(Xreg)))   # wp1 only same if same number
          ldups <- lapply(sort(unique(a)), function(z) dups[which(a == z)])  # list of all duplicate subsets
          
          remove1 <- unlist(sapply(ldups, function(z) z[-1]))
          remove <- c(remove, remove1)   # modify remove and keep columns
          keep <- (1:ncol(Xreg))[-remove]
          
          nonid_remove <- remove1
          nonid_keep <- sapply(ldups, function(z) z[1])
          nonid <- list(keep=nonid_keep, remove=nonid_remove)
          
          nonid_remove <- sapply(ldups, function(z) z[-1])
        } else {
          ldups <- nonid_keep <-  nonid_remove <- NA
        }
        
        # save(ldups, nonid_keep, nonid_remove, file=file.path(outdir, "nonid_data.RData"))
      } else { stop("Invalid input for dupliate column logicals")}
      
      if(remove_dups & !is.null(nonid)){
        nout <- cbind(unlist(ldups), rep(1:length(ldups), times=sapply(ldups, length)), 0)
        nout <- cbind(nout, nout[,1] <= S^2, nout[,1] > S^2)
        
        ia <- matrix(1:S^2, S,S)
        ib <- matrix(1:L^2 + S^2, L,L)
        cts <- t(sapply(nout[,1], function(z) rbind(which(ia == z, arr.ind=T), which(ib == z, arr.ind=T))))
        nout <- cbind(nout, cts, NA, NA)  
        
        nout[nout[,5]==1, 6:7] <- c(treaties[nout[nout[,5]==1, 6]], treaties[nout[nout[,5]==1, 7]])
        nout[nout[,4]==1, 8:9] <- c(countries[nout[nout[,4]==1, 6]], countries[nout[nout[,4]==1, 7]])
        nout[nout[,5]==0, 7:6] <- NA
        # if(sum(nout[,4]) >= 1){  # A indicator
        #   nout <- cbind(nout, cows$StateAbb[match(nout[,8], cows$CCode)], cows$StateAbb[match(nout[,9], cows$CCode)])
        # } else {nout <- cbind(nout, NA, NA)}
        colnames(nout) <- c("column", "nonid_group", "keep", "Aind", "Bind", "b_i", "b_j", "a_i", "a_j")
        nout <- nout[order(as.numeric(nout[,1])),]
        write.table(nout, file.path(outdir, "nonid_columns.txt"), row.names=F, quote=F)
      } else {
        nout <- NA
      }
      
      if(length(remove) > 0){
        ia <- matrix(1:S^2, S,S)
        ib <- matrix(1:L^2 + S^2, L,L)
        
        ra <- t(sapply(remove[remove <= S^2], function(z) which(ia == z, arr.ind=T)))
        ra <- cbind(ra[,2], ra[,1])   # transpose!
        rb <- t(sapply(remove[remove >  S^2 & remove <=  S^2 + L^2], function(z) which(ib == z, arr.ind=T)))
        
        rb <- cbind(rb, treaties[rb[,1]], treaties[rb[,2]], NA, NA)  #, NA, NA)
        ra <- cbind(ra, NA, NA, countries[ra[,1]], countries[ra[,2]])
        # ra <- cbind(ra, cows$StateNme[match(ra[,5], cows$CCode)], cows$StateNme[match(ra[,6], cows$CCode)])
        
        
        allrem <- rbind(ra,rb)
        allrem <- cbind(allrem, 1*(remove[remove <= S^2 + L^2] %in% nonid_remove))
        colnames(allrem) <- c("i", "j", "b_i", "b_j", "a_i", "a_j", "nonid")  # "treaty1", "treaty2", "cow1", "cow2", "state1", "state2",
        
        rbeta <- remove[remove >  S^2 + L^2] - S^2 - L^2
        if(!is.na(rbeta[1])){
          allrem <- cbind(allrem, NA, NA)
          colnames(allrem) <- c("i", "j", "b_i", "b_j", "a_i", "a_j", "nonid", "ibeta", "coef")
          rc <- cbind(NA, NA, NA, NA, NA, NA, NA, NA, rbeta, (names(X)[!(names(X) %in% c("i","j","t","intercept"))])[rbeta])
        }
        
        write.table(nout, file.path(outdir, "allremoved.txt"), row.names=F, quote=F)
      }
      
      # Full Removal
      remove <- union(remove, unlist(nonid_remove))
      keep <- setdiff(1:ncol(Xreg), remove)      
      
      # Penalty factor
      pf <- rep(1, length(keep))
      if(!penbeta & use_cov){   # add some zeros to not penalize beta coefficients if desired
        pf <- rep(1, length(keep))
        beta_cols <- intersect((S^2 + L^2 + 1):(ncol(Xreg)), keep)
        pf[keep %in% beta_cols] <- 0   # NO penalty on betas when fed to lasso procedure
      }
      
      #### Save inputs 
      save(ldups, nonid_keep, nonid_remove, keep, remove, Xreg, Yreg, 
           Y, D, X, outdir, use_cov, type, penalty, sparsedata, writeXreg, readXreg, 
           nameXreg, seed, ncv, verbose, maxit, ncores, thresh,
           cutA, cutB, remove_dups, keep_one_dup, anames, bnames, 
           nout, allrem, withnet, pf,
           file=file.path(outdir, "prefit2.RData"))
      ####
      
      cat("************* Fit for cutA=cutB=",cutA,"*******************\n\n")
      if(strtrim(cvtype,4) == "year"){ yearsin <- Y$year[indices] } else {yearsin <- NA}
      results <- cv_pr_lasso(Yreg, Xreg[,keep], outdir, penalty, seed=seed, ncv=ncv, verbose=verbose, maxit=maxit, ncores=ncores, thresh=thresh, pf=pf, cvtype=cvtype, years=yearsin)  
      
      cat("Done with while-loop, chose cutA=cutB=", cutA,"\n")
      
      checkfit <- results$nl_full > 25 & all(results$nl_cv > .9*results$nl_full)
      if(!checkfit){
        cutA <- cutB <- cutA + 1
      }
    } # end while loop
    
    
  } else {   # else for withnet
    Xreg <- X[,!(names(X) %in% c("i","j","t","intercept"))]
    Xreg[is.na(Xreg)] <- 0
    rows <- X$i + (X$j-1)*S + (X$t-1)*S*L   # unfolded indices
    indices <- order(rows)
    Xreg <- as.matrix(Xreg[indices, ])  # 1's and zeros, in order
    
    keep <- 1:ncol(Xreg)
    p <- ncol(Xreg)
    remove <- NA
    ldups <- nonid_keep<- nonid_remove <- nout <- allrem <- NA
    
    #### Save inputs 
    save(ldups, nonid_keep, nonid_remove, keep, remove, Xreg, Yreg, 
         Y, D, X, outdir, use_cov, type, penalty, sparsedata, writeXreg, readXreg, 
         nameXreg, seed, ncv, verbose, maxit, ncores, thresh,
         cutA, cutB, global_cov, remove_dups, keep_one_dup, anames, bnames, 
         nout, allrem, withnet, 
         file=file.path(outdir, "prefit2.RData"))
    ####
    
    if(penbeta){
      if(strtrim(cvtype,4) == "year"){ yearsin <- Y$year[indices] } else {yearsin <- NA}
      if(verbose){
        cat("At middle: ncol(Xreg) = ", ncol(Xreg), ", length(keep) = ", length(keep), "\n")
      }
      
      results <- cv_pr_lasso(Yreg, Xreg[,keep], outdir, penalty, seed=seed, ncv=ncv, verbose=verbose, maxit=maxit, ncores=ncores, thresh=thresh, cvtype=cvtype, years=yearsin)  # implict even penalization (pf default input)
    } else {
      r <- glm.fit(x=cbind(1,Xreg[,keep]), y=Yreg, family =binomial())
      results$coefs <- r$coefficients
      results$Yhat <- r$fitted.values
      results$fit <- r
      results$cvm_full <- results$lambda_min <- results$imin <- results$nl_full <- results$nl_cv <- NA
      results$prfull <- simple_pr_auc(results$Yhat, Yreg)
    }
  }
  
  
  
  #### Form coefficient results
  allcoefs <- rep(NA, ncol(Xreg))
  allcoefs[keep] <- results$coefs[-1]
  allcoefs <- c(results$coefs[1], allcoefs)  # add in intercept

  if(withnet){
    A <- t(matrix(allcoefs[1:(S^2) + 1], S,S))   # transposed A for interpretability
    # if(verbose){ cat("dim(A)=", dim(A), "class(A)=", class(A), "\n") }
    colnames(A) <- rownames(A) <- countries
    B <- matrix(allcoefs[1:(L^2) + 1 + S^2], L,L)   # A as in coding, not t(A)
    # if(verbose){ cat("dim(B)=", dim(B), "class(B)=", class(B), "\n") }
    colnames(B) <- rownames(B) <- treaties
  } else { A <- B <- NA}
  
  if(use_cov){
    beta <- matrix(c(allcoefs[1], tail(allcoefs, p)), ncol=1)
    # if(verbose){ cat("length(beta)=", length(beta), ", p=", p, "\n") }
    # if(verbose){ cat("number of rownames for beta =", length(c("intercept", names(X)[!(names(X) %in% c("i","j","t","intercept"))])), "\n") }
    rownames(beta) <- c("intercept", names(X)[!(names(X) %in% c("i","j","t","intercept"))]) 
  } else {beta <- NA}
  
  Yhat <- results$Yhat
  fit <- results$fit
  cvms <- results$cvms
  cvm_full <- results$cvm_full
  lambda_min <- results$lambda_min
  imin <- results$imin
  nl_full <- results$nl_full
  nl_cv <- results$nl_cv
  prfull <- cvm_full[imin]
  if(!withnet & !penbeta){ prfull <- results$prfull}
  save(A,B,beta, results, Yhat, fit, cvms, cvm_full, lambda_min, imin, nl_full, nl_cv, prfull, file=file.path(outdir, "unprocessed_results.RData"))
  ####
  
  
  return(list(A=A, B=B, beta=beta, Yhat=results$Yhat, fit=results$fit, cvms=results$cvms, cvm_full=results$cvm_full, lambda_min=results$lambda_min, imin=results$imin, nl_full=results$nl_full, nl_cv=results$nl_cv, prfull=prfull, cutA=cutA))
}





# P/R cross-validation with automatic lambda selection, standardized X
cv_pr_lasso <- function(Yreg, Xreg, outdir, penalty=1, seed=NA, ncv=10, verbose=F, maxit=1e5, ncores=1, thresh=1e-7, pf=rep(1, ncol(Xreg)), cvtype="random", years=NA)  
{
  if(! is.numeric(penalty)){ stop("Penalized regression method requires numeric penalty value (alpha for glmnet(.))") }
  
  dir.create(outdir, showWarnings = F)
  
  
  #### Perform full fit and extract coefficients for lambda
  if(is.numeric(seed)){ set.seed(seed)}  # set seed
  fitfull <- fit <- glmnet(Xreg, y=as.factor(c(Yreg)), alpha=penalty, family='binomial', intercept=T, lambda.min.ratio=1e-12, standardize=T, maxit=maxit, nlambda=100, thresh=thresh, penalty.factor=pf)
  save(fit, file=file.path(outdir, paste0("fullfit.RData")))
  workedfull <- length(fitfull$lambda) 
  if(verbose){
    cat("full fit done, nlambda", workedfull, "\n")
  }
  lambdas <- fitfull$lambda   # save selected lambda sequence
  ####
  
  
  
  
  
  #### Cross-validate
  lambdas <- sort(lambdas, decreasing = T)  # decreasing sequence
  
  if(cvtype == "random"){
    ones_check <- rep(0, ncv)
    if(is.numeric(seed)){ set.seed(seed)}  # set seed
    count <- 0
    while(any(ones_check == 0) & count < 1e4){   # if there are any partitions without 1's, and make loop finite
      count <- count + 1
      # cat("cv partition count", count, "\n")
      rowscramble <- sample(1:nrow(Xreg), size = nrow(Xreg), replace = F)
      cvs <- rowscramble %% ncv + 1
      ones_check <- sapply(1:ncv, function(z) sum(Yreg[cvs == z]))   # number of 1's in Yreg for each cv partition
    }
  } else if (strtrim(cvtype,4) == "year"){
    ncv <- length(unique(years))
    cvs <- match(years, sort(unique(years)))
    
  } else {
    stop("invalid cross-validation type input cvtype")
  }
    
  cvms <- matrix(NA, length(lambdas), ncv)
  assign("last.warning", NULL, envir = baseenv())   # clear warnings
  nlambda <- rep(NA, ncv)  # number of lambda fit
  
  ####
  save(thresh,seed,ncv,maxit,ncores,penalty,pf, file=file.path(outdir, paste0("cv_pr_lasso_inputs.RData")))
  ####
  
  
  if(ncores==1){
    fitlist <- vector("list", ncv)
    for(i in 1:ncv){
      test <- which(cvs == i)
      train <- (1:nrow(Xreg))[-test]
      
      fitlist[[i]] <- fit <- glmnet(Xreg[train,], y=as.factor(c(Yreg[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-12, standardize=T, thresh = thresh, penalty.factor = pf)
      worked <- nlambda[i] <- length(fit$lambda)
      
      save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
      
      if(verbose){
        cat("\n************************************************\n")
        cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
        print(warnings())
        cat("\n************************************************")
        cat("\n")
        assign("last.warning", NULL, envir = baseenv())
      }
    }
    
  } else if (ncores > 1 ){
    writeLines(c(""), file.path(outdir, "log.txt"))
    
    registerDoMC(cores=round(ncores))
    mcoptions <- list(preschedule=FALSE, set.seed=T)
    fitlist <- foreach(i=1:ncv, .options.multicore=mcoptions, .packages=c("glmnet") ) %dopar% {  # .combine =cbind
      test <- which(cvs == i)
      train <- (1:nrow(Xreg))[-test]
      
      fit <- glmnet(Xreg[train,], y=as.factor(c(Yreg[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-12, standardize=T, thresh = thresh, penalty.factor = pf)
      save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
      worked <- length(fit$lambda)
      sink(file.path(outdir, "log.txt"), append=TRUE)   # write out to log file
      cat("\n************************************************\n")
      cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
      print(warnings())
      cat("\n************************************************")
      cat("\n")
      fit
    }
    
  } else {stop("ncores must be numeric >= 1")}
  
  for(i in 1:ncv){
    test <- which(cvs == i)
    train <- (1:nrow(Xreg))[-test]
    
    fit <- fitlist[[i]]
    #glmnet(Xreg[train,keep], y=as.factor(c(Y[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit)
    worked <- nlambda[i] <- length(fit$lambda)
    
    if(worked > 1){   # if fit worked
      Yhats <- predict(fit, newx=Xreg[test,], type="response")
      
      pr_temp <- rep(NA, ncol(Yhats))
      for(j in 1:ncol(Yhats)){
        # pr_temp[j] <- pr_curve(Yhats[,j], Y[test], n=2000)$auc
        # pr_temp[j] <- roc_curve(Yhats[,j], Y[test], n=min(length(Y[test]), 2500))$auc
        # pr_temp[j] <- simple_roc_auc(Yhats[,j], Y[test])
        pr_temp[j] <- simple_pr_auc(Yhats[, j], Yreg[test])
        # pr_temp[j] <- auc(roc(Yhats[,j], Y[test]))
      }
      
      cvms[match(fit$lambda, lambdas), i] <- pr_temp   # save lambda values
      # save(fit, pr_temp, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
    }
    if(verbose){
      cat("processed cv", i, "of", ncv, ", nlambda", worked, "\n")
    }
  }
  
  # Calculate minimum lambda value
  mean_cvms <- apply(cvms, 1, mean, na.rm=T)
  imin <- which(mean_cvms == max(mean_cvms, na.rm=T))[1]
  lambda_min <- lambdas[imin]
  ####
  
  
  #### Calculate results from full fit
  fit <- fitfull
  if(workedfull > 1){
    Yhats <- predict(fitfull, newx=Xreg, type="response")
    cvm_full <- rep(NA, ncol(Yhats))
    for(j in 1:ncol(Yhats)){
      cvm_full[j] <- simple_pr_auc(Yhats[,j], Yreg) # pr_curve(Yhats[,j], Y, n=min(length(Y), 2500))$auc
    }
    Yhat <- predict(fit, newx=Xreg, type="response", s=lambda_min)
    coefs <- coef(fit, s=lambda_min)
    
  } else {
    fit <- Yhat <- coefs <- cvm_full <- NA
  }
  ####
  
  
  return(list(coefs=coefs, Yhat=Yhat, fit=fitfull, Yreg=Yreg, Xreg=Xreg, cvms=cvms, cvm_full=cvm_full, lambda_min=lambda_min, imin=imin, nl_full=workedfull, nl_cv=nlambda))
}



# Process generic results
process_generic <- function(wd=getwd(), selectiveInfLib="~/Dropbox/BiTEN/selectiveInferenceFM", calcses=T, withnet=T)
{
  setwd(wd)
  library("RColorBrewer")
  
  nonet <- !withnet
  outdir <- wd
  penbeta <- T
  
  # #### Inputs
  # final_years <- 1993  #setdiff(1992:1999, c(1993)) #setdiff(1960:1999, c(1960, 1974, 1977, 1994)) # last year of fit window
  # window <- 20
  # outstub <- "cv_pr_split_window"
  # # nonet <- T
  # penbeta <- T
  # ps <- "_lag3"
  # calcses <- F  # calculate standard errors??
  # ####
  
  # 
  # for(nonet in c(T,F)){
  #   
  #   for(yw in final_years){
  # cat("*********************  \n starting (end-window) year", yw, "\n********************* \n\n" )
  
  # outdir <- paste0(outstub, yw-window+1, "_", yw)  #, ps,"_penbeta", as.numeric(penbeta))
  # outdir <- paste0(outstub, yw-window+1, "_", yw, ps,"_penbeta", as.numeric(penbeta))
  # if(nonet){outdir <- paste0(outdir, "_nonet")}
  
  #### Read results from cv_pr*2.R,    
  #
  #### NOTE: A is transposed version for interpretation (NOT as in fitting)  ####
  #
  load(file.path(outdir, "prefit2.RData"))
  # outdir <- paste0(outstub, yw-window+1, "_", yw, ps,"_penbeta", as.numeric(penbeta))
  # if(nonet){outdir <- paste0(outdir, "_nonet")}
  
  load(file.path(outdir, "unprocessed_results.RData"))
  cat("PR for full fit: ", prfull, "\n")
  write.table(round(prfull, 4), file.path(outdir, "cvm_full.txt"), row.names=F, col.names=T, quote=F)
  px <- ncol(X) - 4  # length of beta without intercept
  
  
  if(!nonet){
    S <- ncol(A)
    L <- ncol(B)
    write.table(A, file.path(outdir, "A.txt"), row.names=T, col.names = T, quote=F)   # NOTE TRANSPOSE
    write.table(B, file.path(outdir, "B.txt"), row.names=T, col.names = T, quote=F)
    write.table(beta, file.path(outdir, "beta.txt"), row.names=T, col.names = T, quote=F)
    writeLines("The version of A saved in this directory is the transposed version as in the BLIN model for interpretability,
               NOT the version of A used in the code (which is t(A)). ", con=file.path(outdir, "A_transpose_note.txt"))
    write.table(cvms, file.path(outdir,"cvms.txt"), row.names=F, col.names=F, quote=F)
    write.table(max(cutA,cutB), file.path(outdir,"cA.txt"), row.names=F, col.names=F, quote=F)
  } else {
    write.table(beta, file.path(outdir, "beta.txt"), row.names=T, col.names = T, quote=F)
    write.table(cvms, file.path(outdir,"cvms.txt"), row.names=F, col.names=F, quote=F)
  }
  ####
  
  
  #### Some basic network stuff
  ncv <- ncol(results$cvms)
  if(!nonet){
    Arem <- table(colSums(Xreg[,1:(S^2)]))
    Arem <- Arem[as.numeric(names(Arem)) <= max(cutA, cutB)]
    if(length(Arem)  <= max(cutA, cutB)){ Arem[length(Arem):max(cutA, cutB) + 1] <- 0}
    names(Arem) <- paste0("A", 0:max(cutA, cutB))
    A_na_frac <- mean(is.na(A))
    Adensity <- sum(A != 0, na.rm=T)/S^2
    
    Brem <- table(colSums(Xreg[,S^2 + 1:(L^2)]))
    Brem <- Brem[as.numeric(names(Brem)) <= max(cutA, cutB)]
    if(length(Brem)  <= max(cutA, cutB)){ Brem[length(Brem):max(cutA, cutB) + 1] <- 0}
    names(Brem) <- paste0("B", 0:max(cutA, cutB))
    B_na_frac <- mean(is.na(B))
    Bdensity <- sum(B != 0, na.rm=T)/L^2
    
    outtab <- rbind(Arem, Brem)
    outtab <- cbind(outtab, round(c(A_na_frac, B_na_frac), 5))
    outtab <- cbind(outtab, round(c(Adensity, Bdensity), 5))
    rownames(outtab) <- c("A", "B")
    colnames(outtab) <-  c(paste0("col_sum_", 0:max(cutA, cutB)), "frac_NA", "frac_nonzero")
    
    write.table(outtab, file.path(outdir, "ABsummary.txt"), quote=F)
  }
  ####
  
  
  #### Plots for lambda
  if(penbeta | !nonet){
    colors <- brewer.pal(8, "Set2")
    Lout <- data.frame(cbind(results$fit$lambda, 1*((1:length(results$fit$lambda)) == imin), apply(results$cvms, 1, function(z) sum(!is.na(z)))))
    names(Lout) <- c("lambda", "min", "nfits")
    write.table(Lout, file.path(outdir, "lambdas.txt"), row.names=F, col.names=T, quote=F)
    
    pdf(file.path(outdir, "nlambda.pdf"), height=4,width=4)
    par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    barplot(apply(results$cvms, 2, function(z) sum(!is.na(z))), names.arg =1:ncv, xlab="CV", ylab="nlambda", col=colors[3], border = "white")
    dev.off()
    
    pdf(file.path(outdir, "AUC_cv.pdf"), height=4,width=4)
    par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    boxplot(t(results$cvms[nrow(results$cvms):1,]), range=0, xaxt="n", ylab="Out-of-sample AUC P/R", xlab=expression("log("~lambda~")"))
    axis(1, at=1:nrow(results$cvms), round(log(results$fit$lambda)[nrow(results$cvms):1],1))
    rmin <- which(length(results$fit$lambda):1 == results$imin)  # reversed index
    abline(v=rmin, col="gray50", lwd=2.5)
    dev.off()
    
    pdf(file.path(outdir, "AUC_cv2.pdf"), height=4,width=4)
    par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    plot(log(results$fit$lambda), results$cvm_full, ylab="AUC P/R", xlab=expression("log("~lambda~")"))
    lines(log(results$fit$lambda), apply(results$cvms, 1, mean, na.rm=T), col="red", lwd=2)
    axis(1, at=1:nrow(results$cvms), round(log(results$fit$lambda)[nrow(results$cvms):1],1))
    abline(v=log(results$lambda_min), col="blue", lwd=2)
    legend("topright", c("Full dataset", "Mean AUC in CV", expression("Min CV"~lambda)), col=c("black", "red", "blue"), lwd=2, lty=c(NA,1,1), pch=c(1,NA,NA), bty="n", cex=.6)
    dev.off()
  }
  ####
  
  
  
  
  #### Coeff plots
  A[is.na(A)] <- 0
  B[is.na(B)] <- 0
  if(!nonet){
    pdf(file.path(outdir, "Acoeff.pdf"), height=4,width=4)
    par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    plot(colSums(Xreg[,1:(S^2)]), c(A), xlab="Column sums for A", ylab="Coefficient estimates of A")
    dev.off()
    
    ca <- sort(unique(colSums(Xreg[,1:(S^2)])))
    i_colsums <- lapply(ca, function(z) which(colSums(Xreg[,1:(S^2)]) == z))
    a_nonzero <- sapply(i_colsums, function(z) sum(c(A)[z] != 0, na.rm=T)/length(z))
    pdf(file.path(outdir, "Acoeff_frac.pdf"), height=4,width=4)
    par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    plot(ca, a_nonzero, xlab="Column sums for A", ylab="Fraction of nonzero entries in A", ylim=c(0,.2))
    dev.off()
    
    pdf(file.path(outdir, "Bcoeff.pdf"), height=4,width=4)
    par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    plot(colSums(Xreg[,S^2 + 1:(L^2)]), c(B), xlab="Column sums for B", ylab="Coefficient estimates of B")
    dev.off()
    # plot(table(colSums(Xreg[,intersect(1:(S^2),keep)])), xlab="Column sums for A", ylab="Frequency")
    # plot(table(colSums(Xreg[,intersect(keep, (S^2) + 1:(L^2))])), xlab="Column sums for B", ylab="Frequency")
    
    cb <- sort(unique(colSums(Xreg[,S^2 + 1:(L^2)])))
    i_colsums <- lapply(cb, function(z) which(colSums(Xreg[,S^2 + 1:(L^2)]) == z))
    b_nonzero <- sapply(i_colsums, function(z) sum(c(B)[z] != 0, na.rm=T)/length(z))
    pdf(file.path(outdir, "Bcoeff_frac.pdf"), height=4,width=4)
    par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    plot(cb, b_nonzero, xlab="Column sums for B", ylab="Fraction of nonzero entries in B", ylim=c(0,.2))
    dev.off()
  }
  ####
  
  
  
  
  #### Standard errors using selectiveInference package
  library("selectiveInference", lib.loc = selectiveInfLib)  # after Frank's edits
  # library("selectiveInference", lib.loc = )  # after Frank's edits
  # Xreg <- results$Xreg
  Lcheck <- lambda_min*nrow(Xreg) < max(fit$lambda) & lambda_min*nrow(Xreg) > min(fit$lambda)
  
  if(penbeta & !nonet){
    # theta <- coef(results$fit)[,imin]   # with intercept, but fixedLassoInf knows to remove first entry (in fact intercept required)
    theta <- c(beta[1], c(t(A)), c(B), c(beta[-1]))
  } else { 
    theta <- beta 
  }
  
  
  remove0 <- which(theta[-1] == 0)
  keep0 <- setdiff(1:(length(theta) - 1), remove0)
  keep1 <- intersect(keep, keep0)
  remove1 <- setdiff(1:(length(theta) - 1), keep1)
  
  Xreg1 <- Xreg[,keep1]
  theta1 <- c(theta[1], theta[-1][keep1])
  preds <- length(theta1) > 1
  
  
  if(preds & calcses){
    out <- fixedLassoInf(x=(Xreg1), y=Yreg, beta=theta1, lambda=lambda_min*nrow(Xreg1), alpha=.05, family = "binomial", tol.beta=1e-10)  #lambda_min*nrow(Xreg1)
    cat("Fraction of nonzero coefficients that are significant by p-value:  ", mean(out$pv <= .05), "\n")
    
    
    if(sum(range(Xreg[,keep1] - Xreg1)) != 0){stop("Something went wrong with indexing in terms of original indices")}
    
    se <- rep(NA, ncol(Xreg))
    se[keep1] <- out$sd
    
    if(!nonet){
      sea <- c(t(matrix(se[1:(S^2)], S, S)))   # transpose A results
      seb <- se[S^2 + 1:(L^2)]
      sebeta <- c(0, tail(se, px))
      
      save(out,  remove1, keep1, se, sea, seb, sebeta, theta, file=file.path(outdir, "se_from_selectiveInference.RData"))
      
      r <- range(c(
        range(which(is.na(sea)) - which(is.na(c((A))) | c((A)) == 0)),
        range(which(is.na(seb)) - which(is.na(c(B)) | c(B) == 0)),
        range(which(is.na(sebeta)) - which(is.na(c(beta)) | c(beta) == 0))))   # check if NAs and zeros line up... they do!
      cat("Match 0s and NA check, should be zero: ", r,"\n")
    } else {
      sebeta <- c(0, tail(se, px))
      save(out,  remove1, keep1, se, sebeta, theta, file=file.path(outdir, "se_from_selectiveInference.RData"))
    }
  }
  
  print(warnings())
  assign("last.warning", NULL, envir = baseenv())  # clear warnings
  ####
  
  
  
  
  #### Shape into useful output
  countries <- anames
  treaties <- bnames
  if(preds & calcses){
    beta_names <- rownames(beta)
    alpha <- .05
    px <- length(beta) - 1
    
    if(!nonet){
      ciA <- data.frame(matrix(NA, S^2, 13))   ;  names(ciA) <- c("i", "j", "a_i", "a_j", "estimate", "lower", "upper", "pval", "signif", "se", "alpha", "nonid_keep", "nonid_remove")
      ciB <- data.frame(matrix(NA, L^2, 13))   ;  names(ciB) <- c("i", "j", "b_i", "b_j", "estimate", "lower", "upper", "pval", "signif", "se", "alpha", "nonid_keep", "nonid_remove")
      cibeta <- data.frame(matrix(NA, px + 1, 8))   ;  names(cibeta) <- c("covariate", "estimate", "lower", "upper", "pval", "signif", "se", "alpha")
      
      ciA[,1:2] <- cbind(rep(1:(S), times=S), rep(1:(S), each=S))   # columnwise unfolding
      ciB[,1:2] <- cbind(rep(1:(L), times=L), rep(1:(L), each=L))
      cibeta$covariate <- beta_names
      ciA$cow1 <- countries[ciA$i]  ;  ciA$cow2 <- countries[ciA$j]
      ciB$treaty1 <- treaties[ciB$i]  ;  ciB$treaty2 <- treaties[ciB$j]
      ciA$alpha <- ciB$alpha <- cibeta$alpha <- alpha
      
      ciA$estimate <- c(A)
      ciB$estimate <- c(B)
      cibeta$estimate <- c(beta)
      
      ciA$se <- sea
      ciB$se <- seb
      cibeta$se <- sebeta 
      
      ci <- matrix(NA, ncol(Xreg), 2)
      ci[keep1,] <- out$ci
      ciA[,6:7] <- ci[c(t(matrix(1:(S^2), S, S))),]   # transposed!
      ciB[,6:7] <- ci[S^2 + 1:(L^2),]
      cibeta[,3:4] <- rbind(0, ci[tail(1:ncol(Xreg), px),])
      
      pv <- rep(NA, ncol(Xreg))
      pv[keep1] <- out$pv
      ciA$pval <- pv[c(t(matrix(1:(S^2), S, S)))]   # transposed!
      ciB$pval <- pv[S^2 + 1:(L^2)]
      cibeta$pval <- c(0, tail(pv, px))
      
      cisignif <-  1*!(ci[,1] < 0 & ci[,2] > 0 )   # significance based on confidence intervals
      cat("Fraction of pvalue < alpha and non-0-containing confidence intervals agree is : ", mean(cisignif == 1*(pv < alpha), na.rm=T), "\n")
      cisignif2 <- 1*(pv < alpha)  # significance based on pvalue
      cisignif <- cisignif*cisignif2   # only significant if both significant
      
      ciA$signif <- cisignif[c(t(matrix(1:(S^2), S, S)))]   
      ciB$signif <- cisignif[S^2 + 1:(L^2)]
      cibeta$signif <- c(1, tail(cisignif, px))
      
      ciA$nonid_keep <- ciA$nonid_remove <- ciB$nonid_keep <- ciB$nonid_remove <- 0
      if(length(nonid_keep) > 0 & !is.na(nonid_keep[1])){
        ia <- t(matrix(1:(S^2), S, S))
        ib <- matrix(1:(L^2) + S^2, L, L)
        nonid_remove2 <- unlist(nonid_remove)
        
        ka <- sapply(nonid_keep[nonid_keep <= S^2], function(z) which(z == ia))  # transposed (interpretable) A index
        if(length(ka) > 0){
          ciA$nonid_keep[ka] <- 1
          ra <- sapply(nonid_remove2[nonid_remove2 <= S^2], function(z) which(z == ia))  
          ciA$nonid_remove[ra] <- 1
        }
        
        kb <- sapply(nonid_keep[nonid_keep > S^2], function(z) which(z == ib))  
        if(length(kb) > 0){
          ciB$nonid_keep[kb] <- 1
          rb <- sapply(nonid_remove2[nonid_remove2 > S^2], function(z) which(z == ib))  
          ciB$nonid_remove[rb] <- 1
        }
      } 
      
      write.table(ciA, file.path(outdir, "ciA.txt"))
      write.table(ciB, file.path(outdir, "ciB.txt"))
      write.table(cibeta, file.path(outdir, "cibeta.txt"))
    } else {
      cibeta <- data.frame(matrix(NA, px + 1, 8))   ;  names(cibeta) <- c("covariate", "estimate", "lower", "upper", "pval", "signif", "se", "alpha")
      
      cibeta$covariate <- beta_names
      cibeta$alpha <- alpha
      cibeta$estimate <- c(beta)
      cibeta$se <- sebeta 
      
      ci <- matrix(NA, ncol(Xreg), 2)
      ci[keep1,] <- out$ci
      cibeta[,3:4] <- rbind(0, ci[tail(1:ncol(Xreg), px),])
      
      pv <- rep(NA, ncol(Xreg))
      pv[keep1] <- out$pv
      cibeta$pval <- c(0, tail(pv, px))
      
      cisignif <-  1*!(ci[,1] < 0 & ci[,2] > 0 )   # significance based on confidence intervals
      cat("Fraction of pvalue < alpha and non-0-containing confidence intervals agree is : ", mean(cisignif == 1*(pv < alpha), na.rm=T), "\n")
      cisignif2 <- 1*(pv < alpha)  # significance based on pvalue
      cisignif <- cisignif*cisignif2   # only significant if both significant
      cibeta$signif <- c(1, tail(cisignif, px))
      
      write.table(cibeta, file.path(outdir, "cibeta.txt"))
    }
  }
  ####
  
  
  
  
}


AB_plots <- function(wd=getwd(), withnet=T)
{
  setwd(wd)
  library("RColorBrewer")
  colors <- brewer.pal(8, "Set2")
  nonet <- !withnet
  outdir <- wd
  penbeta <- T
  
  load(file.path(outdir, "prefit2.RData"))
  outdir <- wd
  
  load(file.path(outdir, "unprocessed_results.RData"))
  outdir <- wd
  
  ciA <- read.table(file.path(outdir, "ciA.txt"))
  S <- sqrt(nrow(ciA))
  ciB <- read.table(file.path(outdir, "ciB.txt"))
  L <- sqrt(nrow(ciB))
  
  ia <- c(t(matrix(1:S^2, S, S)))   # transposed A indices
  csa <- colSums(Xreg[,ia])
  csb <- colSums(Xreg[,S^2 + 1:(L^2)])
  
  ia0 <- which(ciA$estimate != 0)
  ia1 <- which(ciA$signif == 1)
  ib0 <- which(ciB$estimate != 0)
  ib1 <- which(ciB$signif == 1)
  
  
  pdf(file.path(outdir, "Acoeff_ests.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(csa[ia0], ciA$estimate[ia0], ylab="estimate", xlab="Number of opportunities for inference", col="gray80", lwd=1.5)
  points(csa[ia1], ciA$estimate[ia1], col=colors[1], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  pdf(file.path(outdir, "Acoeff_scaled.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(csa[ia0], abs(ciA$estimate[ia0]*csa[ia0]), ylab="Absolute total influence", xlab="Number of opportunities for inference", col="gray80", lwd=1.5)
  points(csa[ia1], abs(ciA$estimate[ia1]*csa[ia1]), col=colors[1], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  pdf(file.path(outdir, "Acoeff_pval.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(log10(ciA$pval)[ia0], ciA$estimate[ia0], ylab="estimate", xlab="log10( p-value )", col="gray80", lwd=1.5)
  points(log10(ciA$pval)[ia1], ciA$estimate[ia1], col=colors[1], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  pdf(file.path(outdir, "Acoeff_scaled_pval.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(log10(ciA$pval)[ia0], abs(ciA$estimate[ia0]*csa[ia0]), ylab="Absolute total influence", xlab="log10( p-value )", col="gray80", lwd=1.5)
  points(log10(ciA$pval)[ia1], abs(ciA$estimate[ia1]*csa[ia1]), col=colors[1], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  
  
  # plot(csb[ib0], ciB$estimate[ib0], ylab="estimate", xlab="Number of opportunities for inference", col="gray80", lwd=1.5)
  # points(csb[ib1], ciB$estimate[ib1], col=colors[2], lwd=2)
  # 
  # plot(csb[ib0], abs(ciB$estimate[ib0]*csb[ib0]), ylab="Absolute total influence", xlab="Number of opportunities for inference", col="gray80", lwd=1.5)
  # points(csb[ib1], abs(ciB$estimate[ib1]*csb[ib1]), col=colors[2], lwd=2)
  
  pdf(file.path(outdir, "Acoeff_withCIs.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  jitter <- rnorm(length(ia1),0,.1)
  lowerbounds <- ciA$lower[ia1]
  upperbounds <- ciA$upper[ia1]
  upperbounds[is.infinite(upperbounds)] <- 1e5
  lowerbounds[is.infinite(lowerbounds)] <- -1e5
  plot(csa[ia1] + jitter, ciA$estimate[ia1], col=paste0(colors[1], "75"), lwd=2, 
       ylab="estimate", xlab="Number of opportunities for inference", ylim=10*range(ciA$estimate[ciA$signif==1], na.rm=T))
  
  arrows(csa[ia1] + jitter, lowerbounds, y1=upperbounds, angle=90, code=3, length=.025, col=colors[1]) #, code=3, angle=20, lwd=2)
  # arrows(10, -5, y1=5, angle=90, code=3, length=.025, col=colors[1])
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  m <- length(ia1)
  iar <- sample(which(ciA$estimate != 0), m)   # random A entries
  
  
  pdf(file.path(outdir, "Acoeff_random_withCIs.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  jitter <- rnorm(length(ia1),0,.1)
  lowerbounds <- ciA$lower[iar]
  upperbounds <- ciA$upper[iar]
  upperbounds[is.infinite(upperbounds)] <- 1e5
  lowerbounds[is.infinite(lowerbounds)] <- -1e5
  plot(csa[iar] + jitter, ciA$estimate[iar], col=paste0(colors[1], "75"), lwd=2, 
       ylab="estimate", xlab="Number of opportunities for inference", ylim=10*range(ciA$estimate[ciA$signif==1], na.rm=T))
  
  arrows(csa[iar] + jitter, lowerbounds, y1=upperbounds, angle=90, code=3, length=.025, col=colors[1]) #, code=3, angle=20, lwd=2)
  # arrows(10, -5, y1=5, angle=90, code=3, length=.025, col=colors[1])
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  
  
  pdf(file.path(outdir, "Bcoeff_ests.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(csb[ib0], ciB$estimate[ib0], ylab="estimate", xlab="Number of opportunities for inference", col="gray80", lwd=1.5)
  points(csb[ib1], ciB$estimate[ib1], col=colors[2], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  pdf(file.path(outdir, "Bcoeff_scaled.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(csb[ib0], abs(ciB$estimate[ib0]*csb[ib0]), ylab="Absolute total influence", xlab="Number of opportunities for inference", col="gray80", lwd=1.5)
  points(csb[ib1], abs(ciB$estimate[ib1]*csb[ib1]), col=colors[2], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  pdf(file.path(outdir, "Bcoeff_pval.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(log10(ciB$pval)[ib0], ciB$estimate[ib0], ylab="estimate", xlab="log10( p-value )", col="gray80", lwd=1.5)
  points(log10(ciB$pval)[ib1], ciB$estimate[ib1], col=colors[2], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  pdf(file.path(outdir, "Bcoeff_scaled_pval.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  plot(log10(ciB$pval)[ib0], abs(ciB$estimate[ib0]*csb[ib0]), ylab="Absolute total influence", xlab="log10( p-value )", col="gray80", lwd=1.5)
  points(log10(ciB$pval)[ib1], abs(ciB$estimate[ib1]*csb[ib1]), col=colors[2], lwd=2)
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  pdf(file.path(outdir, "Bcoeff_withCIs.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  jitter <- rnorm(length(ib1),0,.1)
  lowerbounds <- ciB$lower[ib1]
  upperbounds <- ciB$upper[ib1]
  upperbounds[is.infinite(upperbounds)] <- 1e5
  lowerbounds[is.infinite(lowerbounds)] <- -1e5
  plot(csb[ib1] + jitter, ciB$estimate[ib1], col=paste0(colors[2], "75"), lwd=2, 
       ylab="estimate", xlab="Number of opportunities for inference", ylim=10*range(ciB$estimate[ciB$signif==1], na.rm=T))
  
  arrows(csb[ib1] + jitter, lowerbounds, y1=upperbounds, angle=90, code=3, length=.025, col=colors[2]) #, code=3, angle=20, lwd=2)
  # arrows(10, -5, y1=5, angle=90, code=3, length=.025, col=colors[1])
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  m <- length(ib1)
  ibr <- sample(which(ciB$estimate != 0), m)   # random B entries
  
  
  pdf(file.path(outdir, "Bcoeff_random_withCIs.pdf"), height=4,width=6)
  par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  jitter <- rnorm(length(ib1),0,.1)
  lowerbounds <- ciB$lower[ibr]
  upperbounds <- ciB$upper[ibr]
  upperbounds[is.infinite(upperbounds)] <- 1e5
  lowerbounds[is.infinite(lowerbounds)] <- -1e5
  plot(csb[ibr] + jitter, ciB$estimate[ibr], col=paste0(colors[2], "75"), lwd=2, 
       ylab="estimate", xlab="Number of opportunities for inference", ylim=10*range(ciB$estimate[ciB$signif==1], na.rm=T))
  
  arrows(csb[ibr] + jitter, lowerbounds, y1=upperbounds, angle=90, code=3, length=.025, col=colors[2]) #, code=3, angle=20, lwd=2)
  # arrows(10, -5, y1=5, angle=90, code=3, length=.025, col=colors[1])
  abline(h=0, col="black", lwd=1)
  dev.off()
  
  
  ciA$num_inf <- csa
  ciB$num_inf <- csb
  ciA$total_influence <- csa*ciA$estimate
  ciB$total_influence <- csb*ciB$estimate
  
  write.table(ciA, file.path(outdir, "ciA.txt"))
  write.table(ciB, file.path(outdir, "ciB.txt"))
  
  mns <- mean(ciA$estimate[-ia1] > ciA$lower[-ia1] & ciA$estimate[-ia1] < ciA$upper[-ia1], na.rm=T)
  ms <- mean(ciA$estimate[ia1] > ciA$lower[ia1] & ciA$estimate[ia1] < ciA$upper[ia1], na.rm=T)
  absinf <- csa*ciA$estimate
  
  cicheck <- cbind(mns, ms)
  colnames(cicheck) <- c("mean_nonsignif_in_CI", "mean_signif_in_CI")
  return(list(cicheck=cicheck, absinf=absinf))
}


# Post-process BLIN fits to get CIs
fitzeros <- function(wd=getwd(), type1=T, nit=25, verbose=TRUE)
{
  setwd(wd)
  # require("RColorBrewer")
  # require("MatrixModels")    # glm4 
  require("glmnet")    # glm4 
  
  nonet <- F
  # outdir <- wd
  
  load(file.path(wd, "prefit2.RData"))
  load(file.path(wd, "unprocessed_results.RData"))
  
  ciA <- read.table(file.path(wd, "ciA.txt"))
  S <- sqrt(nrow(ciA))
  ciB <- read.table(file.path(wd, "ciB.txt"))
  L <- sqrt(nrow(ciB))
  cibeta <- read.table(file.path(wd, "cibeta.txt"))
  px <- nrow(cibeta)
  
  ia <- c(t(matrix(1:S^2, S, S)))   # transposed A indices
  
  A <- t(matrix(ciA$estimate, S, S))
  B <- matrix(ciB$estimate, L, L)
  beta <- cibeta$estimate
  
  theta <- c(beta[1], c(t(A)), c(B), c(beta[-1]))
  
  remove0 <- which(theta[-1] == 0)
  keep0 <- setdiff(1:(length(theta) - 1), remove0)
  keep1 <- intersect(keep, keep0)
  remove1 <- setdiff(1:(length(theta) - 1), keep1)
  
  Xreg1 <- Xreg[,keep1]
  theta1 <- c(theta[-1][keep1])
  
  if(type1){
    fit1 <- glm(Yreg ~ as.matrix(Xreg1), family="binomial", start=c(theta[1], theta1), model=FALSE, epsilon = 1e-6, maxit=nit)
  } else {
    # add intercept
    Xreg2 <- sparseMatrix(i=1:nrow(Xreg1),j=rep(1, nrow(Xreg1)), x=1, dims=c(nrow(Xreg1), 1+ncol(Xreg1)))
    Xreg2[, 1+1:ncol(Xreg1)] <- Xreg1
    # fit1 <- glm.fit(Xreg2, Yreg, family=binomial(link=logit), start=c(theta[1], theta1))
    # Xreg2 <- model.Matrix(Xreg1, sparse=TRUE)
    # Xreg2 <- sparse.model.matrix(Xreg1)
    # fit1 <- glm4(Yreg ~ Xreg2, family="binomial", sparse=TRUE)
    
    fit1 <- glmnet(Xreg1, as.factor(Yreg), family="binomial")
    fit1 <- speedglm(x=Xreg1, Y=c(Yreg), family=binomial(), sparse=TRUE, formula=Yreg~Xreg1+1)
  }
  
  if(verbose){
    cat("Warnings: \n")
    print(warnings())
  }
  # 
  
  return(fit1)
}


##############################
###  Prediction functions  ###
##############################

# These functions make predictions based on the results of previous fits



# roll forward prediction from some initial time, including random realizations of signatures
#  t0 is first year of PREDICTION
#  stochastic prediction requires number of simulations
roll_forward_predict <- function(t0, nsims, A, B, beta, Y, D, X, lag, region_income, tfinal=NULL, response="ratification_year", seed0=1, verbose=F, NApairs=NULL, write_interval=ceiling(nsims/10), outdir=getwd(), filename=NULL, datadir=NULL)  # countries=NULL, treaties=NULL, 
{
  if(is.null(tfinal)){ tfinal <- max(Y$t)}
  trange <- t0:tfinal
  year_range <- unique(Y$year[Y$t == t0] ) : unique(Y$year[Y$t == tfinal] )
  
  # Initialize arrays
  Ynew <- Y[Y$t >= t0 & Y$t <= tfinal,]
  Yold <- Y[Y$t < t0,]
  Dnew <- D[D$t >= t0 & D$t <= tfinal,]
  Xnew <- X[X$t >= t0 & X$t <= tfinal,]
  S <- max(Dnew$i)  ;   L <- max(Dnew$j)    
  Yhat <- phat <- array(0, c(S,L,length(trange),nsims))
  countries <- sapply(1:S, function(z) unique(Dnew$cowcode[Dnew$i == z]))   # all countries of interest in future
  treaties <- sapply(1:L, function(z) unique(Dnew$treaty[Dnew$j == z]))    # all treaties of interest in future
  dimnames(Yhat)[[1]] <- dimnames(phat)[[1]] <- countries
  dimnames(Yhat)[[2]] <- dimnames(phat)[[2]] <- treaties
  dimnames(Yhat)[[3]] <- dimnames(phat)[[3]] <- year_range
  region_income$i <- Y$i[match(region_income$cowcode, Y$cowcode)]   # add column for i in region_income indicators 
  region_income <- region_income[order(region_income$i),]    # reorder
  
  # Read in big X if X is not big
  if(nrow(X) <  S*L*length(trange)){
    X <- build_big_X(t0, X, S=S, L=L, tfinal=max(X$t), readfile=T)
  }
  
  # Update Yold to include all past actions and possibile signings
  if(min(Yold$year) > 1950){
    if(is.null(datadir)){
      datadir <- "~/Dropbox/BiTEN"
    }
    data <- read_env_treaties(datadir, lag=1, write=F, readfile=T)
    Ytemp <- data$Y    # old data
    rm(data)
    veryoldrats <- Ytemp[Ytemp$year < min(Yold$year),]
    veryoldrats$t <-  veryoldrats$t - max(veryoldrats$t)   # recode t
    veryoldrats$j <- Y$j[match(veryoldrats$treaty, Y$treaty)]  # recode j
    veryoldrats$i <- Y$i[match(veryoldrats$cowcode, Y$cowcode)]  # recode i
    Yold <- rbind(Yold, as.matrix(veryoldrats))
    
    newrats <- Ytemp[Ytemp$year > max(Yold$year) & Ytemp$year <= max(Ynew$year),]   # new ratifications that are not available in current dataset
    newrats$t <-  Y$t[match(newrats$year, Y$year)]  # recode t
    newrats$j <- Y$j[match(newrats$treaty, Y$treaty)]  # recode j
    newrats$i <- Y$i[match(newrats$cowcode, Y$cowcode)]  # recode i
    newrats_string <- apply( newrats[, c("cowcode", "j", "t" )], 1, paste0, collapse=",")
    oldrats_string <- apply( Ynew[, c("cowcode", "j", "t")], 1, paste0, collapse=",")
    suppressWarnings( addrats <- t(sapply(setdiff(newrats_string, oldrats_string), function(z) as.numeric(strsplit(z, ",")[[1]]))) )
    keep_add <- apply(addrats[,2:3], 1, function(z) !any(is.na(z)))  # remove j and t indices that are NA
    keep_add_string <- apply(addrats[keep_add,], 1, paste0, collapse=",")
    
    newrats1 <- newrats[match(keep_add_string, newrats_string ), ]   # unaccounted for signatures for X
    Yold <- rbind(Yold, as.matrix(newrats1))
    remove_rows <- unique( which(is.na(Yold[,c("j", "t")]), arr.ind=T)[,1] )   # remove NA rows FOR J AND T
    Yold <- Yold[-remove_rows,]
  }
  
  
  # Find initial possible signatures in t0 and those i,j pairs that already signed
  possibles <- as.matrix(unique(Ynew[Ynew$t == t0, c("i", "j")]))   # all possible signatures in current year
  already_signed <- as.matrix(unique(Yold[Yold$t < t0 & Yold$ratification_year == 1, c("i", "j")]))   # already signed country/treaty pairs
  
  # Find i,j pairs that enter after t0 and the corresponding time when they do so
  match_matrix <- matrix(rnorm(S*L), S, L)   
  match_matrix[rbind( possibles, already_signed)] <- 0   # unsigned and signed country/treaty pairs
  leftovers <- which(match_matrix != 0, arr.ind=T)   # not signed AND not possible in year t0
  suppressWarnings( tmin <- sapply(1:nrow(leftovers), function(z) min(Ynew$t[Ynew$i == leftovers[z,1] & Ynew$j == leftovers[z,2]]) ) )
  tmin[which(is.infinite(tmin))] <- NA   # set infinte to NAs
  late_entry <- (tmin > t0 )*tmin
  late_entry[which(late_entry <= t0)] <- NA  #NA out non-late entries
  late_entry1 <- cbind(leftovers, late_entry)[!is.na(late_entry),]
  late_entry <- late_entry1
  
  # Find i,j pairs that exit the dataset early
  suppressWarnings( maxt <- sapply(1:nrow(possibles), function(z) max(Ynew$t[Ynew$i == possibles[z,1] & Ynew$j == possibles[z,1]])) )
  suppressWarnings( ratt <- sapply(1:nrow(possibles), function(z) Ynew$t[Ynew$ratification_year == 1 & Ynew$i == possibles[z,1] & Ynew$j == possibles[z,1]]) )
  maxt[is.infinite(maxt)] <- NA
  # ratt[is.infinite(ratt)] <- NA
  ratt1 <- rep(NA, length(ratt))
  ratt1[which(sapply(ratt, length) > 0)] <- unlist(ratt[which(sapply(ratt, length) > 0)])
  early_exit <- cbind(possibles, maxt)[which(maxt < tfinal & is.na(ratt1)), ]  # i,j pairs that exit early
  
  
  # Both early and late
  suppressWarnings( maxt1 <- sapply(1:nrow(late_entry), function(z) max(Ynew$t[Ynew$i == late_entry[z,1] & Ynew$j == late_entry[z,1]])) )
  suppressWarnings( ratt1 <- sapply(1:nrow(late_entry), function(z) Ynew$t[Ynew$ratification_year == 1 & Ynew$i == late_entry[z,1] & Ynew$j == late_entry[z,1]]) )
  maxt1[is.infinite(maxt1)] <- NA
  # ratt1[is.infinite(ratt1)] <- NA
  ratt2 <- rep(NA, length(ratt1))
  ratt2[which(sapply(ratt1, length) > 0)] <- unlist(ratt1[which(sapply(ratt1, length) > 0)])
  both_late_and_early <- cbind(late_entry, maxt1)[which(maxt1 < tfinal & is.na(ratt2)), ]  # i,j pairs that exit early
  if(nrow(both_late_and_early) > 0){
    both_late_and_early <- cbind( both_late_and_early[,c(1:2)], 0, both_late_and_early[,c(3:4)])
    entry_time <- sapply(1:nrow(both_late_and_early), function(z) late_entry[late_entry[,1] == both_late_and_early[z,1] & late_entry[,2] == both_late_and_early[z,2],3])
    both_late_and_early[,3] <- entry_time
    keep <- apply(both_late_and_early, 1, function(z) !(sum(is.na(z)) > 0))   # remove any NAs
    both_late_and_early <- both_late_and_early[keep, ]
    early_exit <- rbind(early_exit, both_late_and_early[, c(1,2,4)])  # add to early exit
  }
  
  # build possible signatures for all future times and simulations
  for(t in trange){
    k <- t - t0 + 1
    iposs <- cbind(possibles[rep(1:nrow(possibles), times=nsims), ], k, rep(1:nsims, each=nrow(possibles)))   # indices in first year that can be signed
    Yhat[,,k,] <- phat[,,k,] <- NA   # all NAs in time period
    Yhat[iposs] <- phat[iposs] <- 0   # possible ratifications
    
    jlate_entry <- which(t >= late_entry1[,3])   # AFTER late entry
    if(length(jlate_entry) > 0){
      # cat("late entry k=", k, "t=", t, "\n")
      ilate <- cbind(late_entry1[rep(jlate_entry, times=nsims), 1:2], k, rep(1:nsims, each=length(jlate_entry)) )
      Yhat[ilate] <- phat[ilate] <- 0   # 0s for i,j pairs that enter late
    }
    
    jearly_exit <- which(t > early_exit[,3])   # AFTER exit
    if(length(jearly_exit) > 0){
      # cat("early exit k=", k, "t=", t,"\n")
      iearly <- cbind(early_exit[rep(jearly_exit, times=nsims), 1:2], k, rep(1:nsims, each=length(jearly_exit)) )
      Yhat[iearly] <- phat[iearly] <- NA   # NAs for i,j pairs that leave early
    }
  }
  
  # Save all ratifications for future updates
  NAinit1 <- already_signed
  #unique(Yold[Yold$ratification_year == 1, c("i", "j")])   # save already signed treaty/country pairs
  if(!is.null(NApairs)){
    NAinit1 <- rbind(NAinit1, NApairs)
  }
  
  # # Set initial non-possible signatures, works wp1
  # possibles <- as.matrix(unique(Ynew[Ynew$t == t0, c("i", "j")]))   # all possible signatures in current year
  # for(j in 1:length(trange)){
  #   iposs <- cbind(possibles[rep(1:nrow(possibles), times=nsims), ], j, rep(1:nsims, each=nrow(possibles)))   # indices in first year that can be signed
  #   Yhat[,,j,] <- phat[,,j,] <- NA   # all NAs in first time period
  #   Yhat[iposs] <- phat[iposs] <- 0   # possible ratifications
  # }
  # 
  # # Now add back in treaties, countries that enter/exit
  # # Add NAs for countries/treaties that enter and leave the data set in question between t0 and tfinal
  # #   use suppressWarnings to deal with countries/treaties NOT in prediction year subset
  # suppressWarnings( crange <- t(sapply(1:S, function(z) range(Ynew$t[Ynew$i == z], finite=T))) )
  # suppressWarnings( treaty_range <- t(sapply(1:L, function(z) range(Ynew$t[Ynew$j == z], finite=T))) )
  # if(sum(is.infinite(crange))){ 
  #   crange[unique(which(is.infinite(crange), arr.ind = T)[,1]),] <- t0+1   # set infinite rows to be NA'd out in Yhat, phat
  # }
  # if(sum(is.infinite(treaty_range))){ 
  #   treaty_range[which(is.infinite(treaty_range), arr.ind = T)[,1],] <- t0+1  # set infinite rows to be NA'd out in Yhat, phat
  # }
  # 
  # ic <- which(crange[,1] > t0 | crange[,2] < tfinal)
  # for(i in ic){
  #   if(crange[i,1] > t0){
  #     i0 <- as.matrix( expand.grid(i, 1:L, crange[i,1]:tfinal - t0 +1, 1:nsims) )  # 0 indices
  #     iNA <- as.matrix( expand.grid(i, 1:L, t0:(crange[i,1]-1) - t0 +1, 1:nsims) )  # NA indices
  #   } else if (crange[i,2] < tfinal) {
  #     i0 <- as.matrix( expand.grid(i, 1:L, crange[i,2]:tfinal - t0 +1, 1:nsims) )  # 0 indices
  #     iNA <- as.matrix( expand.grid(i, 1:L, t0:(crange[i,2]-1) - t0 +1, 1:nsims) )  # NA indices
  #   }
  #   Yhat[i0] <- phat[i0] <- 0  
  #   Yhat[iNA] <- phat[iNA] <- NA 
  # }
  # 
  # jc <- which(treaty_range[,1] > t0 | treaty_range[,2] < tfinal)
  # for(j in jc){
  #   if(treaty_range[j,1] > t0){
  #     i0 <- as.matrix( expand.grid(1:S, j, treaty_range[j,1]:tfinal - t0 +1, 1:nsims) )  # 0 indices
  #     iNA <- as.matrix( expand.grid(1:S, j, t0:(treaty_range[j,1]-1) - t0 +1, 1:nsims) )  # NA indices
  #   } else if (treaty_range[j,2] < tfinal) {
  #     i0 <- as.matrix( expand.grid(1:S, j, treaty_range[j,2]:tfinal - t0 +1, 1:nsims) )  # 0 indices
  #     iNA <- as.matrix( expand.grid(1:S, j, t0:(treaty_range[j,2]-1) - t0 +1, 1:nsims) )  # NA indices
  #   }
  #   Yhat[i0] <- phat[i0] <- 0  
  #   Yhat[iNA] <- phat[iNA] <- NA 
  # }
  # 
  # 
  # # Save all ratifications for future updates
  # NAinit1 <- unique(Yold[Yold$ratification_year == 1, c("i", "j")])   # save already signed treaty/country pairs
  # if(!is.null(NApairs)){
  #   NAinit1 <- rbind(NAinit1, NApairs)
  # }
  
  
  # Prediction
  beta <- as.matrix(beta)
  allcoefs <- c(c(t(as.matrix(A))), c(as.matrix(B)), c(as.matrix(beta))[-1])   # coefs without intercept
  allcoefs[is.na(allcoefs)] <- 0
  keep <- which(allcoefs != 0)
  
  # Make dummy D and X variable arrays to avoid rebuilding
  #   shell is all pairs of i,j such that i<S and j<L for a single time period
  Dshell <- D[0,]
  Dshell[1:(S*L),] <- 0
  Xshell <- X[0,]
  Xshell[1:(S*L),] <- 0
  Dshell$t <- Xshell$t <- 1
  Dshell$i <- Xshell$i <- rep(1:S, times=L)
  Dshell$j <- Xshell$j <- rep(1:L, each=S)
  Dshell$cowcode <- Y$cowcode[match(Dshell$i, Y$i)]
  Dshell$treaty <- Y$cowcode[match(Dshell$j, Y$j)]
  Xshell$intercept <- 1
  Xshell[, -c(1:4)] <- NA   # NAs for all in Xshell
  
  # filename to save
  if(is.null(filename)){
    filename <-  paste0("predict_lag", lag, "_", min(year_range),"_", max(year_range), "_seed", seed0,  ".RData")
  }
  
  
  # Run loop to make predictions, write out results periodicially
  for(i in 1:nsims){
    set.seed(seed0+i-1)   # set seed for repeatability
    NAremove <- NAinit1   # initialize country/treaty pairs to remove from dataset for each simulation
    # Ytemp <- Y[Y$t < t0 & Y$ratification_year==1, c("i", "j", "t", "ratification_year")]   # save all ratifications
    
    for(t in t0:tfinal){    # t is the year of the prediction.
      k <- t - t0 + 1   # index in vavriables to save
      
      if(t == t0){    # initialize possibly lagged autoregressive array and covariate array
        Dpred <- Dnew[Dnew$t == t,] 
        Xpred <- Xnew[Xnew$t == t,]
      }   
      
      # Set impossible signatures to NA (should already be, but double-check)
      Yhat[cbind(NAremove, k, i)] <- phat[cbind(NAremove, k, i)] <- NA
      
      # Build design matrix for appropriate D, X
      Xpred$t <- Dpred$t <- 1
      Xreg <- build_design_additive(Dpred, Xpred, sparsedata=T, write=F, S=S, L=L, tmax=1)
      
      # Subset and pare design matrix to predict
      Xreg <- Xreg[,-(S^2 + L^2 + 1)]   # remove intercept
      keep_mat <- which(!is.na(Yhat[,,k,i]), arr.ind=T)   # i,j pairs to predict, in columnwise order 
      rows <- keep_mat[,1] + (keep_mat[,2]-1)*S   #  + (Y$t-1)*S*L   # unfolded indices
      Xreg <- Xreg[sort(rows), ]   # rows of Xreg that pertain to the entres in Yhat
      Xreg[which(is.na(Xreg))] <- 0   # Set NAs to zero
      
      # Calculate new prediction probabilities and predict
      keep_mat <- keep_mat[order(rows),]   # order entries as calculated
      keep_mat1 <- as.matrix(cbind(keep_mat, k, i))   # including time and sim indices
      # ptilde <- as.matrix(Xreg[, keep] %*% allcoefs[keep] + beta[1])
      ptilde <- as.matrix( as.matrix(Xreg[, keep] %*% allcoefs[keep]) + rep(beta[1], nrow(Xreg)) )
      phat_temp <- 1/(1 + exp( -c(ptilde) ))    # new probabilities
      phat[keep_mat1] <- phat_temp
      yhat_temp <- sapply(phat_temp, function(p) sample(c(0,1), 1, prob=c(1-p,p)))
      Yhat[keep_mat1] <- yhat_temp     # new matrix
      
      # Remove new ratifications from future entries 
      if(t < tfinal){
        NAnew <- which(Yhat[,,k,i] == 1, arr.ind=T)    # new signatures
        # temp <- as.matrix(cbind(NAnew, t, 1))
        # rownames(temp) <- 1:nrow(temp)
        # colnames(temp) <- c("i", "j", "t", "ratification_year")
        # Ytemp <- rbind(Ytemp, temp)   # save new signatures
        for(l in (k+1):dim(Yhat)[3]){
          NAentries <- as.matrix(cbind(NAnew, l, i))
          Yhat[NAentries] <- phat[NAentries] <- NA   # set all future i,j ratifications to NA
        }
        NAremove <- rbind(NAremove, NAnew)   
      }
      
      # Make new X array to predict from
      if(t < tfinal){
        nextyear <- unique(Y$year[Y$t == t+1])   # year for NEXT time period
        Xpred <- update_X_ijt(X, Yhat[,,,i], Yold, t+1, nextyear, region_income)  # covariates Xpred for NEXT time step
      }
      
      # Make new D array to predict from in NEXT time step
      if(t < tfinal){
        Dpred <- Dshell   # shell in which to save
        
        # save indices of which ratifications were signed in previous years
        signs <- NAnew   # signatures from last time period
        
        # If lag is greater than 1, need to augment Yhat with previous signatures
        if(lag > 1){
          if (t < t0+lag-1){   # if need to consult input Y for some signatures
            
            # split indices into those before prediction and after
            lagrange <- lag:1    
            keept <-  (t-lag+1):(t) < t0
            lold <- lagrange[keept]
            lhat <- lagrange[!keept]
            
            # augment with previous signatures from Y (i.e. not simulated)
            temp <- Y[Y$ratification_year==1 & Y$t %in% ((t-lag+1):(t))[keept], c("i", "j")]
            if(length(temp) > 0){
              signs <- rbind(signs, as.matrix(unique(temp)))  # old signatures
            }
            
            for(l in (lhat-1)){   # augment with previous signatures from Yhat, -1 accounts for fact that we are using this for next time period
              signs <- rbind(signs, as.matrix(which(Yhat[,,k-l,i] == 1, arr.ind=T)))   # updated signatures
            }
            
            
          } else {   # can use Yhat exclusively for signatures
            for(l in 1:(lag-1)){   # augment with previous signatures
              signs <- rbind(signs, as.matrix(which(Yhat[,,k-l,i] == 1, arr.ind=T)))
            }
          }
        }
        
        signs <- unique(signs)   # remove any duplicates
        Dpred[signs[,1] + (signs[,2] - 1)*S, response] <- 1   # save signatures regardless of lag
      } 
    }
    
    if(i %% write_interval == 0){
      save(Yhat, phat, S, L, countries, treaties, year_range, trange, response, seed0, NApairs, 
           file=file.path(outdir, filename))
      if(verbose){
        cat("done with sim", i, ";   ")
      }
    }
    
    
  }
  
  # Write out a final time
  save(Yhat, phat, S, L, countries, treaties, year_range, trange, response, seed0, NApairs, 
       file=file.path(outdir, filename))
  if(verbose){
    cat("\n*********************************\n")
    cat("DONE; saved to ", outdir, "\n")
  }
  
  return(list(Yhat=Yhat, phat=phat, S=S, L=L, countries=countries, treaties=treaties, year_range=year_range,trange=trange, response=response, seed0=seed0, NApairs=NApairs))
}





# Calculate prediction error of A and B after using time periods <= tmax
#   over the next tfuture time periods. 
# Return mean-square prediction error and mean absolute error
# type is biten, bilinear, or sadd
prediction_error_trade <- function(A, B, beta, tmax, tfwd, trade_file, states, lag, diff, diag_zero, type="biten", use_cov=T)
{
  
  if(tmax + tfwd > 20){stop("Too many time periods")}
  
  # Read in dataset
  data <- read_trade_data(trade_file, states, tmax=tmax + tfwd, lag, diff, diag_zero)   # read full dataset 
  D <- data$D
  X <- data$X
  Y <- data$Y
  country_names <- data$cnames
  
  S <- dim(Y)[1]
  L <- dim(Y)[2]   # square data
  
  tend <- (tmax - lag - diff + tfwd)   # last time index of data 
  Ynew <- predict_Y_forward(A,B,beta,X[,,1:tend,], D[,,1:tend], tpred=tfwd, type=type, use_cov=use_cov)   # predict over all prediction time periods
  
  sse_sq <- sapply(1:tfwd, function(x) sum((Ynew[,,x] - Y[,,(tmax - lag - diff + x)])^2))   # sse for squared loss all terms
  sse_sq_diag <- sapply(1:tfwd, function(x) sum(diag((Ynew[,,x] - Y[,,(tmax - lag - diff + x)])^2)))   # sse for squared loss diagonal terms
  
  sse_abs <-sapply(1:tfwd, function(x) sum(abs(Ynew[,,x] - Y[,,(tmax - lag - diff + x)])))   # sse for L1 loss all terms
  sse_abs_diag <-sapply(1:tfwd, function(x) sum(diag(abs(Ynew[,,x] - Y[,,(tmax - lag - diff + x)]))))   # sse for L1 loss diagonal terms
  
  mspe_diag <- sse_sq_diag / S
  mspe_offdiag <- (sse_sq - sse_sq_diag) / S / (L-1)
  made_diag <- sse_abs_diag / S
  made_offdiag <- (sse_abs - sse_abs_diag) / S / (L-1)
  
  
  return(list(mspe_diag=mspe_diag, mspe_offdiag = mspe_offdiag, made_diag=made_diag, made_offdiag=made_offdiag, Ynew=Ynew))
}


# Given A and B and data, predict the future Y
# return future Y slices
# first 3 dimensions of X and D must be the same
# tpred is the number of third slices to predict
#   uses the last tpred slices of X and D
# use_cov is flag to use covariates or not
# type is biten, bilinear, or sadd
predict_Y_forward <- function(A, B, beta, X, D, tpred, use_cov=T, type="biten")
{
  
  # Check if dimensions agree
  # if(sum(dim(Y)[c(1,2)] != dim(D)[c(1,2)]) > 0){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(D) != dim(X)[1:length(dim(D))]) > 0){  stop("Dimenions of D and X don't match")}
  if(length(dim(D)) < 3){ stop("Y is not an array")}
  
  
  t0 <- dim(D)[3] - tpred + 1
  t1 <- dim(D)[3]
  
  # if(t1 <1){stop("Need more time periods in D than in Y")}
  
  trange <- t0:t1   # range of time periods to predict upon
  
  if (strtrim(type,3) == "bil"){
    if(use_cov){   # calculate Xbeta based on whether to use covariates or not
      Xbeta <- drop(amprod(X[,,trange,], matrix(beta, nrow=1), length(dim(X[,,trange,]))))   # multiply Xbeta to start
    } else {
      Xbeta <-  0
    }
    Yout <- Xbeta + amprod(amprod(D[,,trange], A, 1), t(B), 2)    # bilinear multiliplicative model
    
  } else if (strtrim(type,3) == "sad" | strtrim(type,3) == "bit") {  # use regression prediction for additive models
    theta <- c(c(A), c(B))  # vector of covariates
    if(use_cov){ theta <- c(theta, c(beta)) }   # append beta coefficients if using them
    Xreg <- build_design_additive(D[,,trange], X[,,trange], type=type, use_cov=use_cov)  # build design matrix
    keep <- !is.na(theta)   # columns of X to keep
    
    Yout <- array(Xreg[,keep] %*% theta[keep], dim(D[,,trange]))
    
  } else { stop("Invalid model type for prediction")}
  
  return(Yout)
}


# Given A and B and data, predict Y
# first 3 dimensions of X and D must be the same
# uses all time dimensions of D
# use_cov is flag to use covariates or not
# type is biten, bilinear, or sadd
#
# 
predict_Y <- function(A, B, beta, X, D, use_cov=T, type="biten")
{
  
  # Check if dimensions agree
  # if(sum(dim(Y)[c(1,2)] != dim(D)[c(1,2)]) > 0){  stop("Dimenions of Y and D don't match")}
  if(sum(dim(D) != dim(X)[1:length(dim(D))]) > 0){  stop("Dimenions of D and X don't match")}
  if(length(dim(D)) < 3){ stop("Y is not an array")}
  
  
  if (strtrim(type,3) == "bil"){
    if(use_cov){   # calculate Xbeta based on whether to use covariates or not
      Xbeta <- drop(amprod(X, matrix(beta, nrow=1), length(dim(X))))   # multiply Xbeta to start
    } else {
      Xbeta <-  0
    }
    Yout <- Xbeta + drop(amprod(amprod(D, A, 1), t(B), 2))    # bilinear multiliplicative model
    
  } else if (strtrim(type,3) == "bit") {  # use regression prediction for additive models
    theta <- c(c(A), c(B))  # vector of covariates
    if(use_cov){ theta <- c(theta, c(beta)) }   # append beta coefficients if using them
    Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov)  # build design matrix
    keep <- !is.na(theta)   # columns of X to keep
    
    Yout <- array(Xreg[,keep] %*% theta[keep], dim(D))
    
  } else if (strtrim(type,3) == "sad") {
    if(use_cov){   # calculate Xbeta based on whether to use covariates or not
      Xbeta <- drop(amprod(X, matrix(beta, nrow=1), length(dim(X))))   # multiply Xbeta to start
    } else {
      Xbeta <-  0
    }
    S <- dim(D)[1]
    L <- dim(D)[2]
    Yout <- Xbeta + drop(amprod(amprod(D, A, 1), matrix(1,L,L), 2) + amprod(amprod(D, matrix(1,S,S), 1), t(B), 2) )    # bilinear multiliplicative model
    
  } else { stop("Invalid model type for prediction")}
  
  return(Yout)
}





# Predict Y using regression representation of additive models
predict_Y_regression <- function(A,B,beta,X,D,use_cov=T,type="sadd")
{
  theta <- c(c(A), c(B))  # vector of covariates
  if(use_cov){ theta <- c(theta, c(beta)) }   # append beta coefficients if using them
  Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov)  # build design matrix
  keep <- !is.na(theta)   # columns of X to keep
  
  Yout <- array(Xreg[,keep] %*% theta[keep], dim(D))
  return(Yout)
}


# Calculate the errors between two 3-mode arrays, along and off diagonal
# all is a flag to aggregate errors over all times together
# 
# Return mspe and made along and off-diagonal
array_errors <- function(Y, Ynew, all=F)
{
  if(sum(dim(Y) != dim(Ynew)) != 0){ stop("Dimensions of Y's do not match")}
  
  tfwd <- dim(Y)[3]  # number of time periods
  S <- dim(Y)[1]
  L <- dim(Y)[2]
  
  sse_sq <- sapply(1:tfwd, function(x) sum((Ynew[,,x] - Y[,,x])^2))   # sse for squared loss all terms
  sse_sq_diag <- sapply(1:tfwd, function(x) sum(diag((Ynew[,,x] - Y[,,x])^2)))   # sse for squared loss diagonal terms
  
  sse_abs <-sapply(1:tfwd, function(x) sum(abs(Ynew[,,x] - Y[,,x])))   # sse for L1 loss all terms
  sse_abs_diag <-sapply(1:tfwd, function(x) sum(diag(abs(Ynew[,,x] - Y[,,x]))))   # sse for L1 loss diagonal terms
  
  mspe_diag <- sse_sq_diag / S
  mspe_offdiag <- (sse_sq - sse_sq_diag) / S / (L-1)
  made_diag <- sse_abs_diag / S
  made_offdiag <- (sse_abs - sse_abs_diag) / S / (L-1)
  
  if(all){   # collapse over time if desired
    mspe_diag <- mean(mspe_diag)
    mspe_offdiag <- mean(mspe_offdiag)
    made_diag <- mean(made_diag)
    made_offdiag <- mean(made_offdiag)
  }
  
  return(list(mspe_diag=mspe_diag, mspe_offdiag = mspe_offdiag, made_diag=made_diag, made_offdiag=made_offdiag))
}


# Predict one step forward for environmental treaties data
# Y sparse array format
# X sparse array format (covariates)
#   predict for all countries/treaty pairs in A and B, and also return subset of only those in Y[last year]
# y0 is year to predict FROM
# must have X with time greater than this year, times should match between Y and X
Ypred_1step_envtreat_v2 <- function(Y, X, A, B, beta, y0, lag=1, use_cov=T)  #  countries, treaties, 
{
  invlogit <- function(x){ (1+exp(-x))^-1}
  
  # if(lag!=1){stop("Not coded for lags other than 1")}
  
  t0 <- Y$t[match(y0,Y$year)]
  if(!(max(X$t) > t0)){stop("Need X data after year y0")}
  
  # Renumber i,j and subset to countries/treaties only in last fit year(s), i.e. across lag
  # Y0 <- Y[Y$year == y0,]
  Y0 <- Y[Y$year %in% (y0 - lag + 1):y0,]
  countries <- as.numeric(rownames(A))
  treaties <- as.numeric(rownames(B))
  S <- length(countries)
  L <- length(treaties)
  Y0 <- Y0[Y0$cowcode %in% countries & Y0$treaty %in% treaties, ]
  Y0$i <- match(Y0$cowcode, countries)
  Y0$j <- match(Y0$treaty, treaties)
  
  # Build matrix to predict from
  Ya <- matrix(0, S, L)
  Ya[as.matrix(Y0[Y0$ratification_year == 1, c("i", "j")])] <- 1   # build matrix of current year to predict from
  
  # if(lag > 1){
  #   Yt <- Y[Y$year %in% (y0 - lag + 1):y0, ]
  #   Yt <- Yt[Y0$cowcode %in% countries & Y0$treaty %in% treaties, ]
  # } else { Yt <-  Y0 }
  # 
  # Ya <- matrix(0, S, L)
  # for(time in (y0 - lag + 1):y0){
  # Ya[as.matrix(Yt[Yt$ratification_year == 1, c("i", "j")])] <- 1   # build matrix of current year to predict from
  # }
  # 
  # Subset/order X and build X array
  if(use_cov){
    p <- length(beta)
    
    # Subset/order X
    ikeep <- Y$i[match(countries, Y$cowcode)]
    jkeep <- Y$j[match(treaties, Y$treaty)]
    timekeep <- t0+1
    Xkeep <- X[X$i %in% ikeep & X$j %in% jkeep & X$t %in% timekeep,]
    
    iold <- match(countries, sort(unique(Y$cowcode)))
    jold <- match(treaties, sort(unique(Y$treaty)))
    
    Xkeep$i <-  match(Xkeep$i, iold)   # new i
    Xkeep$j <-  match(Xkeep$j, jold)   # new i
    
    rows <- Xkeep$i + (Xkeep$j-1)*S + (Xkeep$t-1)*S*L   # unfolded indices
    Xkeep <- Xkeep[order(rows),]
    
    Xa <- array(0, c(S, L, p))
    colkeep <- (1:ncol(Xkeep))[!(names(Xkeep) %in% c("i", "j", "t"))]
    
    for(k in 1:p){
      Xa[ cbind(Xkeep$i, Xkeep$j, k)] <- Xkeep[ ,colkeep[k]]
    }
    
    Xa[is.na(Xa)] <- 0   # zero out NAs
    beta[is.na(beta)] <- 0  
    
    Xb <- drop(amprod(Xa, matrix(beta, nrow=1), 3))
    
  } else {Xb <- 0}
  
  # Renumber i,j indices 
  # Ytrue_a$i <- match(Ytrue_a$cowcode, ckeep)   # new i,j in country order
  # Ytrue_a$j <- match(Ytrue_a$treaty, tkeep)    # new i,j in treaty order
  # 
  # testindices <- as.matrix(Ytrue_a[Ytrue_a$t == min(Ytrue_a$t) & Ytrue_a$ratification_year == 0, c("i", "j")])
  # 
  # rows <- Ytrue_a$i + (Ytrue_a$j-1)*S + (Ytrue_a$t-1)*S*L   # unfolded indices
  # Ytrue_a <- Ytrue_a[order(rows),]
  
  
  
  
  # Predictions
  A[is.na(A)] <- 0
  B[is.na(B)] <- 0
  phat <- invlogit( t(A) %*% Ya + Ya %*% B + Xb )
  
  colnames(phat) <- treaties
  rownames(phat) <- countries
  
  Ypred <- Y0[,-4]
  Ypred$t <- Y0$t + 1
  Ypred$year <- Y0$year + 1
  Ypred$phat <- phat[cbind(as.character(Ypred$cowcode), as.character(Ypred$treaty))]
  
  
  return(list(phat=phat, Ypred=Ypred, countries=countries, treaties=treaties, ikeep=ikeep, jkeep=jkeep, pred_year=t0+1))
  
}


# Predict multiple steps forward for environmental treaties data, WITH ability to account for lag
# Y sparse array format
# X sparse array format (covariates)
#   predict for all countries/treaty pairs in A and B, and also return subset of only those in Y[last year]
# y0 is year to predict FROM
# must have X with time greater than this year, times should match between Y and X
Ypred_envtreat_lag <- function(Y, X, A, B, beta, y0, steps, sims=1e5, lag=1, use_cov=T, seed=NA, verbose=F)  #  countries, treaties, 
{
  invlogit <- function(x){ (1+exp(-x))^-1}
  rbern <- function(p,n=1){ sample(c(1,0), size=n, prob=c(p, 1-p), replace=T) }  # bernoulli r.v.
  
  # if(lag!=1){stop("Not coded for lags other than 1")}
  if(steps < 2){stop("Use Ypred_1step_envtreat_v2(.) for 1 step prediction")}
  if(max(Y$year) - steps <= y0){stop("Not enough years forward in Y for the number of steps")}
  
  t0 <- Y$t[match(y0,Y$year)]
  if(!(max(X$t) > t0)){stop("Need X data after year y0")}
  
  # Renumber i,j and subset to countries/treaties only in last fit year(s), i.e. across lag
  # Y0 <- Y[Y$year == y0,]
  Y0 <- Y[Y$year %in% (y0 - lag + 1):y0,]
  countries <- as.numeric(rownames(A))
  treaties <- as.numeric(rownames(B))
  S <- length(countries)
  L <- length(treaties)
  Y0 <- Y0[Y0$cowcode %in% countries & Y0$treaty %in% treaties, ]
  Y0$i <- match(Y0$cowcode, countries)
  Y0$j <- match(Y0$treaty, treaties)
  
  # Build lists of eligible country/treaty/year sets
  Yt <- Y[Y$cowcode %in% countries & Y$treaty %in% treaties & Y$year %in% (y0 + 1:steps), ]   # desired/eligible predictions
  Yt$i <- match(Yt$cowcode, countries)
  Yt$j <- match(Yt$treaty, treaties)   
  ijtkeep <- apply(cbind(Yt$cowcode, Yt$treaty, Yt$year), 1, paste0, collapse=",")
  
  ij <- cbind(  rep(rep(countries, times=L), times=steps), rep(rep(treaties, each=S), times=steps) )   # columnwise unfolding
  ijtall <- cbind(ij, rep(y0 + 1:steps, each=S*L))  # columnwise array unfolding
  ijtall <- apply(ijtall, 1, paste0, collapse=",")
  
  ijtremove <- setdiff(ijtall, ijtkeep)   # remove country/treaty sets
  Yremove0 <- Yremove00 <- t( sapply(ijtremove, function(z) as.numeric( unlist(strsplit(z, ","))) ) ) 
  Yremove00[,1] <- as.numeric( match(Yremove0[,1], countries) )   # in terms of numerical array indices
  Yremove00[,2] <- as.numeric( match(Yremove0[,2], treaties) )
  Yremove00[,3] <- as.numeric( match(Yremove0[,3], (y0 - lag + 1):(y0 + steps) ))  # 0:(steps + lag) shifts time index forward by lag
  
  # Build matrix to predict from
  Ya <- array(0, c(S, L, steps + lag))   # time indices are 1:steps + 1 and refer to y0 - lag + 1, y0 - lag + 2, ..., y0 + steps
  for(year in (y0 - lag + 1):y0){
    k <- match(year, (y0 - lag + 1):y0)  # array index
    Ya[ cbind( as.matrix(Y0[Y0$ratification_year == 1 & Y0$year == year, c("i", "j")]), k)  ] <- 1   # build matrix of current year to predict from
  }
  
  
  # Subset/order X and build X array
  if(use_cov){
    p <- length(beta)
    
    # Subset/order X
    ikeep <- Y$i[match(countries, Y$cowcode)]
    jkeep <- Y$j[match(treaties, Y$treaty)]
    timekeep <- t0 + 1:steps
    Xkeep <- X[X$i %in% ikeep & X$j %in% jkeep & X$t %in% timekeep,]
    
    iold <- match(countries, sort(unique(Y$cowcode)))
    jold <- match(treaties, sort(unique(Y$treaty)))
    
    Xkeep$i <-  match(Xkeep$i, iold)   # new i
    Xkeep$j <-  match(Xkeep$j, jold)   # new i
    
    rows <- Xkeep$i + (Xkeep$j-1)*S + (Xkeep$t-1)*S*L   # unfolded indices
    Xkeep <- Xkeep[order(rows),]
    
    Xa <- array(0, c(S, L, steps, p))
    colkeep <- (1:ncol(Xkeep))[!(names(Xkeep) %in% c("i", "j", "t"))]
    
    for(k in 1:p){
      Xa[ cbind(Xkeep$i, Xkeep$j, Xkeep$t - t0, k)] <- Xkeep[ ,colkeep[k]]
    }
    
    Xa[is.na(Xa)] <- 0   # zero out NAs
    beta[is.na(beta)] <- 0  
    
    Xb <- drop(amprod(Xa, matrix(beta, nrow=1), 4))
    
  } else {Xb <- 0}
  
  # Predictions
  A[is.na(A)] <- 0
  B[is.na(B)] <- 0
  
  phat <- array(0, c(S,L,steps))   # predicted values
  colnames(phat) <- colnames(Ya) <- treaties 
  rownames(phat) <- rownames(Ya) <- countries
  dimnames(phat)[[3]] <- y0 + 1:steps
  dimnames(Ya)[[3]] <- (y0 - lag + 1):(y0 + steps)
  Dtemp <- apply(Ya[,,1:lag, drop=F], 1:2, sum)   # lagged variable
  phat[,,1] <-  invlogit( t(A) %*% Dtemp + Dtemp %*% B + Xb[,,1] )  # first prediction straightforward!
  ptemp <- phat  # multiple-overwrite probability
  
  
  if(is.numeric(seed)){ set.seed(seed) }
  for(m in 1:sims){
    Yremove <- Yremove00  # reset zeros
    for(k in 2:steps){
      l <- k + lag - 1
      Ya[,,l] <- matrix(sapply(c(ptemp[,,k-1]), rbern), S, L)    # sample from previous time period
      Ya[Yremove] <- 0
      if(sum(Ya[,,l] > 0)){   # if treaties ratified, then add to remove list
        for(j in (l+1):(steps + lag)){
          Yremove <- rbind(Yremove, cbind(which(Ya[,,l] == 1, arr.ind=T), j) )   # no repeat 1's! all future times
        }
      }
      Dtemp <- apply(Ya[,,1:lag + k - 1, drop=F], 1:2, sum)   # lagged variable
      ptemp[,,k] <- invlogit( t(A) %*% Dtemp + Dtemp %*% B + Xb[,,k] )    # probability at next time period
      if(k > 2 & sum(Ya[,,l]) > 0){
        iadd <- which(Ya[,,l] == 1, arr.ind=T)
        phat[cbind(iadd, k-1)] <- phat[cbind(iadd, k-1)] + 1/sims   # add probabilities
      }
      if(k == steps){   # sample and save
        Ya[,,l+1] <- matrix(sapply(c(ptemp[,,k]), rbern), S, L)
        Ya[Yremove] <- 0
        if(sum(Ya[,,l+1]) > 0){
          iadd <- which(Ya[,,l+1] == 1, arr.ind=T)
          phat[cbind(iadd, k)] <- phat[cbind(iadd, k)] + 1/sims   # add probabilities
        }
      }
    }
    if(verbose & ((m %% 100) == 0) ){ cat("Prediction sim ", m, "done \n")}
  }
  
  # Build predicted dataset
  Ypred <- Yt
  Ypred <- Ypred[,names(Ypred) != "ratification_year"]
  Ypred$phat <- phat[ cbind(Ypred$i, Ypred$j, Ypred$t - t0)]
  
  
  return(list(phat=phat, Ypred=Ypred, countries=countries, treaties=treaties, ikeep=ikeep, jkeep=jkeep, pred_year=t0+1))
  
}



# update X covariates that depend on past ratifications, helper function for roll_forward_predict()
#   ijrats are i,j pairs that ratified at t
#   t is year for which X is updating to, and year is corresponding year
#   Y is an array from Yhat and Yold is a matrix of old data before t
update_X_ijt <- function(X, Y, Yold, t, year, region_income)   
{
  if(length(dim(Y)) != 3){ stop("Y must be a 3-mode array")}
  kold <- which(as.numeric(dimnames(Y)[[3]]) < year)
  
  Xnew <- X[X$t == t,]   # subset to t only to return
  Xold <- X[X$t == t-1,]   # previous step
  region_income_old <- region_income[region_income$year == year-1,]    # region and income of countries in t-1 year
  region_income_new <- region_income[region_income$year == year,]    # region and income of countries in t year
  
  Xnew$lagpercentincome <- Xnew$lagthreshold <- Xnew$lagpercentregion <- NA  # set all to NA to begin
  
  allregions <- unique(region_income$region)   # ALL regions of interest
  L <- dim(Y)[2]   # number of treaties
  S <- dim(Y)[1]   # number of countries
  allincomes <- 0:2    # unique incomes are same in every year
  
  treaties <- as.numeric(unique(dimnames(Y)[[2]]))   # unique treaties of interest, in order
  countries <- as.numeric(unique(dimnames(Y)[[1]]))   # unique countries of interest, in order
  cows2old <- lapply(1:L, function(z) unique(Yold$cowcode[Yold$j == z & Yold$t < t]))  # countries that have the opportunity sign each treaty in Yold 
  
  # compute lagthreshold based on ratifications
  ratsold <- sapply(1:L, function(z) sum(Yold$ratification_year[Yold$j == z & Yold$t < t], na.rm=T))  # ratifications for all years before t of each treaty
  ratsnew <- sapply(1:L, function(z) sum(Y[,z,kold], na.rm=T))  # ratifications for all years up to t-1
  ratsold[is.na(ratsold)] <- 0  # set NAs to zero
  ratsnew[is.na(ratsnew)] <- 0
  rats <- ratsold + ratsnew
  Xnew$lagthreshold <- rats[match(Xnew$j, 1:L)]    # updated ratifications, lagged by 1!
  
  cows1_byincome <- lapply(allincomes, function(z) unique(region_income_old$cowcode[region_income_old$income == z]))
  cows1_byregion <- lapply(allregions, function(z) unique(region_income_old$cowcode[region_income_old$region == z]))
  
  
  for(j in 1:L){
    cows2new_i <- unique( which(!is.na(Y[,j,kold, drop=F]), arr.ind=T)[,1] )   # i-country indices that had the opportunity to sign the given treaty
    cows2new <- countries[cows2new_i]   # countries that had the opportunity to ratify the given treaty in Y
    cows2 <- union(cows2old[[j]], cows2new)   # countries from any time < t that had the chance to sign that particular treaty
    cows2 <- union(cows2, c(265, 955))
    if(j == 14){   # account for treaties 40439 and 40434 signed in year 1950
      cows2 <- union(cows2, 390)
    }
    if(j == 16){   # account for treaties 40439 and 40434 signed in year 1950
      cows2 <- union(cows2, c(2, 94))
    }
    
    # Lag percent income updates
    for(i in 1:length(allincomes)){   # income code
      income <- as.numeric(allincomes[i])  # income in question
      cows1 <- cows1_byincome[[i]]
      
      cows <- intersect(cows1, cows2)   # countries in that particular income group in the previous year that had the chance to sign that particular the treaty
      cows <- cows[!is.na(cows)]   # remove any NAs
      cows_i <- (1:S)[match(cows, countries)]
      if(length(cows) > 0){  # save if there are any cows to update
        rowsold <- which(Yold$j == j & Yold$t < t & Yold$cowcode %in% cows)  # pertinent entries in Yold for previous years
        iold <- as.matrix( expand.grid(cows_i,j,kold) )   # pertinent entries in Yarray
        rowsnew <- which(Xnew$j == j & Xnew$t == t & Xnew$i %in% cows_i)  # pertinent rows for this year
        if(length(rowsnew) > 0){
          saveit <- length(rowsnew) > 0 | nrow(iold) > 0
          if(saveit){
            temp <- sum(Yold$ratification_year[rowsold], na.rm=T) + sum(Y[iold], na.rm=T)     # old and new ratifications
            Xnew$lagpercentincome[rowsnew] <- temp / length(cows)*100  # recalculate
          } else {
            Xnew$lagpercentincome[rowsnew] <- 0   # set to zero if no intersection
          }
        }
      }
    }
    
    
    # Lag percent region updates
    for(i in 1:length(allincomes)){   # income code
      cows1 <- cows1_byregion[[i]]
      
      cows <- intersect(cows1, cows2)   # countries in that particular income group in the previous year that had the chance to sign that particular the treaty
      cows <- cows[!is.na(cows)]   # remove any NAs
      cows_i <- (1:S)[match(cows, countries)]
      if(length(cows) > 0){  # save if there are any cows to update
        rowsold <- which(Yold$j == j & Yold$t < t & Yold$cowcode %in% cows)  # pertinent entries in Yold for previous years
        iold <- as.matrix( expand.grid(cows_i,j,kold) )   # pertinent entries in Yarray
        rowsnew <- which(Xnew$j == j & Xnew$t == t & Xnew$i %in% cows_i)  # pertinent rows for this year
        if(length(rowsnew) > 0){
          saveit <- length(rowsnew) > 0 | nrow(iold) > 0
          if(saveit){
            temp <- sum(Yold$ratification_year[rowsold], na.rm=T) + sum(Y[iold], na.rm=T)     # old and new ratifications
            Xnew$lagpercentregion[rowsnew] <- temp / length(cows)*100  # recalculate
          } else {
            Xnew$lagpercentregion[rowsnew] <- 0   # set to zero if no intersection
          }
        }
      }
    }
    
    
  }
  
  return(Xnew)
}




##########################
###  Update Functions  ###
##########################

# These functions perform updates for MLE block coordinate descents


# Update for block coordinate descent 
# asymmetric A and B
# D is matrix of lagged Y (covariates)
# Y is oservations
# FOR A SINGLE SLICE IN TIME!
# Allow for rank1 representation of U and Z for starting out 
update_MLE_asymmetric_t1 <- function(D, Y, U, V, W, Z)
{
  # sizes
  m <- ncol(W)
  p <- nrow(W)
  # k <- ncol(U)
  # n <- nrow(U)
  
  # Useful products
  BT <- tcrossprod(Z, W)   # Z %*% t(W)
  DD <- crossprod(D)    # t(D) %*% D
  
  
  # Simpler inverses when rank is 1, starting up
  # if(start){
  #   UUinv <- drop(1/crossprod(U[,1]))
  #   V <- solve(tcrossprod(D), (D %*% (crossprod(Y, U) - BT %*% crossprod(D, U)))) * UUinv
  #   
  #   # U update
  #   VD <- crossprod(V, D)   # t(V) %*% D)
  #   U <-  t(VD %*% (t(Y) - tcrossprod(BT, D))) / tcrossprod(VD)[1,1]
  #   
  #   # Z update
  #   DW <- D %*% W
  #   Z <- t( crossprod(DW, Y - U %*% crossprod(V, D))) / crossprod(DW, DW)[1,1]
  #   
  #   # W update
  #   W <- solve(DD, crossprod(D, Y %*% Z - U %*% VD %*% Z)) / crossprod(Z)[1,1]
  #   
  # } else {
  # 
  # V update
  V <- solve(tcrossprod(D), (D %*% (crossprod(Y, U) - BT %*% crossprod(D, U)))) %*% solve(crossprod(U))
  
  # U update
  VD <- crossprod(V, D)   # t(V) %*% D)
  U <- t(solve(tcrossprod(VD), VD %*% (t(Y) - tcrossprod(BT, D))))
  
  # Z update
  DW <- D %*% W
  Z <- t( solve(crossprod(DW, DW), crossprod(DW, Y - U %*% crossprod(V, D))))
  
  # W update
  W <- solve(DD, crossprod(D, Y %*% Z - U %*% VD %*% Z)) %*% solve(crossprod(Z))
  # }
  
  return(list(U=U, V=V, W=W, Z=Z))
}


# Update for block coordinate descent 
# asymmetric A and B
# D is matrix of lagged Y (covariates)
# Y is oservations
# for multiple time periods
update_MLE_asymmetric <- function(D, U, V, W, Z, DDT, DTD, DYT, DTY)
{
  # sizes
  m <- ncol(W)
  p <- nrow(W)
  k <- ncol(U)
  n <- nrow(U)
  
  t <- dim(D)[3]   # number of time slices
  
  # W and Z update
  DAD <- Reduce("+", lapply(1:t, function(x) crossprod(D[,,x ], tcrossprod(U, V) %*% D[,,x])  ))
  Z <- t( solve( crossprod(W, DTD %*% W), crossprod(W, DTY - DAD)))
  W <- solve(DTD, DTY - DAD) %*% Z %*% solve(crossprod(Z))
  
  # U and V update
  DBD <- Reduce("+", lapply(1:t, function(x) tcrossprod(D[,,x ] %*% tcrossprod(Z, W), D[,,x])  ))
  V <- solve(DDT, DYT - DBD) %*% U %*% solve(crossprod(U))
  U <- t(solve(crossprod(V, DDT %*% V), crossprod(V, DYT - DBD)))
  
  return(list(U=U, V=V, W=W, Z=Z))
}


# HOFF multiplicative bilinear model: Update for block coordinate descent 
# D is matrix of lagged Y (covariates)
# Y is oservations
# for multiple time periods
update_MLE_bilinear <- function(D, Y, A, B)
{
  
  t <- dim(D)[3]   # number of time slices
  
  # A update
  DBY <- Reduce("+", lapply(1:t, function(x) tcrossprod(D[,,x ] %*% B, Y[,,x] ) ))    # D B Y^T
  DBBD <- Reduce("+", lapply(1:t, function(x) tcrossprod( D[,,x] %*% B)  ))  # D B B^T D^T
  A <- t( solve( DBBD, DBY ))
  
  # B update
  DAY <- Reduce("+", lapply(1:t, function(x) crossprod(A %*% D[,,x ], Y[,,x])  ))    # D^T A^T Y
  DAAD <- Reduce("+", lapply(1:t, function(x) crossprod( A %*% D[,,x]  )))    # D^T A^T A D
  B <- ( solve( DAAD, DAY ))   # no transpose! 
  
  return(list(A=A, B=B))
}


# Update for block coordinate descent 
# asymmetric A and B
# D is matrix of lagged Y (covariates)
# Y is observations
# for multiple time periods
# type is sadd or biten additive model
update_MLE_additive <- function(A, B, D, DDT, DTD, DYT, DTY, type="biten")
{
  S <- dim(D)[1]
  L <- dim(D)[2]
  tmax <- dim(D)[3]
  
  if(type == "sadd"){
    Jl <- matrix(1,L,L)
    Js <- matrix(1,S,S)
    DBD <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, Js %*% D[,,x] %*% B)))    # D J B^T D^T J
    A <- t(solve(DDT, DYT - DBD))
    DAD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% A %*% D[,,x] %*% Jl)))    # D^T J A D J
    B <- (solve(DTD, DTY - DAD))   # no transpose! Use B instead of B^T
    
  } else if (tolower(strtrim(type,1)) == "b") {
    DBD <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], D[,,x] %*% B)))    # D B^T D^T
    A <- t(solve(DDT, DYT - DBD))
    DAD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], A %*% D[,,x])))    # D^T A D 
    B <- (solve(DTD, DTY - DAD))  # no transpose! Use B instead of B^T
    
  } else { stop( "Invalid type " ) }
  
  return(list(A=A, B=B))
}





##########################
###  Helper Functions  ###
##########################

# These functions read data, plot data, build useful matrices, and perform lower-level operations


# Read in trade data
# return list of X, Y, and D arrays with 0's along diagonal
# filename is location of trade .csv file
# tmax is maximum time period to consider, evaluate for 1:tmax
# states is subset of states to consider
#    states may be character vector or index vector
# lag is number of lags to consider
# diff is logical to difference responses or not
# 
read_trade_data <- function(filename, states=1:58, tmax=20, lag=3, diff=T, diag_zero=T, scale=F)
{
  
  trade.data <- read.table(filename, header=T)
  # head(trade.data)
  
  ################
  #create Y array#
  ################
  
  S <- max(trade.data$i) # number of exporting countries
  L <- max(trade.data$j) # number of importing countries
  T <- tmax  #max(trade.data$t) # number of time periods of trade
  
  # Y <- array(NA, c(S,L,tmax,1)) # create full data array
  # 
  # # fill in array structure with the trade data
  i_t <- which(trade.data$t <= tmax)
  # for(i in i_t){
  #   Y[trade.data$i[i], trade.data$j[i], trade.data$t[i], 1] <- trade.data$ltrade[i]
  # }
  # Y <- drop(Y)
  Y <- array(NA, c(S,L,tmax,1)) # create full data array
  Y[cbind(trade.data$i[i_t], trade.data$j[i_t], trade.data$t[i_t], 1)] <- trade.data$ltrade[i_t]
  
  dimnames(Y)[[1]] <- unique(levels(trade.data$exp))
  dimnames(Y)[[2]] <- unique(levels(trade.data$exp))
  
  
  
  ## code to subset out a portion of the data to analyze
  if(class(states) == "character"){
    sub <- sapply(states, function(x) which(dimnames(Y)[[1]] == x))
  } else if (class(states) == "integer"){
    sub <- states
  }
  # sub <- 1:58 #use all of the data for trade
  length(sub)
  dimnames(Y)[[1]][sub]
  
  S.sub <- length(sub) 
  L.sub <- length(sub)
  
  i_na <- which(is.na(Y), arr.ind=T)  # save out NA entries
  
  ################
  #create X array#
  ################
  p <- 7
  X <- array(NA, c(S,L,tmax,p))
  
  # X[,,,1] <- 1 # first covariate is an intercept
  # 
  # for(i in i_t){
  #   X[trade.data$i[i], trade.data$j[i], trade.data$t[i],2:p ] <- as.matrix(trade.data[i,8:13], nrow=1)
  # }
  
  X <- array(NA, c(S,L,tmax,p))
  X[,,,1] <- 1
  for(i in 2:p){
    X[cbind(trade.data$i[i_t], trade.data$j[i_t], trade.data$t[i_t], i)] <- trade.data[i_t,i + 6]
  }
  
  
  # Fix diagonals before processing -- assumes square Y, first two dimensions of X
  if(diag_zero){  # set diagonals to zero
    X[is.na(X)] <- 0
    for(i in 1:dim(X)[3]){
      diag(X[,,i,1]) <- 0
    }
    Y[is.na(Y)] <- 0
  } else {  # put averages of rows/columns in diagonal
    avgs <- sapply(1:dim(Y)[3], function(t) sapply(1:dim(Y)[1], function(z) mean(c(Y[z,-z,t,1], Y[-z,z,t,1]) ) ))  # averages of all diagonals
    Y[cbind( rep(1:dim(Y)[1], tmax), rep(1:dim(Y)[1], tmax), rep(1:tmax, each=dim(Y)[1]), 1 ) ] <- c(avgs)    # replace all at once!
    for(i in 1:p){
      avgs <- sapply(1:dim(X)[3], function(t) sapply(1:dim(X)[1], function(z) mean(c(X[z,-z,t,i], X[-z,z,t,i]) ) ))   # averages of diagonals for ith covariate
      X[cbind( rep(1:dim(X)[1], tmax), rep(1:dim(X)[1], tmax), rep(1:tmax, each=dim(X)[1]), i ) ] <- c(avgs)    # replace all averages at once
    }
  } 
  
  
  
  #create Y^t-1 and the difference that will go with A and B
  if (diff){
    Y_prev <- Y[,,(lag+1):(T-1),]
    Y_prevlag <- Y[,,1:(T-lag-1),] #we will have T-lag-1 timeperiods now
    if(scale) {
      normmat <- Y_prev
      normmat_lag <- Y_prevlag
    } else {
      normmat <- 1 #array(1, dim(Y_prev))
      normmat_lag <- 1 #array(1, dim(Y_prevlag))
    }
    Y <- Y[,,(lag+2):T,]
    X <- X[,,(lag+2):T,]
    D <- (Y - Y_prev)/normmat   #response
    D_lag <- (Y_prev - Y_prevlag)/normmat_lag   # previous changes which act as triggers of other events
    Z <- array(D_lag, c(S,L,T-lag-1,1))  # previous changes
    Y <- array(D, c(S,L,T-lag-1,1))      # outcome in model
    T <- T-lag-1
    
  } else if (!diff & !scale){
    D <- Y[,,(lag+1):T,]
    X <- X[,,(lag+1):T,]
    D_lag <-  Y[,,1:(T-lag),]
    Z <- array(D_lag, c(S,L,T-lag,1))  # previous changes
    Y <- array(D, c(S,L,T-lag,1))      # outcome in model    
    T <- T - lag
    
  } else if (!diff & scale) {
    D <- Y[,,(lag+2):T,] / Y[,,(lag+1):(T-1),]
    X <- X[,,(lag+2):T,]
    D_lag <-  Y[,,1:(T-lag-1),] /  Y[,,2:(T-lag),] 
    Z <- array(D_lag, c(S,L,T-lag-1,1))  # previous changes
    Y <- array(D, c(S,L,T-lag-1,1))      # outcome in model    
    T <- T - lag - 1
  }
  
  
  
  dimnames(Y)[[1]] <- unique(levels(trade.data$exp))
  dimnames(Y)[[2]] <- unique(levels(trade.data$exp))
  
  
  ###take subset of data/we used all data but could take subset here
  Y_sub <- Y[sub,sub,,]
  Y_sub <- array(Y_sub, c(S.sub,L.sub,T,1))
  Z_sub <- Z[sub,sub,,]
  Z_sub <- array(Z_sub, c(S.sub,L.sub,T,1))
  X_sub <- X[sub,sub,,]
  
  
  return(list(X=drop(X_sub), Y=drop(Y_sub), D=drop(Z_sub), cnames=unique(levels(trade.data$exp))[sub]))
}


# Read in forum data
# return list of Y and D arrays 
# filename is location of trade .csv file
# tmax is maximum time period to consider, evaluate for 1:tmax
# states is subset of states to consider
# lag is number of lags to consider
# diff is logical to difference responses or not
# 
read_forum_data <- function(filename, tmax=20, lag=3, diff=T, scale=F)  #, diag_zero=T, scale=F)
{
  
  trade.data <- read.table(filename, header=T)
  # head(trade.data)
  
  ################
  #create Y array#
  ################
  
  S <- max(trade.data$i) # number of exporting countries
  L <- max(trade.data$j) # number of importing countries
  T <- tmax  #max(trade.data$t) # number of time periods of trade
  
  i_t <- which(trade.data$t <= tmax)
  
  Y <- array(NA, c(S,L,tmax,1)) # create full data array
  Y[cbind(trade.data$i[i_t], trade.data$j[i_t], trade.data$t[i_t], 1)] <- trade.data$y[i_t]
  
  # dimnames(Y)[[1]] <- unique(levels(trade.data$exp))
  # dimnames(Y)[[2]] <- unique(levels(trade.data$exp))
  
  
  ## code to subset out a portion of the data to analyze
  # if(class(states) == "character"){
  #   sub <- sapply(states, function(x) which(dimnames(Y)[[1]] == x))
  # } else if (class(states) == "integer"){
  #   sub <- states
  # }
  # # sub <- 1:58 #use all of the data for trade
  # length(sub)
  # dimnames(Y)[[1]][sub]
  
  # S.sub <- length(sub) 
  # L.sub <- length(sub)
  
  
  
  #create Y^t-1 and the difference that will go with A and B
  if (diff){
    Y_prev <- Y[,,(lag+1):(T-1),]
    Y_prevlag <- Y[,,1:(T-lag-1),] #we will have T-lag-1 timeperiods now
    if(scale) {
      normmat <- Y_prev
      normmat_lag <- Y_prevlag
    } else {
      normmat <- 1 #array(1, dim(Y_prev))
      normmat_lag <- 1 #array(1, dim(Y_prevlag))
    }
    Y <- Y[,,(lag+2):T,]
    # X <- X[,,(lag+2):T,]
    D <- (Y - Y_prev)/normmat   #response
    D_lag <- (Y_prev - Y_prevlag)/normmat_lag   # previous changes which act as triggers of other events
    Z <- array(D_lag, c(S,L,T-lag-1,1))  # previous changes
    Y <- array(D, c(S,L,T-lag-1,1))      # outcome in model
    T <- T-lag-1
    
  } else if (!diff & !scale){
    D <- Y[,,(lag+1):T,]
    # X <- X[,,(lag+1):T,]
    D_lag <-  Y[,,1:(T-lag),]
    Z <- array(D_lag, c(S,L,T-lag,1))  # previous changes
    Y <- array(D, c(S,L,T-lag,1))      # outcome in model    
    T <- T - lag
    
  } else if (!diff & scale) {
    D <- Y[,,(lag+2):T,] / Y[,,(lag+1):(T-1),]
    # X <- X[,,(lag+2):T,]
    D_lag <-  Y[,,1:(T-lag-1),] /  Y[,,2:(T-lag),] 
    Z <- array(D_lag, c(S,L,T-lag-1,1))  # previous changes
    Y <- array(D, c(S,L,T-lag-1,1))      # outcome in model    
    T <- T - lag - 1
  }
  
  
  
  # dimnames(Y)[[1]] <- unique(levels(trade.data$exp))
  # dimnames(Y)[[2]] <- unique(levels(trade.data$exp))
  
  
  ###take subset of data/we used all data but could take subset here
  # Y_sub <- Y[sub,sub,,]
  # Y_sub <- array(Y_sub, c(S.sub,L.sub,T,1))
  # Z_sub <- Z[sub,sub,,]
  # Z_sub <- array(Z_sub, c(S.sub,L.sub,T,1))
  # X_sub <- X[sub,sub,,]
  
  
  # return(list(X=drop(X_sub), Y=drop(Y_sub), D=drop(Z_sub), cnames=unique(levels(trade.data$exp))[sub]))
  return(list(Y=drop(Y), D=drop(Z)))
}



# Build X for all ijt triples starting at t0 and ending at tfinal
#  contains rows for all ijt triples, not just thoes in the dataset
build_big_X <- function(t0, X, S=NULL, L=NULL, tfinal=NULL, readfile=F, wd=NULL)
{
  filename <- "Xbig.RData"
  if(readfile){
    if(!is.null(wd)){
      if(filename %in% list.files(wd)){
        load(filename)
        cat("Read-in big X \n")
        return(Xbig)
      } else {
        warning("Did not find Xbig.RData in wd... building manually")
      }
    } else {
      if(filename %in% list.files()){
        load(filename)
        cat("Read-in big X \n")
        return(Xbig)
      } else {
        warning("Did not find Xbig.RData in current working directory... building manually")
      }
    }
  }
  
  # Check data
  if(!all(c("i", "j", "t") %in% names(X))){stop("Need columns named i,j, and t")}
  if(!((t0-1) %in% unique(X$t))){stop("X must contain at least time period t-1")}
  
  # Build data size
  if(is.null(S)){S <- max(X$i)}
  if(is.null(L)){L <- max(X$j)}
  if(is.null(tfinal)){tfinal <- max(X$t)}
  Xsave <- X[X$t >= t0, ]   # only t greater than t0-1
  dt <- tfinal + 1 - t0   # number of time periods to save
  
  # Initialize new X
  Xnew <- X[0,]  # new X with same columns
  Xnew[1:(S*L*dt),] <- NA    # row size of new X
  Xnew[, c("i", "j", "t")] <- which(array(0, c(S,L,dt)) == 0, arr.ind=T) # ijt indices
  Xnew$t <- Xnew$t + t0 - 1   # increase t index to appropriate range, MUST DO THESE STEPS as relies on i,j,t presence
  Xnew$intercept <- 1   # intercept
  
  # Save existing entries in X
  rowsX <- Xsave$i + (Xsave$j-1)*S + (Xsave$t - t0)*S*L   # locations of X values in Xnew
  Xnew[rowsX,] <- Xsave   # save all old X values in larger data frame
  
  
  # entries in X that only vary by i,t
  irange <- c("polity", "lopen", "lnso2pc", "lgrgdpc", "memberships", "lnGDP")   
  # entries in X that only vary by j,t
  jrange <- c("hard_law2",  "global", "gdmix",  "ass_all", "ass_dev")   
  # entries in X that only vary by i,j,t
  ijtrange <- c("t_new", "t2_new", "t3_new")   # don't actually use these two, just notes for me
  ijonly <- c("leg_ap")
  
  # fill in i and j indices
  for(t in t0:tfinal){   # loop through time periods
    
    # i indices
    for(i in 1:S){    
      ii <- which(Xnew$i == i & Xnew$t == t)   # relevant rows
      itemp <- unique(Xnew[ii,irange])   # all rows should be the same except those that are all NA
      if(all(is.na(itemp))){
        itemp <- matrix(NA, 1, length(irange))   # if everything is NA, keep them that way
      } else {
        NArows <- apply(itemp, 1, function(z) all(is.na(z)))   # otherwise remove the NA row
        itemp <- itemp[!NArows,]
      }
      if(nrow(itemp) != 1){warning("There is some issue with filling in values that are constant among i,t")}   # if saving more than one row something went wrong
      
      NAii <- which(apply( Xnew[ii,irange], 1, function(z) all(is.na(z))))
      Xnew[ii[NAii],irange] <- itemp   # saving borrowed across all relevant pairs! only to those rows that are all NA
    }
    
    # j indices
    for(j in 1:L){    
      jj <- which(Xnew$j == j & Xnew$t == t)   # relevant rows
      jtemp <- unique(Xnew[jj,jrange])   # all rows should be the same except those that are all NA
      if(all(is.na(jtemp))){
        jtemp <- matrix(NA, 1, length(jrange))   # if everything is NA, keep them that way
      } else {
        NArows <- apply(jtemp, 1, function(z) all(is.na(z)))   # otherwise remove the NA row
        jtemp <- jtemp[!NArows,]
      }
      if(nrow(jtemp) != 1){warning("There is some issue with filling in values that are constant among j,t")}   # if saving more than one row something went wrong
      
      NAjj <- which(apply( Xnew[jj, jrange], 1, function(z) all(is.na(z))))
      Xnew[jj[NAjj], jrange] <- jtemp   # saving borrowed across all relevant pairs! only to those rows that are all NA    
    }
  }
  
  # fill in ij indices
  ijmin <- ijmax <- NULL
  for(i in 1:S){
    for(j in 1:L){
      ij <- which(Xnew$i == i & Xnew$j == j)
      NAij <- which(is.na(Xnew$t_new[ij]))
      ij0 <- which(Xnew$t_new[ij] == 0)
      
      # Work on zeros
      if(length(NAij) > 0){   # if there are possibly entries to replace
        if(length(ij0) > 0){  # if there is a zero 
          if(min(Xnew$t[ij[NAij]]) < Xnew$t[ij[ij0]]){   # if therre are NA entries before the zero
            ijreplace <- which( Xnew$t[ij[NAij]] < Xnew$t[ij[ij0]] )   # indices to replace
            Xnew$t_new[ij[ijreplace]] <- 0   # zeros before the zero
          }
        }
      }
      
      # Work on NAs beyond zeros
      NAij <- which(is.na(Xnew$t_new[ij]))
      if(length(NAij) > 0){   # if there are still entries to replace
        if(length(NAij) != length(ij)){   # if not all NAs to replace
          tm <- max(which(!is.na(Xnew$t_new[ij])))  # max index
          Xnew$t_new[ij[(tm+1):length(ij)]] <- Xnew$t_new[ij[tm]] + 1:length( (tm+1):length(ij) )
        } else {
          pasttemp <- X$t[X$i==i & X$j==j]  # go back to X for the ones greater then zero
          if(sum(is.na(pasttemp)) != length(pasttemp)){   # if not all past values are NAs
            t1 <- max(pasttemp, na.rm=T)   # reference time point
            maxt <- max(Xnew$t[ij], na.rm=T)
            ijreplace <- which(Xnew$t[ij] > maxt)   # replace those NAs for which t is larger
            Xnew$t_new[ij[ijreplace]] <- (X$t_new[X$i == i & X$j == j & X$t== t1] + Xnew$t[ij] - t1)[ijreplace]
          }
        }
      }
      
      
      # leg_ap
      latemp <- unique(X$leg_ap[X$i==i & X$j == j])  # possible legistlative approvals
      if(length(latemp) > 0){
        latemp1 <- latemp[!is.na(latemp)]  # non-NAs
        if(length(latemp1) > 1){
          warning("There exists at least one leg_ap variable that depends on t")
          cat("leg_ap: i",i,"j",j,"\n")
        }
        Xnew$leg_ap[ij] <- latemp1   # write regardless of whether there is more than 1
      }
      
      
    }
  }
  
  # Assign squared and cubed times
  Xnew$t2_new <- Xnew$t_new^2
  Xnew$t3_new <- Xnew$t_new^3
  
  
  return(Xnew)
}


# Function to load state interaction data
# allow for aggregation over variable weeks
read_state_interaction_data <- function(use_cov=T, agg = 1, lag=1, diff=F, remakeVAR=T)
{
  
  if(use_cov){
    load("YDX")
  } else { 
    load("YX")
    X<-X[,,,1,1,] 
    D <- X
    X <- array(NA, c(dim(D), 5))
    Xb <- 0
  }
  
  if(agg > 1 & is.numeric(agg)){
    tmax <- dim(Y)[4]
    aggmat <- matrix(1:tmax, nrow=agg)
    if(agg == 2){
      aggdims <- dim(aggmat)[2] - 1
    } else{
      aggdims <- floor(tmax/agg) + ((tmax%%agg) >= (agg/2)) # aggregate over last one only if at least half the data
    }
    Ynew <- Dnew <- array(NA, c(dim(Y)[1:3], aggdims-lag))
    Xnew <- array(NA, c(dim(X)[1:3], aggdims-lag, dim(X)[5]))
    aggmat[aggmat[,aggdims] < aggmat[1,aggdims], aggdims] <- 0   # zero out indices not to be included
    
    for(y in lag:aggdims){
      Ynew[,,,y-lag] <- apply(Y[,,,aggmat[,y]], 1:3, mean)   # aggregate
      Dnew[,,,y-lag] <- apply(Y[,,,aggmat[,(y-lag):(y-1)]], 1:3, mean) # next entry
      if(use_cov){
        Xnew[,,,y-lag,] <- apply(X[,,,aggmat[,y],], c(1:3,5), mean)   # average value
      }
    }
    
  } else if(agg == "year") {  # aggregate over years
    years <- as.numeric(unique(strtrim(dimnames(Y)[[4]], 4)))
    y_indices <- lapply(years, function(x) which(strtrim(dimnames(Y)[[4]], 4) == x))
    Ynew <- Dnew <- array(NA, c(dim(Y)[1:3], length(years)-1))
    Xnew <- array(NA, c(dim(X)[1:3], length(years)-1, dim(X)[5]))
    for(y in (lag + 1):length(y_indices)){
      Ynew[,,,y-lag] <- apply(Y[,,,y_indices[[y]]], 1:3, mean)   # aggregate
      Dnew[,,,y-lag] <- apply(Y[,,,unlist(y_indices[(y-lag):(y-1)])], 1:3, mean) # next entry
      if(use_cov){
        Xnew[,,,y-lag,] <- apply(X[,,,y_indices[[y]],], c(1:3,5), mean)   # average value
      }
    }
  } else if (agg == 1){
    Ynew <- Y[,,,lag:dim(Y)[4]]
    Xnew <- X[,,,lag:dim(X)[4],]
    Dnew <- D[,,,lag:dim(Y)[4]]
    if(remakeVAR | lag > 1){
      for(t in (lag+1):dim(Dnew)[4]){
        Dnew[,,,t] <- apply(Ynew[,,,(t-1):(t-lag)], 1:3, sum)   # sum previous lags
      }
      D[,,,1] <- apply(Ynew[,,,1:(lag-1)], 1:3, sum)
      Dnew[adiag(Dnew)] <- 0  # zeros on diagonals
    } 
    
  } else { stop("Something wrong with aggregation input")}
  
  if(diff){  # differencing of data
    Y <- Ynew
    D <- Dnew
    Xnew <- Xnew[,,,2:dim(Xnew)[4],]  # subset X
    Ynew <- Y[,,,2:dim(Y)[4]]
    Dnew <- D[,,,2:dim(D)[4]]
    for(t in 2:dim(Y)[4]){
      Ynew[,,,t-1] <- Y[,,,t] - Y[,,,t-1]  # difference Y, D
      Dnew[,,,t-1] <- D[,,,t] - D[,,,t-1]
    }
  }
  
  return(list(X=Xnew, Y=Ynew, D=Dnew, dnames=dimnames(X)))
}


# Function to load state interaction data
# allow for aggregation over variable weeks
read_state_interaction_data_multi <- function(use_cov=FALSE, agg=1, lagvec=c(1,1,1), diff=FALSE, remakeVAR=TRUE, readdir="~/Dropbox/BiTEN/gitBiTEN/hoff_code")
{
  lagvec <- as.numeric(lagvec)
  maxlag <- max(lagvec)
  ulag <- sort(unique(lagvec))
  agg <- as.numeric(agg)
  if(agg != 1){
    stop("only coding for agg = 1 right now")
  }
  
  setwd(readdir)
  Dlist <- vector("list", length(lagvec))
  for(i in ulag){
    data <- read_state_interaction_data(use_cov=use_cov, agg=agg, lag=i, diff=diff, remakeVAR=remakeVAR)
    j <- which(i == lagvec)
    for(k in j){
      Dlist[[k]] <- data$D
    }
    
    
    if(i == maxlag){
      Y <- data$Y
    }
  }
  
  # Trim D
  tmaxes <- sapply(Dlist, function(z) dim(z)[length(dim(z))])
  tmin <- min(tmaxes)
  t0 <- tmaxes - tmin + 1
  Dlist_out <- lapply(1:length(Dlist), function(z) Dlist[[z]][,,,t0[z]:tmaxes[z]])
  
  return(list(X=NA, Y=Y, D=Dlist_out))
}



# Read environmental treaties data from wd and write Y, D, X (sparse representation) to .RData file
read_env_treaties <- function(wd, lag, outdir=wd, write=F, readfile=T)
{
  # setwd(wd)
  
  # Check if processed already. If so, just load the data
  filename <- paste0("envtreat_lag",lag,".RData")
  if(filename %in% list.files() & readfile){   # check if already processed
    load(filename) 
    return(list(D=D, X=X, Y=y, yones=yones))
    
  } else {
    load(file.path(wd, "EventsCleanedFinal_TCY_02082018.RData"))
    # require(DataCombine)
    # require(car)
    
    countries <- sort(unique(dat_new$cowcode))
    treaties <- sort(unique(dat_new$treaty))
    
    #### Remove countries and treaties with zero ratificaitons
    signing <- unique(dat_new$cowcode[dat_new$ratification_year == 1])
    dat_new <- dat_new[dat_new$cowcode %in% signing, ]   # only countries that ratified treaties
    
    signing <- unique(dat_new$treaty[dat_new$ratification_year == 1])
    dat_new <- dat_new[dat_new$treaty %in% signing, ]   # only countries that ratified treaties
    
    countries <- sort(unique(dat_new$cowcode))
    treaties <- sort(unique(dat_new$treaty))
    ####
    
    
    #### Process covariates (code copied from Ben's file)
    # dat_orig <- dat_new
    # dat_new <- dat_orig
    # dat <- dat_new
    # Clean up the data
    
    # dat$polity <- with(dat, as.numeric(levels(polity))[polity]) 
    # dat$hard_law2 <- with(dat,  as.numeric(levels(hard_law2))[hard_law2])
    # dat$leg_ap <- with(dat, as.numeric(levels(leg_ap))[leg_ap])
    # dat$lnso2pc <- with(dat, as.numeric(levels(lnso2pc))[lnso2pc])
    # dat$lgrgdpc <- with(dat, as.numeric(levels(lgrgdpc))[lgrgdpc])
    # dat$ass_all <- with(dat, as.numeric(levels(ass_all))[ass_all])
    # dat$ass_dev <- with(dat, as.numeric(levels(ass_dev))[ass_dev])
    # dat$global <- with(dat, as.numeric(levels(global))[global])
    # dat$gdmix <- with(dat, as.numeric(levels(gdmix))[gdmix])
    # dat$lopen <- with(dat, as.numeric(levels(lopen))[lopen])
    # dat$lnGDP <- with(dat, as.numeric(levels(lnGDP))[lnGDP])
    # table(dat$regionpercent)
    # summary(dat$incomepercent)
    # 
    # # I need to lag
    # summary(dat$threshold_new)
    # summary(dat$incomepercent)
    # summary(dat$regionpercent)
    # # First, create a dyadID
    # 
    # id <- dat$cowcode*100000 + dat$treaty
    # dat$id <- id
    # 
    # library(DataCombine)
    # # newdata <- dat[order(dat$id, dat$year),]
    # # If no missingness, I can recode NAs as 0s for lags
    # table(is.na(dat$threshold_new))
    # dat_new <- slide(data = dat, 
    #                  Var = "threshold_new",
    #                  GroupVar = "id",
    #                  slideBy = -1,
    #                  NewVar = "lagthreshold",
    #                  keepInvalid = TRUE)
    # 
    # library(car)
    # dat_new$lagthreshold <- recode(dat_new$lagthreshold, "NA = 0")
    # table(dat_new$lagthreshold)
    # table(is.na(dat_new$lagthreshold))
    # dat <- dat_new
    # 
    # # income
    # dat_new <- slide(data = dat, 
    #                  Var = "incomepercent",
    #                  GroupVar = "id",
    #                  slideBy = -1,
    #                  NewVar = "lagpercentincome",
    #                  keepInvalid = TRUE)
    # table(is.na(dat_new$lagpercentincome))
    # dat_new$lagpercentincome <- recode(dat_new$lagpercentincome, "NA = 0")
    # dat <- dat_new
    # 
    # # region
    # dat_new <- slide(data = dat, 
    #                  Var = "regionpercent",
    #                  GroupVar = "id",
    #                  slideBy = -1,
    #                  NewVar = "lagpercentregion",
    #                  keepInvalid = TRUE)
    # table(is.na(dat_new$lagpercentregion))
    # dat_new$lagpercentregion <- recode(dat_new$lagpercentregion, "NA = 0")
    # dat <- dat_new
    
    # Add i,j,t indices
    X <- dat_new[,c("cowcode", "treaty", "year", "polity","hard_law2", "leg_ap", "lnso2pc", "lgrgdpc", "ass_all", "ass_dev", "global", "gdmix", "lopen", "memberships", "lagthreshold", "lagpercentregion", 
                "lagpercentincome", "lnGDP", "t_new", "t2_new", "t3_new")]   #  
    X <- cbind(1,X)
    names(X)[1] <- "intercept"
    
    X$t <- X$year - 1949
    X$i <- match(X$cowcode, countries)
    X$j <- match(X$treaty, treaties)
    X <- X[,c("i", "j", "t", "intercept", "polity","hard_law2", "leg_ap", "lnso2pc", "lgrgdpc", "ass_all", "ass_dev", "global", "gdmix", "lopen", "memberships", "lagthreshold", "lagpercentregion", 
           "lagpercentincome", "lnGDP", "t_new", "t2_new", "t3_new")]
      
    X <- X[X$t > lag, ]   # only times greater than lag
    X$t <- X$t - lag   # reset first time to 1
    ####
  
    
    #### Reorder X
    S <- max(X$i)
    L <- max(X$j)
    tmax <- max(X$t)
    
    rows <- X$i + (X$j-1)*S + (X$t-1)*S*L   # unfolded indices

    X <- X[order(rows),]
    ####
    
    
    #### y
    # dat_new <- dat
    
    # countries <- sort(unique(dat_new$cowcode))  # countries and treaties after subsetting
    # treaties <- sort(unique(dat_new$treaty))
    
    # ones <- which(dat_new$ratification_year == 1)
    # iones <- dat_new[ones, c("cowcode", "treaty", "year")]
    # iones$t <- iones$year - 1949
    # 
    # iones$i <- match(iones$cowcode, countries)
    # iones$j <- match(iones$treaty, treaties)
    # 
    # yones <- iones[, c("i", "j", "t", "cowcode", "treaty", "year")]   # reorder columns
    # yones <- yones[yones$t > lag,]  # subset for lag
    # yones$t <- yones$t - lag   # reset first time to 1
    # 
    # ikeep <- dat_new[, c("ratification_year", "cowcode", "treaty", "year")]
    # ikeep$t <- ikeep$year - 1949
    # ikeep$i <- match(ikeep$cowcode, countries)
    # ikeep$j <- match(ikeep$treaty, treaties)
    # y <- ikeep[, c("i", "j", "t", "ratification_year", "cowcode", "treaty", "year")]   # reorder columns
    
    y <- dat_new[, c("cowcode", "treaty", "year", "ratification_year")]
    y$t <- y$year - 1949
    y$i <- match(y$cowcode, countries)
    y$j <- match(y$treaty, treaties)

    S <- max(y$i)
    L <- max(y$j)
    tmax <- max(y$t)
    
    rows <- y$i + (y$j-1)*S + (y$t-1)*S*L   # unfolded indices
    
    y <- y[order(rows),]
    
    y <- y[,c("i", "j", "t", "ratification_year", "cowcode", "treaty", "year")]   # reorder columns
    
    
    #### Make D
    D <- y[y$t <= max(y$t) - lag,]   # lagged variables
    y <- y[y$t > lag,]  # subset for lag
    y$t <- y$t - lag   # reset first time to 1
    ####
    
    
    #### Post-process to get rid of unneeded data points (after treaty signed)
    #  actually doesn't find anything, but good check
    # check <- NULL
    # for(k in 1:nrow(yones)){
    #   remove <- which(y$i == yones$i[k] & y$j == yones$j[k] & y$t > yones$t[k])  # ratifications i,j 
    #   c(check, remove)
    #   if(length(remove) > 0){
    #     cat("removed ", remove)
    #     y <- y[-remove,]
    #     X <- X[-remove,]
    #     D <- D[-remove,]
    #   }
    # }
    ####
    
    if(write){
      save(D,X,y, file=file.path(wd, paste0("envtreat_lag",lag,".RData")))
    } 
    
    return(list(D=D, X=X, Y=y))
    
  }
  
}


# Add global/regional indicator covariate to X matrix
add_global <- function(X, treaties, readdir=getwd())
{
  # library("foreign")
  if(max(X$j) != length(treaties)){stop("Treaties not same size as j")}
  
  # tkeep <- read.dta(file.path(readdir, "Treaty_IDs.dta"))   # treaty IDs
  tkeep <- read.table(file.path(readdir, "Treaty_IDs.txt"), header=T)   # treaty IDs
  
  
  Xtreaties <- treaties[X$j]
  j <- match(Xtreaties, tkeep$treaty)
  X <- cbind(X, tkeep$global_treaty[j])
  names(X)[ncol(X)] <- "global_ind"
  
  return(X)
}

# Add EU indicator covariate to X matrix
add_eu_dummy <- function(data, readdir=getwd(), rand=F, seed=1)
{
  eus <- read.table(file.path(readdir, "eudummy.txt"), header=T)
  

  eudummy <- rep(0, nrow(data$Y))
  for(i in 1:nrow(eus)){
    ones <- data$Y$cowcode == eus$cow[i] & data$Y$year > eus$year[i]
    eudummy[ones] <- 1 
  }
  
  if(rand){   # if rand, scramble EU
    if(is.numeric(seed)){ set.seed(seed)}
    eudummy <- sample(eudummy)
  }
  
  data$X$eudummy <- eudummy
  
  return(data)
}


# Subset environmental treaties data
#   remove_y is the year(s) to remove from Y, i.e. removing 1951 from Y will remove 1950 from D (but not 1951!)
subset_env_treaties <- function(Y, D, X, cutoff_c, cutoff_t, remove_c=NULL, remove_t=NULL, remove_y=NULL, lag=1)
{
  if(sum(X[,c("i", "j", "t")] != Y[,c("i", "j", "t")]) != 0){
    stop("X and Y not in same order")
  }
  
  if( (!is.na(cutoff_c) | !is.na(cutoff_t)) | !is.null(remove_t) | !(is.null(remove_y)) | !(is.null(remove_c)) ){
    Y1 <- Y[Y$ratification_year==1,]
    cows <- sort(unique(Y$cowcode))
    remove_cowcode <- cows[sapply(cows, function(z) sum(Y$ratification_year[Y$cowcode == z]) <= cutoff_c)]
      #as.numeric(names(sort(table(Y1$cowcode)[table(Y1$cowcode) <= cutoff_c])))
    remove_cowcode <- c(remove_cowcode, remove_c)
    treats <- sort(unique(Y$treaty))
    remove_treaty <- treats[sapply(treats, function(z) sum(Y$ratification_year[Y$treaty == z]) <= cutoff_t)]
      #as.numeric(names(sort(table(Y1$treaty)[table(Y1$treaty) <= cutoff_t])))
    remove_treaty <- c(remove_treaty, remove_t)
    remove_year <- remove_y
    
    keep_cowcode <- sort(setdiff(unique(Y$cowcode), remove_cowcode))
    keep_treaty <- sort(setdiff(unique(Y$treaty), remove_treaty))
    keep_year <- sort(setdiff(unique(Y$year), remove_year))
    
    keepX <- (Y$cowcode %in% keep_cowcode) & (Y$treaty %in% keep_treaty) & (Y$year %in% keep_year)
    X <- X[keepX, ]
    Y <- Y[keepX,]
    keepD <- (D$cowcode %in% keep_cowcode) & (D$treaty %in% keep_treaty) & ((D$year + lag) %in% keep_year)
    D <- D[keepD,]
    # Old way of subsetting
    # Y <- Y[Y$cowcode %in% keep_cowcode, ]
    # Y <- Y[Y$treaty %in% keep_treaty, ]
    # Y <- Y[Y$year %in% keep_year, ]
    # D <- D[D$cowcode %in% keep_cowcode, ]
    # D <- D[D$treaty %in% keep_treaty, ]
    # D <- D[(D$year + lag) %in% keep_year, ]   # D is lagged

    Y$i <- X$i <- match(Y$cowcode, keep_cowcode)
    Y$j <- X$j <- match(Y$treaty, keep_treaty)
    Y$t <- X$t <- match(Y$year, keep_year)
    
    D$i <- match(D$cowcode, keep_cowcode)
    D$j <- match(D$treaty, keep_treaty)
    D$t <- match(D$year + lag, keep_year)
    
    countries <- keep_cowcode
    treaties <- keep_treaty
    years <- keep_year
    
  } else {
    countries <- sort(unique(Y$cowcode))
    treaties <- sort(unique(Y$treaty))
    years <- sort(unique(Y$year))
    
  } 
  return(list(Y=Y, D=D, X=X, countries=countries, treaties=treaties, yearsY=years))
}


# Subset environmental treaties data
#   lag is range of time dependence
#   i.e. year t depends on year(s) (t-lag):(t-1)
#   does not lag covariates
lag_env_treaties <- function(Y,D,X,lag)
{
  
  # Add up D's
  if(lag > 1){
    Dnew <- NULL
    # Dnew$ratification_year <- 0
    # Dpairs <- apply(Dnew[,c("cowcode", "treaty")], 1, function(z) paste0(z, collapse=","))
    for(i in min(Y$t):max(Y$t)){
      Yt <- Y[Y$t==i,]
      trange <- (i-lag):(i-1) + 1   # range of time periods previously from D
      Dt <- D[D$t %in% trange, ]
      # tpairs <- apply(Dt[,c("cowcode", "treaty")], 1, function(z) paste0(z, collapse=","))  # country/treaty pairs in final fit year
      # tpairs1 <- tpairs[Dt$ratification_year==1]
      Dt$t <- i  # renumber t 
      Dt$year <- max(Dt$year)  # renumber t
      
      # Reduce down into single time period (i.e. aggregate over lag)
      d0 <- unique(Dt[, c("i", "j", "t")])    # unique triplets in range of lags
      indices <- lapply(1:nrow(d0), function(z) which(Dt$i == d0[z,"i"] & Dt$j == d0[z,"j"] & Dt$t == d0[z,"t"]))   # rows of unique triplets
      Dt1 <- Dt[sapply(indices, function(z) z[1]),]   # only unique triplets
      Dt1$ratification_year <- sapply(indices, function(z) 1*(any(Dt$ratification_year[z] == 1)))   # if there is a ratification in any of the triplets, set to 1
      Dt <- Dt1
      
      Dnew <- rbind(Dnew, Dt)
    }
    D <- Dnew
    
    # Subset Y,D,X with lag
    Y <- Y[Y$year >= min(D$year) + lag, ]
    D <- D[D$year >= min(Y$year) - 1, ]
    X <- X[X$t %in% unique(Y$t), ]
    
    # Renumber t's
    if(any(sort(unique(Y$t)) != sort(unique(D$t)))){ stop("Something went wrong with indexing")}
    t0 <- min(D$t)
    D$t <- D$t - t0 + 1
    Y$t <- Y$t - t0 + 1
    X$t <- X$t - t0 + 1
  }
  
  return(list(Y=Y, D=D, X=X))
}


# Query treaties to keep/remove based on Tobias' recommendations
which_treaties <- function(treaty_type, readdir=getwd())   # added to MLE functions
{
  # library("foreign")
  if(!(strtrim(treaty_type, 3) %in% c("fra", "glo", "rob" ))){stop("Invalide treaty type")}
  
  # tkeep <- read.dta(file.path(readdir, "Treaty_IDs.dta"))
  tkeep <- read.table(file.path(readdir, "Treaty_IDs.txt"), header=T)   # treaty IDs
  j <- which(strtrim((names(tkeep)), 3) == strtrim(treaty_type, 3))
  
  keep <- tkeep$treaty[tkeep[,j] == 1]
  remove <- tkeep$treaty[tkeep[,j] == 0]
  
  return(list(keep=keep, remove=remove))
}




# Lag generic data
#   lag is range of time dependence
#   i.e. year t depends on year(s) (t-lag):(t-1)
#   does not lag covariates
lag_generic <- function(Y,D,X,lag, timename="year", response="ratification_year")
{
  
  # Add up D's
  if(lag > 1){
    Dnew <- NULL
    # Dnew$ratification_year <- 0
    # Dpairs <- apply(Dnew[,c("cowcode", "treaty")], 1, function(z) paste0(z, collapse=","))
    for(i in min(Y$t):max(Y$t)){
      Yt <- Y[Y$t==i,]
      trange <- (i-lag):(i-1) + 1   # range of time periods previously from D
      Dt <- D[D$t %in% trange, ]
      # tpairs <- apply(Dt[,c("cowcode", "treaty")], 1, function(z) paste0(z, collapse=","))  # country/treaty pairs in final fit year
      # tpairs1 <- tpairs[Dt$ratification_year==1]
      Dt$t <- i  # renumber t 
      Dt[, timename] <- max(Dt[, timename])  # renumber t
      
      # Reduce down into single time period (i.e. aggregate over lag)
      d0 <- unique(Dt[, c("i", "j", "t")])    # unique triplets in range of lags
      indices <- lapply(1:nrow(d0), function(z) which(Dt$i == d0[z,"i"] & Dt$j == d0[z,"j"] & Dt$t == d0[z,"t"]))   # rows of unique triplets
      Dt1 <- Dt[sapply(indices, function(z) z[1]),]   # only unique triplets
      Dt1[, response] <- sapply(indices, function(z) 1*(any(Dt[z, response] == 1)))   # if there is a ratification in any of the triplets, set to 1
      Dt <- Dt1
      
      Dnew <- rbind(Dnew, Dt)
    }
    D <- Dnew
    
    # Subset Y,D,X with lag
    Y <- Y[Y[, timename] >= min(D[, timename]) + lag, ]
    D <- D[D[, timename] >= min(Y[, timename]) - 1, ]
    X <- X[X$t %in% unique(Y$t), ]
    
    # Renumber t's
    if(any(sort(unique(Y$t)) != sort(unique(D$t)))){ stop("Something went wrong with indexing")}
    t0 <- min(D$t)
    D$t <- D$t - t0 + 1
    Y$t <- Y$t - t0 + 1
    X$t <- X$t - t0 + 1
  }
  
  return(list(Y=Y, D=D, X=X))
}




# Subset generic data
#   remove_y is the year(s) to remove from Y, i.e. removing 1951 from Y will remove 1950 from D (but not 1951!)
subset_generic <- function(Y, D, X, cutoff_c, cutoff_t, aname="cowcode", bname="treaty", timename="year", response="ratification_year",
                           remove_a=NULL, remove_b=NULL, remove_t=NULL, lag=1)
{
  if(sum(X[,c("i", "j", "t")] != Y[,c("i", "j", "t")]) != 0){
    stop("X and Y not in same order")
  }
  
  remove_y <- remove_t
  remove_t <- remove_b
  remove_c <- remove_a
  
  if( (!is.na(cutoff_c) | !is.na(cutoff_t)) | !is.null(remove_t) | !(is.null(remove_y)) | !(is.null(remove_c)) ){
    Y1 <- Y[Y[, response]==1,]
    cows <- sort(unique(Y[, aname]))
    remove_cowcode <- cows[sapply(cows, function(z) sum(Y[Y[, aname] == z, response]) <= cutoff_c)]
    #as.numeric(names(sort(table(Y1$cowcode)[table(Y1$cowcode) <= cutoff_c])))
    remove_cowcode <- c(remove_cowcode, remove_c)
    treats <- sort(unique(Y[, bname]))
    remove_treaty <- treats[sapply(treats, function(z) sum(Y[Y[, bname] == z, response]) <= cutoff_t)]
    #as.numeric(names(sort(table(Y1$treaty)[table(Y1$treaty) <= cutoff_t])))
    remove_treaty <- c(remove_treaty, remove_t)
    remove_year <- remove_y
    
    keep_cowcode <- sort(setdiff(unique(Y[, aname]), remove_cowcode))
    keep_treaty <- sort(setdiff(unique(Y[, bname]), remove_treaty))
    keep_year <- sort(setdiff(unique(Y[, timename]), remove_year))
    
    
    keepX <- (Y[, aname] %in% keep_cowcode) & (Y[, bname] %in% keep_treaty) & (Y[, timename] %in% keep_year)
    X <- X[keepX, ]
    Y <- Y[keepX,]
    keepD <- (D[, aname] %in% keep_cowcode) & (D[, bname] %in% keep_treaty) & ((D[, timename] + lag) %in% keep_year)
    D <- D[keepD,]
    # Old way of subsetting
    # Y <- Y[Y$cowcode %in% keep_cowcode, ]
    # Y <- Y[Y$treaty %in% keep_treaty, ]
    # Y <- Y[Y$year %in% keep_year, ]
    # D <- D[D$cowcode %in% keep_cowcode, ]
    # D <- D[D$treaty %in% keep_treaty, ]
    # D <- D[(D$year + lag) %in% keep_year, ]   # D is lagged
    
    Y$i <- X$i <- match(Y[, aname], keep_cowcode)
    Y$j <- X$j <- match(Y[, bname], keep_treaty)
    Y$t <- X$t <- match(Y[, timename], keep_year)
    
    D$i <- match(D[, aname], keep_cowcode)
    D$j <- match(D[, bname], keep_treaty)
    D$t <- match(D[, timename] + lag, keep_year)
    
    countries <- keep_cowcode
    treaties <- keep_treaty
    years <- keep_year
    
  } else {
    countries <- sort(unique(Y[, aname]))
    treaties <- sort(unique(Y[, bname]))
    years <- sort(unique(Y[, timename]))
    
  } 
  return(list(Y=Y, D=D, X=X, countries=countries, treaties=treaties, yearsY=years))
}


# Read in generic data
# dat is datafrome with only Y and X values
#  aname = "cowcode" or other names of A matrix
#  bname = "treaty" or other names of B matrix
#  timename = "year" or other time index
#
# Returns list of Y,D,X 
#
read_generic <- function(dat, aname, bname, timename, response)
{
  lag <- 1   # default lag
  # load(readfile)
  
  # if(!("dat" %in% ls())){stop("Data must be named dat in readfile")}
  
  dat
  
  aname <- as.character(aname)
  bname <- as.character(bname)
  timename <- as.character(timename)
  response <- as.character(response)
  
  dat <- data.frame(dat)
  
  #### remove unratified countries / treaties
  signing <- unique(dat[dat[, response] == 1, aname])
  dat <- dat[dat[, aname] %in% signing, ]   # only countries that ratified treaties
  
  signing <- unique(dat[dat[, response] == 1, bname])
  dat <- dat[dat[, bname] %in% signing, ]   # only countries that ratified treaties
  
  countries <- sort(unique(dat[, aname]))
  treaties <- sort(unique(dat[, bname]))
  mintime <- min(dat[, timename])
  ####
  
  
  #### X initial
  X <- dat
  X <- cbind(1,X)
  names(X)[1] <- "intercept"
  
  X$t <- X[, timename] - mintime + 1
  X$i <- match(X[, aname], countries)
  X$j <- match(X[, bname], treaties)
  ####
  
  
  #### Reorder X
  S <- max(X$i)
  L <- max(X$j)
  tmax <- max(X$t)
  
  rows <- X$i + (X$j-1)*S + (X$t-1)*S*L   # unfolded indices
  
  X <- X[order(rows),]
  ####
  
  
  
  #### Split into y and X
  y <- X[, c("i", "j", "t", response, aname, bname, timename)]
  X <- X[, c("i", "j", "t", setdiff(names(X), names(y)))]
  ####
  
  
  #### Make D
  D <- y[y$t <= max(y$t) - lag,]   # lagged variables
  y <- y[y$t > lag,]  # subset for lag
  y$t <- y$t - lag   # reset first time to 1
  X <- X[X$t > lag,]  # subset for lag
  X$t <- X$t - lag   # reset first time to 1
  ####
  
  
  return(list(D=D, X=X, Y=y))
}










# Function to build regression design matrix for additive models
# D is array to regress Y on for A and B, additive models only!
# X is array of covariates to regress Y upon for beta
# type is "biten" or "sadd"
# use_cov is boolean flag for using covariates, i.e. X/betas, or not
# 
#
# Returns design matrix
#    includes intercept if there is one in X (i.e. all covariates in X are included)
#    returns S*L*tmax \times S^2 + L^2 + ncol(X) matrix, with columnwise-vectorized order of rows
build_design_additive <- function(D, X, type="biten", use_cov=T, sparsedata=F, write=T, response="ratification_year", S=NULL, L=NULL, tmax=NULL)
{
  
  if(!sparsedata){
    # Check size
    if(length(dim(D)) != 3 ){ stop("D is not a 3-mode array") }
    if(sum(dim(D) != dim(X)[1:length(dim(D))]) > 0 & use_cov){  stop("Dimenions of D and X don't match")}
    
    # Find sizes
    if(is.null(S) & is.null(L) & is.null(tmax)){
      S = nrow(D[,,1])
      L = ncol(D[,,1])
      tmax = dim(D)[3]
    }

    
    # Build X matrix
    if(use_cov){   # start with beta columns base on use_cov flag
      p <- dim(X)[4]
      Xreg <- matrix(0, S*L*tmax, S^2 + L^2 + p)  # initialize X
      Xreg[, S^2 + L^2 + 1:p] <- t(mat(X, 4))  # beta columns
    } else {  
      p <-  0
      Xreg <- matrix(0, S*L*tmax, S^2 + L^2 + p)  # initialize X
    } 
    
    if(strtrim(type, 3) == "sad"){
      Js <- matrix(1, S, S)
      Jl <- matrix(1, L, L)
    } else if (strtrim(type, 3) == "bit"){
      Js <- diag(S)
      Jl <- diag(L)
    } else { stop("Invalid model type") }
    
    for(t in 1:tmax){  # vec A then vec B
      Xreg[ 1:(S*L) + (t - 1)*S*L, 1:S^2] <- kronecker(t(D[,,t] %*% Jl), diag(S))   # A columns
      Xreg[ 1:(S*L) + (t - 1)*S*L, S^2 + 1:L^2] <- kronecker(diag(L), Js %*% D[,,t])    # B columns
    }
  
  } else if (sparsedata){
    
    # Check if i,j,t in column names
    if( !("i" %in% names(D)) | !("j" %in% names(D)) | !("t" %in% names(D))){stop("D must have column names i,j, and t")}
    if( !("i" %in% names(X)) | !("j" %in% names(X)) | !("t" %in% names(X))){stop("X must have column names i,j, and t")}
    
    if(strtrim(type, 3) == "sad"){stop("SADD sparse not implemented")} 
    else if (strtrim(type, 3) == "bit"){
      
      filename <- paste0("Xreg_", type, "")
      
      # Find sizes
      if(is.null(S) & is.null(L) & is.null(tmax)){
        S <- max(D$i)
        L <- max(D$j)
        tmax <- max(D$t)
      }
      
      numcols <- S^2 + L^2 
      if(use_cov){
        p <- ncol(X) - 3   # remove i,j,t values
        numcols <- numcols + p
        
        X <- X[X$i <= S & X$j <= L & X$t <= tmax,]
      }
      
      # Xreg <- sparseMatrix(i=1,j=1,x=0, dims=c(tmax*S*L, numcols))   # initialize
      onerows <- which(as.vector(D[,response]) == 1)    # rows of D that have 1s in response
      reg1s <- matrix(0, length(onerows)*(S+L), 2)   # indices in Xreg that are 1s
      count <- 0
      for(k in onerows){
        count <- count+1
        # Xreg[cbind((D$j[k]-1)*S + 1:S + (D$t[k]-1)*S*L, (D$i[k]-1)*S + 1:S)] <- rep(1, S)   # A portion
        # Xreg[cbind((0:(L-1))*S + D$i[k] + (D$t[k]-1)*S*L, S^2 + (0:(L-1))*L + D$j[k])] <- rep(1, L)   # B portion
        i <- D$i[k]    ;   j <- D$j[k]   ;   t <- D$t[k]
        reg1s[1:S + (S+L)*(count-1),] <- cbind(S*L*(t-1) + S*(j -1) + 1:S, S*(i-1) + 1:S)
        reg1s[S + 1:L + (S+L)*(count-1),] <- cbind(S*L*(t-1) + i + S*(0:(L-1)), S^2 + j + L*(0:(L-1)))
      }
      Xreg <- sparseMatrix(i=reg1s[,1],j=reg1s[,2], x=1, dims=c(tmax*S*L, numcols))   # initialize
      if(use_cov){
        keep <- which(!(names(X) %in% c("i","j","t")))  # columns to keep that aren't i,j, or t
        Xreg[X$i + (X$j-1)*S + (X$t-1)*S*L, S^2 + L^2 + 1:p] <- as.matrix(X[, keep])
      }
      
    } else { stop("Invalid model type") }
    
  } else { stop("sparsedata must be true/false")}
  
  return(Xreg)
}




# Function to build regression design matrix for additive models
# D is array to regress Y on for A and B, additive models only!
# X is array of covariates to regress Y upon for beta
# type is "biten" or "sadd"
# use_cov is boolean flag for using covariates, i.e. X/betas, or not
# 
#
# Returns design matrix
#    includes intercept if there is one in X (i.e. all covariates in X are included)
#    returns S*L*tmax \times S^2 + L^2 + ncol(X) matrix, with columnwise-vectorized order of rows
build_design_additive_multi <- function(D,  S=NULL, L=NULL, tmax=NULL)
{
  D1 <- D[[1]]
  D2 <- D[[2]]
  D3 <- D[[3]]
  
  Ddims <- sapply(D, dim)
  if(!all(apply(Ddims, 1, vec.all.equal))){ stop("Not all lagged matrices in D have the same dimension")}
  Ddims <- dim(D[[1]])
  rm(D)
  
  X <- NA
  type <- "biten"
  use_cov <- FALSE 
  sparsedata=FALSE 
  write=FALSE
  response=NA 
  
    # Check size
    # if(length(dim(D)) != 4 ){ stop("D is not a 4-mode array") }
    # if(sum(dim(D) != dim(X)[1:length(dim(D))]) > 0 & use_cov){  stop("Dimenions of D and X don't match")}
    
    # Find sizes
    if(is.null(S) | is.null(L) | is.null(tmax)){
      S <- Ddims[1]
      L <- Ddims[2]
      R <- Ddims[3]
      tmax <- Ddims[4]
    }
    
    Xreg <- sparseMatrix(i=1, j=1, x=0, dims=c(S*L*R*tmax, S^2 + L^2 + R^2))  # initialize X
    
    
    if (strtrim(type, 3) == "bit"){
      Js <- diag(S)
      Jl <- diag(L)
      Jr <- diag(R)
    } else { stop("Invalid model type") }
    
    for(t in 1:tmax){  # vec A then vec B
      Xtemp <- NULL
      for(r in 1:R){
        temp <- cbind( kronecker(t(D1[,,r,t] %*% Jl), diag(S)), kronecker(diag(L), Js %*% D2[,,r,t]) )
        Xtemp <- rbind(Xtemp, temp)
      }
      Rcols <- sapply(1:R, function(z) c(D3[,,z,t]))
      Xtemp <- cbind(Xtemp, kronecker(Jr, Rcols))
      
      Xreg[ 1:(R*S*L) + (t - 1)*S*R*L, ] <- Xtemp
    }
    

  return(Xreg)
}




# Function to build regression design matrix for additive models
# D is array to regress Y on for A and B, additive models only!
# X is array of covariates to regress Y upon for beta
# type is "biten" or "sadd"
# use_cov is boolean flag for using covariates, i.e. X/betas, or not
# 
#
# Returns design matrix
#    includes intercept if there is one in X (i.e. all covariates in X are included)
#    returns S*L*tmax \times S^2 + L^2 + ncol(X) matrix, with columnwise-vectorized order of rows
build_design_additive_multi2 <- function(D, X=NA, verbose=FALSE)
{
  
  
  Ddims <- sapply(D, dim)
  if(!all(apply(Ddims, 1, vec.all.equal))){ stop("Not all lagged matrices in D have the same dimension")}
  Ddims <- dim(D[[1]])
  
  if(length(Ddims) > 4 | length(Ddims) < 3){ stop("Only coded for 2 or 3 modes of interest") }
  
  X <- NA
  type <- "biten"
  use_cov <- FALSE 
  sparsedata=FALSE 
  write=FALSE
  response=NA 
  
  
  # Find sizes
  tmax <- Ddims[length(Ddims)]
  ms <- Ddims[-length(Ddims)]
  
  S <- ms[1]
  L <- ms[2]
  if(length(Ddims) == 4){
    R <- ms[3]
    nentries <- R*tmax*S*L*(S+L) + (R*S*L)*tmax*R
    tripletX <- matrix(0, nentries, 3)
  }
  # Initialize X
  # Xreg <- sparseMatrix(i=1, j=1, x=0, dims=c(tmax*prod(ms), sum(ms^2)))   # initialize X
  # tripletX <- NULL
  
  
  if(length(Ddims) == 4){
    R <- ms[3]
    k <- prod(ms[1:2])
    rt_mat <- matrix(1:(R*tmax), R, tmax) 
    count <- 0
    for(t in 1:tmax){  # vec A then vec B and so on
      # Xtemp <- NULL
      for(r in 1:R){
        for(j in 1:max(S, L)){
          if(j <= S){
            rowsX <- (rt_mat[r,t]-1)*k + seq(1, k - S + 1, by=L) + j-1
            colsX <- seq(1, S^2 - S + 1, by=S) + j-1
            temp1 <- cbind(rep(rowsX, times=L), rep(colsX, each=S), c(t(D[[1]][,,r,t])))
            # tripletX <- rbind(tripletX, temp1)
            tripletX[count + 1:nrow(temp1), ] <- temp1
            count <- count + nrow(temp1)
          }
          if(j <= L){
            rowsX <- (rt_mat[r,t]-1)*k + (j-1)*S + 1:S
            colsX <- S^2 + (j-1)*L + 1:L
            temp2 <- cbind(rep(rowsX, times=S), rep(colsX, each=L), c(D[[2]][,,r,t]))
            # tripletX <- rbind(tripletX, temp2)
            tripletX[count + 1:nrow(temp2), ] <- temp2
            count <- count + nrow(temp2)
          }
          
        }
        # Xreg2[1:k + (rt_mat[r,t]-1)*k, ms[1]^2 + 1:(ms[2]^2)] <- kronecker( diag(Ddims[2]), D[[2]][,,r,t])
        # Xreg2[1:k + (rt_mat[r,t]-1)*k, 1:(ms[1]^2)] <- kronecker(t(D[[1]][,,r,t]), diag(Ddims[1]))
        # Xreg[1:k + (rt_mat[r,t]-1)*k, sum(ms[1:2]^2) + (r-1)*R + 1:R] <- matrix(D[[3]][,,,t], ncol=R)
        rowsX <- 1:k + (rt_mat[r,t]-1)*k
        colsX <- sum(ms[1:2]^2) + (r-1)*R + 1:R
        temp3 <- cbind(rep(rowsX, times=R), rep(colsX, each=k), c(matrix(D[[3]][,,,t], ncol=R)))
        # tripletX <- rbind(tripletX, temp3)
        tripletX[count + 1:nrow(temp3), ] <- temp3
        count <- count + nrow(temp3)
        
      }
      
      if(as.logical(verbose)){
        cat("t=", t, "out of ", tmax, "\n")
      }
      # Rcols <- sapply(1:Ddims[3], function(z) c(D[[3]][,,z,t]))
      # Xtemp <- cbind(Xtemp, kronecker( diag(Ddims[3]), Rcols))
      # Xreg[ 1:prod(Ddims[1:3]) + (t - 1)*prod(Ddims[1:3]), sum(ms[1:2]^2) + 1:(ms[3]^2)] <- kronecker( diag(Ddims[3]), Rcols)
    }
    
    Xreg <- sparseMatrix(i=tripletX[,1], j=tripletX[,2], x=tripletX[,3], dims=c(tmax*prod(ms), sum(ms^2)))   # initialize X
    
    
  } else if (length(Ddims) == 3) {
    for(t in 1:tmax){  # vec A then vec B and so on
      Xtemp <- NULL
      temp <- cbind( kronecker(t(D[[1]][,,t]), diag(Ddims[1]) ), kronecker( diag(Ddims[2]), D[[2]][,,t]) )
      Xtemp <- rbind(Xtemp, temp)
      
      # Rcols <- sapply(1:R, function(z) c(D3[,,z,t]))
      # Xtemp <- cbind(Xtemp, kronecker(Jr, Rcols))
      
      Xreg[ 1:prod(Ddims[1:2]) + (t - 1)*prod(Ddims[1:2]), ] <- Xtemp
    }
    
  }
  
  return(Xreg)
}




# Function to generate data for bilinear and additive models
# (S,L,tmax) is size of Y and D
# tau is error sd (iid normal, mean zero)
# muAB, sigmaAB are mean and sd for entries in U,V,W,Z s.th. A=UV^T and B=WZ^T
# rankA and rankB are the ranks for U,V and W,Z, respectively
# can make D different size with m,n, but not recommended
# gen_type is biten, bilinear, or sadd
# use_cov is boolean flag for using covariates, i.e. X/betas, or not
# genAB is boolean flag for generating AB or just using zeros
# if seed is numeric, set a seed
# binary is a flag for binary data
# 
# Returns Y, D, X, beta, A, B, and 2*LL (true values)
generate_data <- function(S, L, tmax, tau, 
                          muAB, sigmaAB, rankA, rankB, 
                          sigmaD, genAB=T, gen_type="biten", use_cov=T, seed=NA, 
                          binary=F, sparse=F)
{
  
  # Set seed if desired
  if(is.numeric(seed)){ set.seed(seed) }
  m <- S   # legacy sizes for different size D
  n <- L   
  
  # Generate D
  D <- array(rnorm(m*n*tmax, 0, sigmaD), c(m,n,tmax))
  if(binary){
    D <- 1*(D>0)  # threshold if binary
  }
  
  # Generate X
  if(use_cov){
    X1 <- array(1, c(S,L,tmax,1))
    X2 <- array(sample(c(0,1), S*L*tmax, replace=T), c(S,L,tmax,1))
    X3 <- array(rnorm(S*L*tmax), c(S,L,tmax,1))
    X <- abind(X1,X2,X3)
    p <- dim(X)[4]
    beta_true <- matrix(rep(1,p), nrow=1)
    Xbeta <- drop(amprod(X, beta_true, 4))
  } else {
    X <- Xbeta <- 0
    beta_true <- NA
  }
  
  # Generate A and B^T
  U_true <- matrix(rnorm(S*rankA, muAB, sigmaAB), S, rankA)
  V_true <- matrix(rnorm(S*rankA, muAB, sigmaAB), S, rankA)
  W_true <- matrix(rnorm(L*rankB, muAB, sigmaAB), L, rankB)
  Z_true <- matrix(rnorm(L*rankB, -muAB, sigmaAB), L, rankB)
  
  A_true <- tcrossprod(U_true, V_true)
  BT_true <- tcrossprod(Z_true, W_true)
  
  if(is.numeric(sparse)){   # set sparse % of elements to zero
    Aind <- matrix(sample(c(0,1), S^2, replace=T, prob=c(1-sparse, sparse)), S, S)
    Bind <- matrix(sample(c(0,1), L^2, replace=T, prob=c(1-sparse, sparse)), L, L)
    A_true <- Aind*A_true
    BT_true <- Bind*BT_true
  }
  
  # Generate Y's
  if(binary){
    # E <- array(rnorm(S*L*tmax, 0, tau*pi/sqrt(3)), c(S,L,tmax))
    E <- array(rlogis(S*L*tmax, 0, tau), c(S,L,tmax))
  } else {
    E <- array(rnorm(S*L*tmax, 0, tau), c(S,L,tmax))
  }
  
  if (strtrim(gen_type,3) == "bil"){
    Yout <- Xbeta + amprod(amprod(D, A_true, 1), BT_true, 2) + E   # bilinear multiliplicative model
    
  } else if (strtrim(gen_type,3) == "sad") {  
    Jl <- matrix(1, L, n)
    Js <- matrix(1, S, m)
    Yout <- Xbeta + amprod(amprod(D, A_true, 1), Jl, 2) + amprod(amprod(D, Js, 1), BT_true, 2) + E
    
  } else if (strtrim(gen_type,3) == "bit") { 
    A_true <- A_true*S^1.5/rankA
    BT_true <- BT_true*L^1.5/rankB   # increase variability for biten models
    Yout <- Xbeta + amprod(D, A_true, 1) + amprod(D, BT_true, 2) + E
    
  } else { stop("Invalid model type for prediction")}
  
  LLtrue <- -length(Yout)*log( sum( E^2) ) - length(Yout)   # Calc true 2*LL 
  
  if(binary){
    p <- c((1+exp(-Yout + E))^{-1})   # probabilities
    Yout <- 1*(Yout>0)  # threshold if binary
    LLtrue <- 2*sum(c(Yout)*log(p) + (1-c(Yout))*log(1-p))
  }
  
  return(list(Y=Yout, D=D, X=X, E=E, beta=beta_true, A=A_true, B=t(BT_true), LL2=LLtrue))
}


# Initialize estimates of U,V,W,Z
# does a decent job but doesn't appear to be quite enough yet
# for some reason fitting B first seems to work better. Bug???
#
# Estimates starting point for RankA and RankB matrices
initialize_mats <- function(D,Y,rankA=1,rankB=1, first="B")
{
  
  # if((rankA > 1) + (rankB > 1) > 0){
  #   stop("Not implemented for any rank > 1 yet")
  # }
  
  if(first == "A"){
    estA <- rank_1_estimate(tarray(D), tarray(Y))
    U <- estA$v  # since calulate A^T
    V <- estA$u  
    
    estB <- rank_1_estimate(D, Y -  amprod(D, tcrossprod(U, V),2))
    W <- estB$u
    Z <- estB$v
    
  } else if (first == "B"){
    estB <- rank_1_estimate(D, Y)
    W <- estB$u
    Z <- estB$v
    
    estA <- rank_1_estimate(tarray(D), tarray(Y - amprod(D, tcrossprod(Z,W), 2)))
    U <- estA$v  # since calulate A^T
    V <- estA$u 
  } else {
    stop("Choose A or B to fit first")
  }
  
  if(rankA > 1){
    splits <- c(0, round(1:rankA * length(U)/rankA ))
    Uout <- matrix(0, length(U), rankA)
    for (i in 1:rankA + 1){
      Uout[(splits[i-1]+ 1):splits[i], i-1] <- U[(splits[i-1]+ 1):splits[i]]
    }
    
    splits <- c(0, round(1:rankA * length(V)/rankA ))
    Vout <- matrix(0, length(V), rankA)
    for (i in 1:rankA + 1){
      Vout[(splits[i-1]+ 1):splits[i], i-1] <- V[(splits[i-1]+ 1):splits[i]]
    }
    U <- Uout
    V <- Vout
  }
  
  
  if(rankB > 1){
    splits <- c(0, round(1:rankB * length(W)/rankB ))
    Wout <- matrix(0, length(W), rankB)
    for (i in 1:rankB + 1){
      Wout[(splits[i-1]+ 1):splits[i], i-1] <- W[(splits[i-1]+ 1):splits[i]]
    }
    
    splits <- c(0, round(1:rankB * length(Z)/rankB ))
    Zout <- matrix(0, length(Z), rankB)
    for (i in 1:rankB + 1){
      Zout[(splits[i-1]+ 1):splits[i], i-1] <- Z[(splits[i-1]+ 1):splits[i]]
    }
    W <- Wout
    Z <- Zout
  }
  
  return(list(U=U,V=V,W=W,Z=Z))
}


# Initialize estimate of B assuming A = I. 
# Use regression formulation. 
# Take in Y and D.
#
# Returns B in LxL 
init_bilinear <- function(Y, D)
{
  # Check sizes of Y and D
  if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
  
  # Find sizes
  S <- nrow(D[,,1])
  L <- ncol(D[,,1])
  tmax <- dim(D)[3]
  
  # Full design matrix
  Xreg <- build_design_additive(D, X=0, type="biten", use_cov=F)
  Xb <- Xreg[, S^2 + 1:L^2]
  
  fit <- lm(c(Y) ~ -1 + Xb)
  B <- matrix(fit$coef, nrow=L)
  
  return(B)
}



# Matrix to estimate uv^T = B for Y \propto D %*% B, rank 1
# estimate u and v succesively via OLS
# Assume that dim(D) == dim(Y)
rank_1_estimate <- function(D, Y)
{
  
  if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
  
  if(length(dim(Y)) == 2){
    S <- nrow(D)
    L <- ncol(D)
    
    A_reg <- D[rep(1:S, times=L), ]
    fitA <- lm(c(Y) ~ -1 + A_reg)
    a <- fitA$coef
    A_reg2 <- matrix(0, S*L, L)
    A_reg2[cbind(1:(S*L), rep(1:L, each=S))] <- D %*% a 
    fitA2 <-  lm(c(Y) ~ -1 + A_reg2)
    a2 <- fitA2$coefficients
    
    u <- a
    v <- a2
    
  } else if (length(dim(Y)) == 3){
    S <- nrow(D[,,1])
    L <- ncol(D[,,1])
    tmax <- dim(D)[3]
    
    # unfold arrays to matrix case, probably should clean this up
    A_reg <- t(mat(D,2))[rep(1:S, times=L*tmax), ]
    A_reg2 <- matrix(0, S*L*tmax, L)
    # for(t in 1:tmax){
    #   A_reg[1:n^2 + (t-1)*n^2, ] <- D[rep(1:n, times=n),,t]
    # }
    fitA <- lm(c(Y) ~ -1 + A_reg)
    a <- fitA$coef
    for(t in 1:tmax){
      A_reg2[cbind(1:(S*L) + (t-1)*S*L, rep(1:L, each=S) )] <- D[,,t] %*% a
    }
    fitA2 <-  lm(c(Y) ~ -1 + A_reg2)
    a2 <- fitA2$coefficients
    
    u <- a
    v <- a2
  }
  
  
  return(list(u=u, v=v))
}


# array transpose
tarray <- function(A)
{
  if(length(dim(A)) == 2){
    Aout <- t(A)
  } else if(length(dim(A)) == 3) {
    Aout <- array(NA, dim(A)[c(2,1,3)])
    for(i in 1:dim(A)[3]){
      Aout[,,i] <- t(A[,,i])
    }
  }
  
  return(Aout)
}


# Make arcplot of networks
arcplot_fm2 <- function(mplot, groups, colors, linespace=.05, maxwt=NULL, maxdeg=NULL, horizontal=T, cex.lab=1.0)   # positive and negative together
{
  mplot[is.na(mplot)] <- 0
  if(any(rownames(mplot) != colnames(mplot))){stop("Row names and column names don't match")}
  # if(any(c(mplot) < 0)){stop("All entries in mplot must be positive")}
  if(length(colors) < length(unique(groups))){"need at least as many colors as groups"}
  
  
  reorder <- order(groups)
  mplot <- mplot[reorder, reorder]
  groups <- groups[reorder]
  # linespace <- .05   # additional label spacing
  
  nodecolors <- colors[match(groups, unique(groups))]
  
  edgelist <- cbind(rep(colnames(mplot), times=length(rownames(mplot))),   
                    rep(rownames(mplot), each=length(colnames(mplot))))  # columnwise vectorization
  weights <- c(mplot)
  weights[is.na(weights)] <- 0
  weights[weights == 0] <- 0.0
  if(is.null(maxwt)){ maxwt <- max(abs(weights))}
  weights <- weights/maxwt
  
  above <- weights > 0  
  
  nodesizes <- rowSums(mplot!=0, na.rm=T) + colSums(mplot!=0, na.rm=T)
  if(is.null(maxdeg)){ maxdeg <- max(abs(nodesizes))}
  nodesizes <- nodesizes/maxdeg
  
  colorstrans <- paste0(colors, "65")
  arccolors <- colorstrans[  match( groups[ match(edgelist[,1], rownames(mplot)) ],  unique(groups) )  ]  # group color of 1st edge
  arccolors[weights == 0] <- NA
  
  arcplot(edgelist,  cex.labels = cex.lab,  # ordering = new_ord,  #labels = vlabels,
          lwd.arcs = 14*abs(weights), 
          lwd.nodes=1, show.nodes = TRUE, col.nodes = NA, bg.nodes = nodecolors, labels=rownames(mplot),
          cex.nodes = 2.4*nodesizes + .6, pch.nodes = 21,  line = 0.1 +linespace,
          col.arcs = arccolors, col.labels="gray50", above=above, 
          horizontal=horizontal)
}




# Make arcplot of networks
arcplot_fm3 <- function(mplot, groupmat, colors=NULL, linespace=.05, maxwt=NULL, maxdeg=NULL, 
                        horizontal=T, cex.lab=1.0, knode=1, allnodesizes=NULL)   # positive and negative together
{
  mplot[is.na(mplot)] <- 0
  if(any(rownames(mplot) != colnames(mplot))){stop("Row names and column names don't match")}
  # if(any(c(mplot) < 0)){stop("All entries in mplot must be positive")}
  # if(length(colors) < length(unique(groups))){"need at least as many colors as groups"}
  
  
  
  # Ordering correction
  groups <- groupmat[match(rownames(mplot), groupmat[,1]), 2]
  reorder <- order(groups)
  mplot <- mplot[reorder, reorder]
  groups <- groups[reorder]
  
  # Remove zeros
  remove <- which( rowSums(mplot != 0) + colSums(mplot != 0) == 0)
  keep <- setdiff(1:nrow(mplot), remove)
  mplot <- mplot[keep, keep]
  groups <- groups[keep]
  
  
  if(is.null(colors)){
    L <- length(unique(groups))
    colors <- brewer.pal(2*L - 1, "Spectral")[seq(1, 2*L - 1, by = 2)]
  } else {
    colors <- colors[match(unique(groups), sort(unique(groupmat[,2])))]
  }
  
  
  # Actuual edgelist
  i1 <- which(mplot != 0, arr.ind=TRUE)
  edgelist <- cbind( rownames(mplot)[i1[,1]], colnames(mplot)[i1[,2]] )
  
  
  nodes <- unique(c(t(edgelist)))   # internal arcplot ordering
  nodegroups = groupmat[match(nodes, groupmat[,1]),2]
  # groups = groupmat[,2]
  nodecolors <- colors[match(nodegroups, unique(groups))]
  

  weights <- mplot[i1]
  weights[is.na(weights)] <- 0
  weights[weights == 0] <- 0.0
  if(is.null(maxwt)){ maxwt <- max(abs(weights))}
  weights <- weights/maxwt
  
  above <- weights > 0  
  
  # nodesizes <- (rowSums(mplot!=0, na.rm=T) + colSums(mplot!=0, na.rm=T))  
  if(is.null(allnodesizes)){
    nodesizes <- rowSums(abs(mplot))  #+ colSums(abs(mplot)))  
  } else {
    nodesizes <- allnodesizes[match(nodes, colnames(allnodesizes))]
  }

  if(is.null(maxdeg)){ maxdeg <- max(abs(nodesizes))}
  nodesizes <- nodesizes/maxdeg
  
  colorstrans <- paste0(colors, "65")
  arccolors <- colorstrans[  match( groups[ match(edgelist[,1], rownames(mplot)) ],  unique(groups) )  ]  # group color of 1st edge
  # arccolors[weights == 0] <- NAC
  
  arcplot(edgelist,  cex.labels = cex.lab,  # ordering = new_ord,  #labels = vlabels,
          lwd.arcs = 14*abs(weights), 
          lwd.nodes=1, show.nodes = TRUE, col.nodes = NA, bg.nodes = nodecolors, ordering=rownames(mplot),
          cex.nodes = 2.4*nodesizes*knode + .6, pch.nodes = 21,  line = 0.1 +linespace,
          col.arcs = arccolors, col.labels="gray50", above=above, 
          horizontal=horizontal)
}


# Plot the matrix as previously, with colored regions
plot_matrix <- function(A, plottitle, filename, outdir, country_names, write=F)
{
  
  
  colors.all <- data.frame(country=country_names)
  colors.all$region <- countrycode(colors.all$country, "country.name", "region")
  colors.all$continent <- countrycode(colors.all$country, "country.name", "continent")
  
  colors.all$region <- factor(colors.all$region, levels=c("Eastern Africa", "Northern Africa",
                                                          "Caribbean", "Central America", "Northern America", "South America",
                                                          "Eastern Asia", "Southern Asia", "South-Eastern Asia", "Western Asia",
                                                          "Northern Europe", "Southern Europe", "Western Europe",
                                                          "Australia and New Zealand"))    
  
  # cbind(as.character(colors.all$region), as.numeric(colors.all$region))
  
  color.choice.region <- c("pink", "magenta", "slateblue3", "skyblue2",
                           "navy", "royalblue", "red", "coral", "firebrick4",
                           "tomato2", "green", "lightgreen", "darkgreen",
                           "yellow3")
  
  color.choice.cont <- c("pink", "blue", "red", "green", "yellow3")
  
  colors.all$region.col <- color.choice.region[as.numeric(colors.all$region)]
  colors.all$cont.col <- color.choice.cont[as.factor(colors.all$continent)]
  
  # pdf("OUTPUT/ColorMapsAll.pdf")
  
  rownames(A) <- country_names
  colnames(A) <- country_names
  
  flip <- function(a){ t(a)[,nrow(a):1]}  #for image.plot to be in same order of matrix
  #http://blog.snap.uaf.edu/2012/06/08/matrix-rotation-for-image-and-contour-plots-in-r/
  
  ###### by continent
  
  #A
  #sort by continent index
  ind.cont <- sort(colors.all$continent, index.return=T)$ix
  
  #for coloring axis labels
  #http://stackoverflow.com/questions/18839731/vary-colors-of-axis-labels-in-r-based-on-another-variable
  
  par(mar=c(6,6,3,2), cex.axis=.5)
  if(write==T){
    pdf(file.path(outdir, filename), width=7, height=6)
  }
  par(mar=c(6,6,3,2), cex.axis=.5)
  image.plot(flip(A[ind.cont,ind.cont]), axes=F)
  title(plottitle)
  S.sub <- L.sub <- dim(Y)[1]
  loc.S <- seq(0,1,length=S.sub)
  loc.L <- seq(0,1,length=L.sub)
  for(j in 1:S.sub){
    axis(2,labels=colnames(flip(A[ind.cont,ind.cont]))[j], at=loc.S[j], las=1,
         col.axis=rev(colors.all$cont.col[ind.cont])[j])
  }
  for(j in 1:L.sub){
    axis(1,labels=rownames(flip(A[ind.cont,ind.cont]))[j], at=loc.L[j], las=3,
         col.axis=(colors.all$cont.col[ind.cont])[j])
  }
  if(write==T){
    dev.off()
  }
  
  
}

# Bare bones linear model solver, coefficients may not be unique
solve_lm <- function(x,y)
{
  remove <- which(is.na(y))  # remove NAs
  
  if(length(remove) > 0){
    coef <- ginv(crossprod(x[-remove, ])) %*% crossprod(x[-remove, ],y[-remove])
  } else {
    coef <- ginv(crossprod(x)) %*% crossprod(x,y)
  }
  
  pred <- x %*% coef
  return(list(coef=coef, pred=pred))
}

# Check for duplicate columns in matrix M. Return matrix of pairs of columns. Works with probability 1. 
checkcols <- function(M) 
{   # return entries with same columns
  b = rnorm(nrow(M))
  a = as.vector(t(M)%*%b)
  diffs <- c()
  for(i in 1:(length(a) - 1)){   # for-loop for memory size
    for(j in (i+1):length(a)){
      if(a[i] - a[j] == 0){ diffs <- rbind(diffs, c(i,j)) }
    }
  }
  return(diffs)
}

# plot a scatter of the vectorization of y_hat against y
plot_scatter <- function(y_hat, y, plottitle, filename, outdir, write=F)
{
  
  if(write==T){
    pdf(file.path(outdir, filename), width=6, height=6)
  }
  par(mar=c(4,4,4,.5), cex.axis=.5)
  plot(c(y),c(y_hat),xlab="vec(Y)", ylab=expression("vec(hat(Y))"), main = plottitle)
  abline(0,1,col="red")
  # S.sub <- L.sub <- dim(Y)[1]
  # loc.S <- seq(0,1,length=S.sub)
  # loc.L <- seq(0,1,length=L.sub)
  # for(j in 1:S.sub){
  #   axis(2,labels=colnames(flip(A[ind.cont,ind.cont]))[j], at=loc.S[j], las=1,
  #        col.axis=rev(colors.all$cont.col[ind.cont])[j])
  # }
  # for(j in 1:L.sub){
  #   axis(1,labels=rownames(flip(A[ind.cont,ind.cont]))[j], at=loc.L[j], las=3,
  #        col.axis=(colors.all$cont.col[ind.cont])[j])
  # }
  if(write==T){
    dev.off()
  }
  
  
}

# find diagonal indices of k-mode array, where "diagonal" is i=j for first two indices of input array Y
adiag <- function(Y)
{
  if(length(dim(Y)) < 3){ stop("Y is not a 3-mode array")}
  if(dim(Y)[1] != dim(Y)[2]){ stop("Y is not square in first two dimensions")}
  
  n <- dim(Y)[1]
  rest <- as.matrix(expand.grid(lapply(dim(Y)[-c(1,2)], function(x) 1:x)))   # combinations of all dimensions beyond first two
  r <- length(dim(Y)) - 2   # number of additional dimensions
  first2 <- cbind(rep(1:n, times=nrow(rest)), rep(1:n, times=nrow(rest)))
  rest <- matrix(rep(rest, each=n), ncol=r)
  
  return(cbind(first2, rest))
}


# ROC curve of predictions p, labels l, and number of iterations n
roc_curve <- function(p, l, n=1000)
{
  if(length(p) != length(l)){ stop("Lengths of predictions and labels must be the same.")}
  s <- sort(p)[seq(1, length(p), length.out=n)]
  s <- c(0,s,1)
  tpr <- fpr <- dA <- rep(NA, length(s))
  pos <- sum(l==1)
  neg <- sum(l==0)
  for(i in 1:length(s)){
    cut <- s[i]
    yhat <- 1*(p>=cut)
    tpr[i] <- sum(yhat == 1 & l == 1)/pos  # true positive rate
    fpr[i] <- sum(yhat == 1 & l == 0)/neg  # false positive rate
    if(i>1){
      dA[i] <- (tpr[i] + tpr[i-1])/2*(fpr[i] - fpr[i-1])
    }
  }
  auc <- sum(-dA, na.rm=T)  # area under cruve
  return(list(tpr=tpr, fpr=fpr, auc=auc))
}

# simple, direct auc
#  https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
simple_roc_auc <- function(p,l)
{
  x = p;  y = l
  x1 = x[y==1]; n1 = length(x1); 
  x2 = x[y==0]; n2 = length(x2);
  r = rank(c(x1,x2))  
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n2)
  auc
}

# simple, direct auc of precision/recall curve
simple_pr_auc <- function(p,l)
{
  if(sum(l)==0){stop("Labels are all zero. Need at least one 1. ")}
  
  l <- l[order(p, decreasing=T)]
  p <- p[order(p, decreasing=T)]
  tp <- cumsum(l)
  n <- length(l)
  n1 <- sum(l)
  n0 <- n - n1
  np <- 1:n
  
  prec <- tp/np
  # prec[is.infinite(prec)] <- 0
  rec <- tp/n1
  
  prec <- prec[order(rec)]
  rec <- rec[order(rec)]
  
  dx <- c(rec[1], rec[2:n] - rec[1:(n-1)])
  if(rec[1] !=0 ){
    prec_traps <- c(prec[1], .5*(prec[1:(n-1)] + prec[2:n]))
  } else {
    prec_traps <- c(prec[1]/2, .5*(prec[1:(n-1)] + prec[2:n]))
  }
  
  auc <- sum(prec_traps*dx)
  auc
}



# Precision-recall curve of predictions p, labels l
#   n does nothing now but is a legacy input for the old function
pr_curve <- function(p, l, n=1000)
{
  # old
  
  # if(length(p) != length(l)){ stop("Lengths of predictions and labels must be the same.")}
  # s <- sort(p)[seq(1, length(p), length.out=n)]
  # s <- c(0,s,1)
  # prec <- rec <- dA <- rep(NA, length(s))
  # pos <- sum(l==1)
  # neg <- sum(l==0)
  # for(i in 1:length(s)){
  #   cut <- s[i]
  #   yhat <- 1*(p>=cut)
  #   rec[i] <- sum(yhat == 1 & l == 1)/pos  # true positive rate
  #   prec[i] <- sum(yhat == 1 & l == 1)/sum(yhat==1) # false positive rate
  #   if(i>1){
  #     dA[i] <- (prec[i] + prec[i-1])/2*(rec[i] - rec[i-1])
  #   }
  # }
  # auc <- sum(-dA, na.rm=T)  # area under cruve
  
  #new
  d <- cbind(p,l)
  d <- d[order(d[,1], decreasing = T),]   # data ordered by predicted probability
  
  n1 <- sum(l)   # number of 1's
  tp <- cumsum(d[,2])  # true positives
  # fp <- 1:length(l) - tp 
  
  rec <- tp/n1
  prec <- tp/( 1:length(l) )  #  (tp + fp)
  auc <- trap_int(rec, prec)
  
  return(list(prec=prec, rec=rec, auc=auc))
}


# Trapezoidal integral
trap_int <- function(x,y)
{
  # y <- y[order(x)]
  # x <- x[order(x)]
  
  dx <- x[-1] - x[-length(x)]
  my <- sapply(2:length(y), function(i) mean(c(y[i], y[i-1])))
  
  sum(dx*my)
}



# Error rates of interest at given fraction of 1's (f)
error_rates <- function(p,l,f)
{
  # input massaging
  if(sum(l)==0){stop("Labels are all zero. Need at least one 1. ")}
  if(length(l)!=length(p)){stop("Lengths of predictions and labels must be the same")}
  f <- as.numeric(f)[1]     
  
  fracs <- (1:length(p))/length(p)  # all fraction 
  i <- which.min(abs(fracs - f))   # index of closest fraction
  t <- sort(p, decreasing=T)[i] - 1e-10
  
  # summaries
  n1 <- sum(l)
  n0 <- length(l) - n1
  tp <- sum(l[p > t])  # true positives
  fp <- sum(p > t) - tp  # false positives
  fn <- sum(l[p <= t])   # false negatives
  
  # error rates and data
  prec <- tp/(tp + fp)
  rec <- tp/n1
  fpr <- fp/n0
  fnr <- fn/n1
  fdr <- 1 - prec
  mean1 <- mean(p[l == 1])
  mean0 <- mean(p[l == 0])
  
  return(list(prec=prec, rec=rec, fpr=fpr, fnr=fnr, fdr=fdr, mean1=mean1, mean0=mean0))
}





#' Matricization
#' 
#' Matricize an array
#'
#' This functions matricizes an array along a given mode. 
#'
#' @param A an array
#' @param k a mode of \code{A}, along which to matricize
#' @keywords arrays matrices
#' @export
#' @examples
#' A<-rsan(4,3,2)
#' mat(A,2) 
mat<-function(A,k)
{
  Ak<-t(apply(A,k,"c"))
  if(nrow(Ak)!=dim(A)[k])  { Ak<-t(Ak) }
  Ak
}


#' Array-matrix product
#'
#' Multiply an array by a matrix along a given mode
#'
#' This function multiplies a matricized array by another 
#' matrix, and then reforms the result into a new array. 
#'
#' @param A a real valued array 
#' @param M a real matrix
#' @param k an integer, a mode of \code{A}
#' @author Peter Hoff
#' @keywords arrays
#' @export
#' @examples
#' A<-rsan(c(5,4,3))
#' B<-rsan(c(2,5))
#' amprod(A,B,1)
amprod<-function(A,M,k)
{
  K<-length(dim(A))
  AM<-M%*%mat(A,k)
  AMA<-array(AM, dim=c(dim(M)[1],dim(A)[-k]) )
  aperm(AMA,  match(1:K,c(k,(1:K)[-k]) ) )
}



# Multipartite BLIN mean 
blin_mean_multi <- function(X,A,modes=1:length(A))
{
  XA <- Reduce("+", lapply(modes, function(z) amprod(X[[z]], A[[z]], z)))
  XA
}


# Check if all elements in a vector are equal
vec.all.equal <- function(x) {diff(range(x)) < .Machine$double.eps ^ 0.5}




# Old or future functions below


# # Function to update W, iterating between W and WW^T updates
# update_W_symmetric <- function(D, Y, U, V, W)
# {
#   mat1 <- crossprod( crossprod(V, D), crossprod(U, D))
#   DD <- crossprod(D)
#   WW <- tcrossprod(W)
#   mat2 <- DD %*% WW %*% W 
#   mat2 <- mat2 + WW %*% DD %*% W
#     
#   W <- solve(t(mat1) + mat1, mat2)
#   return(W)
# }
# 
# 
# update_MLE_symmetric <- function(D, Y, U, V, W)
# {
#   # Call W update from within
# }


# # Given A and B and data, predict the future Y
# # return future Y slices
# # first 3 dimensions of X and D must be the same, and greater than Y
# predict_Y_forum <- function(A,B,Y,D)
# {
#   
#   # Check if dimensions agree
#   if(sum(dim(Y)[c(1,2)] != dim(D)[c(1,2)]) > 0){  stop("Dimenions of Y and D don't match")}
#   # if(sum(dim(D) != dim(X)[1:length(dim(D))]) > 0){  stop("Dimenions of D and X don't match")}
#   if(length(dim(Y)) < 3){ stop("Y is not an array")}
#   
#   # S <- dim(Y)[1]
#   # L <- dim(Y)[2]
#   t0 <- dim(Y)[3]
#   t1 <- dim(D)[3] - dim(Y)[3]
#   
#   if(t1 <1){stop("Need more time periods in D than in Y")}
#   
#   trange <- (t0+1):(t0+t1)   # range of time periods to predict upon
#   
#   # Yout <- drop(amprod(X[,,trange,], matrix(beta, nrow=1), length(dim(X[,,trange,]))))   # multiply Xbeta to start
#   Yout <- amprod(D[,,trange], A, 1) + amprod(D[,,trange], t(B), 2)   # add AD and DB
#   
#   return(Yout)
# }







# # Perform cross-validation using Precision recall AUC as loss measure
# cv_pr_old <- function(Y, D, X, outdir, use_cov=T, type="biten", penalty=1, sparsedata=F, writeXreg=T, readXreg=T, 
#                       nameXreg="Xreg.txt", seed=NA, ncv=10, lambdas=10^seq(-5,-7, length.out=100), verbose=F, maxit=1e5, ncores=1)
# {
#   if(is.numeric(seed)){ set.seed(seed)}
#   if(! is.numeric(penalty)){ stop("Penalized regression method requires numeric penalty value (alpha") }
#   
#   dir.create(outdir, showWarnings = F)
#   
#   if(!sparsedata){
#     
#     stop("Not implemented for non-sparse data")
#     # # Check sizes
#     # if(length(dim(Y)) != 3 ){ stop("Y is not a 3-mode array") }
#     # if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
#     # if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
#     # 
#     # # Find sizes, assuming dim(D) = dim(Y)
#     # S <- nrow(D[,,1])
#     # L <- ncol(D[,,1])
#     # tmax <- dim(D)[3]
#     # 
#     # # Build X matrix
#     # Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata)  # build design matrix
#     # 
#     # if(use_cov){
#     #   p <- dim(X)[4]
#     # }
#     
#   } else if (sparsedata){
#     
#     # Check if i,j,t in column names
#     if( !("i" %in% names(Y)) | !("j" %in% names(Y)) | !("t" %in% names(Y))){stop("Y must have column names i,j, and t")}
#     
#     # Find sizes
#     S <- max(D$i)
#     L <- max(D$j)
#     tmax <- max(D$t)
#     
#     # Sort Y
#     rows <- Y$i + (Y$j-1)*S + (Y$t-1)*S*L   # unfolded indices
#     indices <- order(rows)
#     Y <- Y$ratification_year[indices]  # 1's and zeros, in order
#     
#     # Build X matrix (full)
#     if(readXreg & nameXreg %in% list.files()){
#       cat("\n Reading in Xreg \n")
#       load(nameXreg)  # should be already sorted/subsetted
#     } else {
#       if(readXreg){cat("\n nameXreg file not found. Generating Xreg \n")}
#       Xreg <- build_design_additive(D, X, type=type, use_cov=use_cov, sparsedata = sparsedata)  # build design matrix
#       cat("Design built, dimensions ", dim(Xreg), "\n")
#       rows <- rows[indices]  # sorted rows to keep
#       Xreg <- Xreg[rows, ]  # subset
#       Xreg <- Xreg[,-(S^2 + L^2 + 1)]  # remove intercept
#     }
#     
#     if(use_cov){
#       p <- ncol(X) - 3 - 1   # -3 for i,j,t, -1 for intercept removed
#     }
#     
#   } else { stop("sparsedata must be true/false")}
#   
#   
#   Xreg[which(is.na(Xreg), arr.ind=T)] <- 0   # set NAs to zero
#   if(writeXreg){  # write out
#     save(Xreg, file=nameXreg)
#   }
#   
#   # Bookkeep zero columns
#   remove <- which(colSums(Xreg) == 0)   
#   keep <- (1:ncol(Xreg))[-remove]
#   
#   # Cross-validate
#   lambdas <- sort(lambdas, decreasing = T)  # decreasing sequence
#   cvs <- 1:nrow(Xreg) %% ncv + 1
#   cvms <- matrix(NA, length(lambdas), ncv)
#   assign("last.warning", NULL, envir = baseenv())   # clear warnings
#   
#   if(ncores==1){
#     fitlist <- vector("list", ncv)
#     for(i in 1:ncv){
#       test <- which(cvs == i)
#       train <- (1:nrow(Xreg))[-test]
#       
#       fitlist[[i]] <- fit <- glmnet(Xreg[train,keep], y=as.factor(c(Y[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-10, standardize=F)
#       worked <- length(fit$lambda)
#       
#       save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
#       
#       if(verbose){
#         cat("\n************************************************\n")
#         cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
#         print(warnings())
#         cat("\n************************************************")
#         cat("\n")
#         assign("last.warning", NULL, envir = baseenv())
#       }
#     }
#     
#   } else if (ncores > 1 ){
#     writeLines(c(""), file.path(outdir, "log.txt"))
#     
#     mcoptions <- list(preschedule=FALSE, set.seed=T)
#     fitlist <- foreach(i=1:ncv, .options.multicore=mcoptions, .packages=c("glmnet") ) %dopar% {  # .combine =cbind
#       test <- which(cvs == i)
#       train <- (1:nrow(Xreg))[-test]
#       
#       fit <- glmnet(Xreg[train,keep], y=as.factor(c(Y[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit, lambda.min.ratio=1e-12, standardize=F)
#       save(fit, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
#       worked <- length(fit$lambda)
#       sink(file.path(outdir, "log.txt"), append=TRUE)   # write out to log file
#       cat("\n************************************************\n")
#       cat("fit cv", i, "of", ncv, ", nlambda", worked, "\n")
#       print(warnings())
#       cat("\n************************************************")
#       cat("\n")
#       fit
#     }
#     
#   } else {stop("ncores must be numeric >= 1")}
#   
#   for(i in 1:ncv){
#     test <- which(cvs == i)
#     train <- (1:nrow(Xreg))[-test]
#     
#     fit <- fitlist[[i]]
#     #glmnet(Xreg[train,keep], y=as.factor(c(Y[train])), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, maxit=maxit)
#     worked <- length(fit$lambda)
#     
#     if(worked > 1){   # if fit worked
#       Yhats <- predict(fit, newx=Xreg[test,keep], type="response")
#       
#       pr_temp <- rep(NA, ncol(Yhats))
#       for(j in 1:ncol(Yhats)){
#         # pr_temp[j] <- pr_curve(Yhats[,j], Y[test], n=2000)$auc
#         # pr_temp[j] <- roc_curve(Yhats[,j], Y[test], n=min(length(Y[test]), 2500))$auc
#         # pr_temp[j] <- simple_roc_auc(Yhats[,j], Y[test])
#         pr_temp[j] <- simple_pr_auc(Yhats[, j], Y[test])
#         # pr_temp[j] <- auc(roc(Yhats[,j], Y[test]))
#       }
#       
#       cvms[match(fit$lambda, lambdas), i] <- pr_temp   # save lambda values
#       # save(fit, pr_temp, file=file.path(outdir, paste0("pr_cv", i, ".RData")))
#     }
#     if(verbose){
#       cat("processed cv", i, "of", ncv, ", nlambda", worked, "\n")
#     }
#   }
#   
#   # Calculate minimum lambda value
#   mean_cvms <- apply(cvms, 1, mean, na.rm=T)
#   imin <- which(mean_cvms == max(mean_cvms, na.rm=T))
#   lambda_min <- lambdas[imin]
#   
#   # Perform full fit and extract coefficients for lambda_min
#   fit <- glmnet(Xreg[,keep], y=as.factor(c(Y)), alpha=penalty, family='binomial', intercept=T, lambda=lambdas, lambda.min.ratio=1e-10, standardize=F)
#   worked <- length(fit$lambda) 
#   if(verbose){
#     cat("full fit done, nlambda", worked, "\n")
#   }
#   if(worked > 1){
#     Yhats <- predict(fit, newx=Xreg[,keep], type="response")
#     cvm_full <- rep(NA, ncol(Yhats))
#     for(j in 1:ncol(Yhats)){
#       cvm_full[j] <- simple_pr_auc(Yhats[,j], Y) # pr_curve(Yhats[,j], Y, n=min(length(Y), 2500))$auc
#     }
#     Yhat <- predict(fit, newx=Xreg[,keep], type="response", s=lambda_min)
#     coefs <- coef(fit, s=lambda_min)
#     
#     # inflate coefs with NAs for non-estimable quantities
#     oldcoefs <- coefs   
#     coefs <- rep(NA, ncol(Xreg))
#     coefs[keep] <- oldcoefs[-1]   # add in NAs for un-estimable quantities, -1 for intercept
#     coefs <- c(oldcoefs[1], coefs)  # intercept
#     
#     A <- matrix(coefs[1:S^2 + 1], S, S)   # + 1 for intercept
#     B <- matrix(coefs[S^2 + 1 + 1:L^2], L, L)
#     if(use_cov){
#       beta <- coefs[c(1,1+S^2 + L^2 + 1:p)]
#       beta <- matrix(beta, nrow=length(beta))
#     } else { beta <- NA }
#   } else {
#     fit <- A <- B <- Yhat <- NA
#   }
#   
#   return(list(A=A, B=B, beta=beta, Yhat=Yhat, fit=fit, Yreg=Y, Xreg=Xreg, Xkeep=keep, Xremove=remove, cvms=cvms, cvm_full=cvm_full, lambda_min=lambda_min))
# }
#
#
#
#
# # Given A and B and data, predict the future Y
# # return future Y slices
# # first 3 dimensions of X and D must be the same, and greater than Y
# predict_Y_forum <- function(A,B,Y,D)
# {
#   
#   # Check if dimensions agree
#   if(sum(dim(Y)[c(1,2)] != dim(D)[c(1,2)]) > 0){  stop("Dimenions of Y and D don't match")}
#   # if(sum(dim(D) != dim(X)[1:length(dim(D))]) > 0){  stop("Dimenions of D and X don't match")}
#   if(length(dim(Y)) < 3){ stop("Y is not an array")}
#   
#   # S <- dim(Y)[1]
#   # L <- dim(Y)[2]
#   t0 <- dim(Y)[3]
#   t1 <- dim(D)[3] - dim(Y)[3]
#   
#   if(t1 <1){stop("Need more time periods in D than in Y")}
#   
#   trange <- (t0+1):(t0+t1)   # range of time periods to predict upon
#   
#   # Yout <- drop(amprod(X[,,trange,], matrix(beta, nrow=1), length(dim(X[,,trange,]))))   # multiply Xbeta to start
#   Yout <- amprod(D[,,trange], A, 1) + amprod(D[,,trange], t(B), 2)   # add AD and DB
#   
#   return(Yout)
# }






# Function to build regression design matrix for additive models
# D is array to regress Y on for A and B, additive models only!
# X is array of covariates to regress Y upon for beta
# type is "biten" or "sadd"
# use_cov is boolean flag for using covariates, i.e. X/betas, or not
#
#
# Returns design matrix
#    includes intercept if there is one in X (i.e. all covariates in X are included)
#    returns S*L*tmax \times S^2 + L^2 + ncol(X) matrix, with columnwise-vectorized order of rows
build_design_additive_multi1 <- function(D, X=NA)
{

  Ddims <- sapply(D, dim)
  if(!all(apply(Ddims, 1, vec.all.equal))){ stop("Not all lagged matrices in D have the same dimension")}
  Ddims <- dim(D[[1]])

  if(length(Ddims) > 4 | length(Ddims) < 3){ stop("Only coded for 2 or 3 modes of interest") }

  X <- NA
  type <- "biten"
  use_cov <- FALSE
  sparsedata=FALSE
  write=FALSE
  response=NA


  # Find sizes
  tmax <- Ddims[length(Ddims)]
  ms <- Ddims[-length(Ddims)]

  # Initialize X
  Xreg <- sparseMatrix(i=1, j=1, x=0, dims=c(tmax*prod(ms), sum(ms^2)))   # initialize X

  if(length(Ddims) == 4){
    R <- Ddims[3]
    for(t in 1:tmax){  # vec A then vec B and so on
      Xtemp <- NULL
      for(r in 1:Ddims[3]){
        temp <- cbind( kronecker(t(D[[1]][,,r,t]), diag(Ddims[1]) ), kronecker( diag(Ddims[2]), D[[2]][,,r,t]) )
        Xtemp <- rbind(Xtemp, temp)
      }
      Rcols <- sapply(1:Ddims[3], function(z) c(D[[3]][,,z,t]))
      Xtemp <- cbind(Xtemp, kronecker( diag(Ddims[3]), Rcols))

      Xreg[ 1:prod(Ddims[1:3]) + (t - 1)*prod(Ddims[1:3]), ] <- Xtemp
    }

  } else if (length(Ddims) == 3) {
    for(t in 1:tmax){  # vec A then vec B and so on
      Xtemp <- NULL
      temp <- cbind( kronecker(t(D[[1]][,,t]), diag(Ddims[1]) ), kronecker( diag(Ddims[2]), D[[2]][,,t]) )
      Xtemp <- rbind(Xtemp, temp)

      # Rcols <- sapply(1:R, function(z) c(D3[,,z,t]))
      # Xtemp <- cbind(Xtemp, kronecker(Jr, Rcols))

      Xreg[ 1:prod(Ddims[1:2]) + (t - 1)*prod(Ddims[1:2]), ] <- Xtemp
    }

  }

  return(Xreg)
}





# Compute off-diagonal and on-diagonal MSEs 
ABcompare <- function(A,B,Atrue,Btrue,model,modeltrue, computediag = TRUE)
{
  A0 <- A  ;  diag(A0) <- NA  ;  A0 <- A0 / sum(A0, na.rm=T)
  B0 <- B  ;  diag(B0) <- NA  ;  B0 <- B0 / sum(B0, na.rm=T)
  
  Atrue0 <- Atrue  ;  diag(Atrue0) <- NA  ;  Atrue0 <- Atrue0 / sum(Atrue0, na.rm=T)
  Btrue0 <- Btrue  ;  diag(Btrue0) <- NA  ;  Btrue0 <- Btrue0 / sum(Btrue0, na.rm=T)
  
  mseA <- mean( ( A0 - Atrue0 )^2, na.rm=T)
  mseB <- mean( ( B0 - Btrue0 )^2, na.rm=T)
  
  if(computediag){
    diagAB <- outer(diag(A), diag(B), diagfun(model))
    diagABtrue <- outer(diag(Atrue), diag(Btrue), diagfun(modeltrue))
    
    mseAB <- mean( (diagAB - diagABtrue)^2, na.rm=T)
  } else {
    mseAB <- TRUE
  }
  return(list(mseA=mseA, mseB=mseB, mseAB=mseAB))
}

# Compute off-diagonal MSEs using linear regression
ABmse <- function(A,B,Atrue,Btrue)
{
  diag(Atrue) <- NA  
  diag(Btrue) <- NA  
  if(length(dim(A)) > 2){
    k <- dim(A)[3]
  } else {
    k <- 1
  }
  
  y <- c(rep(c(Atrue), times=k), rep(c(Btrue), times=k))
  X <- matrix(0, length(y), 2)
  X[1:length(A), 1] <- c(A)
  X[1:length(B) + length(A),2] <- c(B)
  
  fit <- lm(y ~ X - 1)
  mseAB <- mean(resid(fit)^2)
  
  return(mseAB)
}

# Function to compute diagonal based on model
diagfun <- function(model)
{
  if(strtrim(model, 3) == "bit" | strtrim(model, 3) == "sad"){
    return("+")
  } else if (strtrim(model, 3) == "bil" ){
    return("*")
  } else {
    return(NA)
  }
}


