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



# Additive model solution via regression for both BiTEN and smoothed additive
# Perform regression of Y onto D for a given type
# D is a 3-mode array of covariates
# Y is a 3-mode array of responses of the same size as D
# X is a 4-mode array of covariates
# type is "sadd" or "biten", type of regression to perform
# use_cov is a boolean flag to use covariates/X or not
# test_rank is a flag to test the rank of X (can be slow based on size of D, Y)
# returns A, B, beta, Yhat (array)
additive_regression <- function(Y, D, X, type="biten", use_cov=T, test_rank=F, penalty=NA, whichlambda="1se")
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

# 
# # Fit BLIN model using iterative updates of A and B
# fit_MLE_blin_iterative <- function(Y, D, X, verbose=T, printout=10, tol=1e-10, init="random", sigma_init=1, use_cov=F, maxit=1e9)
# {
#   
#   if(use_cov){ stop("Not including covariates yet")}
#   
#   # Get sizes and check
#   if(sum(dim(Y) != dim(D)) > 0){  stop("Dimenions of Y and D don't match")}
#   if(sum(dim(Y) != dim(X)[1:length(dim(Y))]) > 0 & use_cov){  stop("Dimenions of Y and X don't match")}
#   if(length(dim(Y)) < 3){ stop("Y is not an array")}
#   
#   i_na <- is.na(Y)  # save NAs
#   Y[i_na] <- 0  # set to zero
#   
#   S <- dim(Y)[1]
#   L <- dim(Y)[2]
#   tmax <- dim(Y)[3]
#   if(use_cov){
#     p <- dim(X)[4]
#   }
#   
#   
#   # matrices to regress upon
#   # Xa <- matrix(0, tmax*S*L, S^2)
#   # Xb <- matrix(0, tmax*S*L, L^2)
#   # for(t in 1:tmax){
#   #   Xa[1:(S*L) + (t-1)*S*L, ] <- kronecker(t(D[,,t]), diag(S))
#   #   Xb[1:(S*L) + (t-1)*S*L, ] <- kronecker(diag(L), D[,,t])
#   # }
#   # XXa <- crossprod(Xa)
#   # XXb <- crossprod(Xb)
#   XXa <- solve( Reduce("+", lapply(1:tmax, function(t) tcrossprod(D[,,t],D[,,t])) ) )
#   XXb <- solve( Reduce("+", lapply(1:tmax, function(t) crossprod(D[,,t],D[,,t])) ) )
#   
#   
#   
#   # Initialize
#   if(strtrim(init,3) == "ran") {
#     if(verbose == T){
#       cat("Randomly initializing A,B, beta \n")
#     }
#     A <- matrix(rnorm(S^2, 0, sigma_init), S, S)
#     BT <- matrix(rnorm(L^2, 0, sigma_init), L, L)
#     if(use_cov){
#       beta <- rnorm(p+1, 0, sigma_init)
#     }
#     
#   } else if (strtrim(init, 3) == "reg") {
#     if(verbose == T){
#       cat("Initializing A,B,beta with regression \n\n")
#     }
#     
#     avec <- solve(XXa, crossprod(Xa, Y))
#     bvec <- solve(XXb, crossprod(Xb, Y))
#     
#     A <- matrix(avec, S, S)
#     BT <- t(matrix(bvec, L, L))
#     
#     
#   } else if (strtrim(init, 1) == "I") {
#     if(verbose == T){
#       cat("Initializing A, B as identity \n\n")
#     }
#     if(sum(dim(Y) != dim(D)) > 0){  stop("Cannot initialize identity when dimensions of Y and D are not the same")}
#     
#     A <- diag(S)
#     BT <- diag(L)
#     
#   } else { stop("Invalid initialization type")}
#   
#   
#   Ainit <- A
#   BTinit <- BT
#   
#   
#   # if(use_cov){
#   #   betaOLS <- lm(c(Y) ~ -1 + t(mat(X, 4)) )$coef  # best fit ignoring A,B structure
#   # } else { betaOLS <- NA}
#   
#   # Find optimal values
#   change <- 100
#   count <- 0
#   
#   while(change > tol & count < maxit){
#     
#     if(count == 0) {   # save initial LL and beta
#       # beta_init <- beta
#       LLinit <- LL <- -length(Y)*log( sum( ( Y - amprod(D, A, 1) - amprod(D, BT, 2) )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
#     }
#     
#     # Ytilde <- Y - amprod(D, A, 1)
#     # bvec <- solve(XXb, crossprod(Xb, Ytilde))
#     # BTnew <- t(matrix(bvec, L, L))
#     # 
#     # Ytilde <- Y - amprod(D, BTnew, 2)
#     # avec <- solve(XXa, crossprod(Xa, Ytilde))
#     # Anew <- matrix(avec, S, S)
#     
#     Ytilde <- Y - amprod(D, A, 1)
#     XYb <- Reduce("+", lapply(1:tmax, function(t) crossprod(D[,,t], Ytilde[,,t])))
#     BTnew <- t( XXb %*% XYb )
#     
#     Ytilde <- Y - amprod(D, BTnew, 2)
#     XYa <- Reduce("+", lapply(1:tmax, function(t) tcrossprod(Ytilde[,,t], D[,,t])))
#     Anew <- XYa %*% XXa
#     
#     
#     changeAB <- max(abs(c(c(A - Anew), c(BT-BTnew))))  # doesn't make much difference stopping A,B vs using LL
#     A <- Anew
#     BT <- BTnew
#     
#     Yhat <- amprod(D, A, 1) + amprod(D, BT, 2)  # new estimate
#     Y[i_na] <- Yhat[i_na]  # update NAs
#     
#     # LLnew <- -sum( (Ytilde - (amprod(D, A, 1) + amprod(D, BT, 2)) )^2)
#     LLnew <- -length(Y)*log( sum( (Y - Yhat )^2) ) - length(Y)   #+ 2*k + 2*m + length(beta)
#     change <- abs(LLnew - LL) 
#     LL <- LLnew
#     
#     count <- count + 1
#     if(count%%printout == 0 & verbose == T){
#       cat("Iteration: ", count, " \t Criterion: ", change, "\t 2LL: ", LL,"\n")
#     }
#   }
#   
#   if(verbose == T){
#     cat("\n************************************************ \n")
#     
#     # cat("True 2LL*tau^2: \t", LLtrue, "\n")
#     cat("Initial 2LL: \t", LLinit, "\n")
#     cat("Final 2LL: \t", LL, "\n")
#     
#     cat("\n************************************************ \n \n")
#     
#     # cat("OLS beta coefficients are: \t\t", betaOLS, "\n")
#     # cat("Est. beta coefficients are: \t\t", beta, "\n")
#   }
#   
#   return(list(A=A, B=t(BT), Yhat = Yhat, LLt2 = LL, LLt2_init=LLinit))   # beta= beta, betaOLS = betaOLS, 
# }




##############################
###  Prediction functions  ###
##############################

# These functions make predictions based on the results of previous fits


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








##########################
###  Update Functions  ###
##########################

# These functions perform updates for MLE block coordinate descents



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


# # HOFF multiplicative bilinear model: Update for block coordinate descent 
# # D is matrix of lagged Y (covariates)
# # Y is oservations
# # for multiple time periods
# update_MLE_bilinear <- function(D, Y, A, B)
# {
#   
#   t <- dim(D)[3]   # number of time slices
#   
#   # A update
#   DBY <- Reduce("+", lapply(1:t, function(x) tcrossprod(D[,,x ] %*% B, Y[,,x] ) ))    # D B Y^T
#   DBBD <- Reduce("+", lapply(1:t, function(x) tcrossprod( D[,,x] %*% B)  ))  # D B B^T D^T
#   A <- t( solve( DBBD, DBY ))
#   
#   # B update
#   DAY <- Reduce("+", lapply(1:t, function(x) crossprod(A %*% D[,,x ], Y[,,x])  ))    # D^T A^T Y
#   DAAD <- Reduce("+", lapply(1:t, function(x) crossprod( A %*% D[,,x]  )))    # D^T A^T A D
#   B <- ( solve( DAAD, DAY ))   # no transpose! 
#   
#   return(list(A=A, B=B))
# }
# 
# 
# # Update for block coordinate descent 
# # asymmetric A and B
# # D is matrix of lagged Y (covariates)
# # Y is observations
# # for multiple time periods
# # type is sadd or biten additive model
# update_MLE_additive <- function(A, B, D, DDT, DTD, DYT, DTY, type="biten")
# {
#   S <- dim(D)[1]
#   L <- dim(D)[2]
#   tmax <- dim(D)[3]
#   
#   if(type == "sadd"){
#     Jl <- matrix(1,L,L)
#     Js <- matrix(1,S,S)
#     DBD <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x] %*% Jl, Js %*% D[,,x] %*% B)))    # D J B^T D^T J
#     A <- t(solve(DDT, DYT - DBD))
#     DAD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], Js %*% A %*% D[,,x] %*% Jl)))    # D^T J A D J
#     B <- (solve(DTD, DTY - DAD))   # no transpose! Use B instead of B^T
#     
#   } else if (type == "biten") {
#     DBD <- Reduce("+", lapply(1:tmax, function(x) tcrossprod(D[,,x], D[,,x] %*% B)))    # D B^T D^T
#     A <- t(solve(DDT, DYT - DBD))
#     DAD <- Reduce("+", lapply(1:tmax, function(x) crossprod(D[,,x], A %*% D[,,x])))    # D^T A D 
#     B <- (solve(DTD, DTY - DAD))  # no transpose! Use B instead of B^T
#     
#   } else { stop( "Invalid type " ) }
#   
#   return(list(A=A, B=B))
# }





##########################
###  Helper Functions  ###
##########################

# These functions read data, plot data, build useful matrices, and perform lower-level operations


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
arcplot_fm2 <- function(mplot, groups, colors, linespace=.05, maxwt=NULL, maxdeg=NULL, horizontal=T)   # positive and negative together
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
  
  arcplot(edgelist,  cex.labels = 0.8,  # ordering = new_ord,  #labels = vlabels,
          lwd.arcs = 14*weights, 
          lwd.nodes=1, show.nodes = TRUE, col.nodes = NA, bg.nodes = nodecolors, labels=rownames(mplot),
          cex.nodes = 2.4*nodesizes + .6, pch.nodes = 21,  line = 0.1 +linespace,
          col.arcs = arccolors, col.labels="gray50", above=above, 
          horizontal=horizontal)
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


