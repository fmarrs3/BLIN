# BiTEN Project: MLE estimation
# Frank Marrs
# 09/13/16
#
# This file runs MLE estimation of A and B matrices on Y and D based on inputs
# 
# Run 10-fold cross-validation 
#

# Single fit of a given fit type in "type" input
misspec_model_fit <- function(type, Y, X, D, k, m, use_cov, gen_type, S, L, tmax, tau, sim, seed, rankA, rankB, Atrue, Btrue)
{
  
  # Inputs
  # type <- "biten"  # type of model to fit, "bilinear", "biten", or "sadd"
  # use_cov <- F   # should covariates (X/betas) be included in fit?
  if(type == "biten"){
    init <- "ran"   # "random", "regression", or "I" (identity) initialization?
  } else {
    init <- "I"
  }
  
  
  # Set up cross-validation
  set.seed(1)   # set seed to use same partition each time
  n_cv <- min(dim(Y)[3], 10)   # number of cross-validations
  i_cv <- sample(1:dim(Y)[3]%%n_cv +1)   # cross-validation indices
  # Ypred <- Y*0   # prediction values initialization
  
  errors <- matrix(NA, n_cv, 6)   # out-of-sample and in-sample errors 
  fittime <- rep(NA, n_cv)
  pA <- rA <- pB <- rB <- rep(NA, n_cv)
  
  
  # Run MLE fit with CV
  for(i in 1:n_cv){
    test <- which(i_cv == i)
    train <- which(i_cv != i)
    
    if(use_cov){
      Xfit <- X[,,train,, drop=F]
      Xtest <- array(X[,,test,, drop=F], c(dim(X)[1:2], length(test), dim(X)[4])) 
    } else {
      Xfit <- Xtest <- 0
    }
    
    t0 <- proc.time()
    if(type == "biten"){
      results <- fit_MLE_array(Y[,,train], D[,,train], Xfit, k, m, verbose=F, init=init, sigma_init=1, use_cov=use_cov)
      A <- results$A   # fit of matrix A
      B <- results$B   # fit of matrix B
      beta <- results$beta
      Yhat <- results$Yhat   # estimated Y for given data
      # convoluted entry of X and D so they don't drop in size
      Ypred <- drop(predict_Y(A, B, beta,
                                 Xtest, 
                                 D[,,test,drop=F], 
                                 type=type, use_cov=use_cov))  # save prediction on test set
      obs_sparsityA <- obs_sparsityB <- PRA <- PRB <- NA
      comp <- ABcompare(A,B,Atrue,Btrue,type,gen_type)
      mseA <- comp$mseA  ;  mseB <- comp$mseB  ;  mseAB <- comp$mseAB

      
    # } else if (type == "bilinear" & use_cov){
    #   # results <- fit_MLE_array_bilinear(Y[,,train], D[,,train], Xfit, verbose=F, init=init, sigma_init=1, use_cov=use_cov)
    #   # A <- results$A   # fit of matrix A
    #   # B <- results$B   # fit of matrix B
    #   # Yhat <- results$Yhat   # estimated Y for given data
    #   # 
    #   results <- mlm.ALS.cov(Y[,,train], D[,,train], Xfit, verbose=F, tol = 1e-5, imax=1000)
    #   Ypred <- drop(tprod(D[,,test], results$B) + drop(amprod( Xtest, results$beta, 4)) )
    #   Yhat <- tprod(D[,,train], results$B) + drop(amprod(Xfit, results$beta, 4))
    #   A <- results$B[[1]]
    #   B <- t(results$B[[2]])
      
    } else if (type == "bilinear" ) {   # & !use_cov 
      B <- mlm.ALS(Y[,,train], D[,,train], verbose=F, tol = 1e-6, imax=1000)
      Ypred <- drop(tprod(D[,,test,drop=F], B[[1]]))
      Yhat <- tprod(D[,,train], B[[1]])
      A <- B[[1]][[1]]
      B <- t(B[[1]][[2]])
      # results <- fit_MLE_array_bilinear(Y[,,train], D[,,train], Xfit, verbose=F, printout=100, tol=1e-6, init="I", sigma_init=1, use_cov=use_cov, imax=1e4)
      # Yhat <- results$Yhat
      # A <- results$A
      # B <- results$B
      # beta <- results$beta
      # if(use_cov){ Xb <- drop(amprod( Xtest, results$beta, 4)) } else { Xb <- 0 }
      # Ypred <- drop(tprod(D[,,test], list(results$A, t(results$B)))) + Xb
      obs_sparsityA <- obs_sparsityB <- PRA <- PRB <- NA
      comp <- ABcompare(A,B,Atrue,Btrue,type,gen_type)
      mseA <- comp$mseA  ;  mseB <- comp$mseB  ;  mseAB <- comp$mseAB
      
      
    } else if (type == "sadd"){
      # results <- fit_MLE_array_additive(Y[,,train], D[,,train], Xfit, verbose=F, init=init, sigma_init=1, use_cov=use_cov)
      results <- additive_regression(Y[,,train], D[,,train], Xfit, type=type, use_cov=use_cov, test_rank=F)
      A <- results$A   # fit of matrix A
      B <- results$B   # fit of matrix B
      Yhat <- results$Yhat   # estimated Y for given data
      # convoluted entry of X and D so they don't drop in size
      Ypred <- drop(predict_Y_regression(A, B, beta,
                                 Xtest, 
                                 D[,,test,drop=F],
                                 type=type, use_cov=use_cov)) # save prediction on test set
      obs_sparsityA <- obs_sparsityB <- PRA <- PRB <- NA
      comp <- ABcompare(A,B,Atrue,Btrue,type,gen_type)
      mseA <- comp$mseA  ;  mseB <- comp$mseB  ;  mseAB <- comp$mseAB
      
      
    } else if (type == "biten_full"){
      results <- additive_regression(Y[,,train], D[,,train], Xfit, test_rank = F, type="biten", use_cov=use_cov)
      A <- results$A   # fit of matrix A
      B <- results$B   # fit of matrix B
      Yhat <- results$Yhat   # estimated Y for given data
      # convoluted entry of X and D so they don't drop in size
      Ypred <- drop(predict_Y_regression(A, B, beta,
                                 Xtest, 
                                 D[,,test,drop=F], 
                                 type=type, use_cov=use_cov))  # save prediction on test set
      obs_sparsityA <- obs_sparsityB <- PRA <- PRB <- NA
      comp <- ABcompare(A,B,Atrue,Btrue,type,gen_type)
      mseA <- comp$mseA  ;  mseB <- comp$mseB  ;  mseAB <- comp$mseAB
      
      
    } else if (type == "biten_sparse"){
      results <- additive_regression(Y[,,train], D[,,train], Xfit, test_rank = F, type="biten", use_cov=use_cov, penalty = 1, whichlambda="min")
      A <- results$A   # fit of matrix A
      B <- results$B   # fit of matrix B
      Yhat <- results$Yhat   # estimated Y for given data
      # convoluted entry of X and D so they don't drop in size
      Ypred <- drop(predict_Y_regression(A, B, beta,
                                         Xtest, 
                                         D[,,test,drop=F], 
                                         type=type, use_cov=use_cov))  # save prediction on test set
      obs_sparsityA <- mean(A!=0)
      obs_sparsityB <- mean(B!=0)
      A1 <- 1*(A!=0)
      B1 <- 1*(B!=0)
      Atrue1 <- 1*(Atrue != 0)
      Btrue1 <- 1*(Btrue != 0)
      pA[i] <- sum(A1 == Atrue1)/sum(A1)  # precision
      pB[i] <- sum(B1 == Btrue1)/sum(B1)  # precision
      rA[i] <- sum(A1 == Atrue1)/sum(Atrue1)  # recall
      rB[i] <- sum(B1 == Btrue1)/sum(Btrue1)  # recall
      comp <- ABcompare(A,B,Atrue,Btrue,type,gen_type)
      mseA <- comp$mseA  ;  mseB <- comp$mseB  ;  mseAB <- comp$mseAB
      
      
    } else if (type == "mu"){
      mean_model <- mean(Y[,,train])
      Yhat <- array(mean_model, dim(Y[,,train]))
      Ypred <- array(mean_model, dim(Y[,,test]))
      A <- mean_model
      B <- NA
      obs_sparsityA <- obs_sparsityB <- PRA <- PRB <- NA
      mseA <- mseB <- mseAB <- NA
      
    } else { stop("Invalid input type")}
    
    fittime[i] <- (proc.time() - t0)[3]
    insample <- c(sqrt(mean((Y[,,train] - Yhat)^2)), mean(abs(Y[,,train] - Yhat)))  # Save in-sample errors
    outsample <- c(sqrt(mean((Y[,,test] - Ypred)^2)), mean(abs(Y[,,test] - Ypred)))  # Save out-of-sample errors
    R2 <- 1 - c(mean((Y[,,test] - Ypred)^2) / mean(Y[,,test]^2) , mean((Y[,,train] - Yhat)^2) / mean(Y[,,train]^2) ) 
    errors[i,] <- c(outsample, insample, R2)  
  
  }
  
  # names(out) <- c("gen_type", "fit_type", "S", "L", "tmax", "tau", "sim", "seed", "rankA", "rankB", "init", "CV", "rmse_out", "made_out", "rmse_in", "made_in")
  out_new <- cbind(gen_type, type, S, L, tmax, tau, sim, seed, rankA, rankB, init, 1:n_cv, errors, fittime, use_cov, obs_sparsityA, obs_sparsityB, pA, rA, pB, rB, mseA, mseB, mseAB)
  return(list(out=out_new, Yhat=Yhat, Ypred=Ypred, A=A,B=B))
}


# Compute off-diagonal and on-diagonal MSEs 
ABcompare <- function(A,B,Atrue,Btrue,model,modeltrue)
{
  A0 <- A  ;  diag(A0) <- NA  ;  A0 <- A0 / sum(A0, na.rm=T)
  B0 <- B  ;  diag(B0) <- NA  ;  B0 <- B0 / sum(B0, na.rm=T)
  
  Atrue0 <- Atrue  ;  diag(Atrue0) <- NA  ;  Atrue0 <- Atrue0 / sum(Atrue0, na.rm=T)
  Btrue0 <- Btrue  ;  diag(Btrue0) <- NA  ;  Btrue0 <- Btrue0 / sum(Btrue0, na.rm=T)
  
  mseA <- mean( ( A0 - Atrue0 )^2, na.rm=T)
  mseB <- mean( ( B0 - Btrue0 )^2, na.rm=T)
  
  diagAB <- outer(diag(A), diag(B), diagfun(model))
  diagABtrue <- outer(diag(Atrue), diag(Btrue), diagfun(modeltrue))
  
  mseAB <- mean( (diagAB - diagABtrue)^2, na.rm=T)
  
  return(list(mseA=mseA, mseB=mseB, mseAB=mseAB))
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


