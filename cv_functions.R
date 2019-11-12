

#### Run simulation study
run_cv <- function(nsims=1e2,       # number of simulations
                   qs=c(0,.5,.9),   # fractions of zero entries
                   outdir=file.path(getwd(), "results"))  # where to write
{
  outdir_coef <- outdir
  dir.create(outdir, showWarnings = F)
  dir.create(outdir_coef, showWarnings = F)
  
  tau <- 1
  S <- 10
  L <- S
  use_cov <- FALSE
  gen_cov <- FALSE   # generate data with covariates?
  # k <- m <- 1   # ranks ofBiTEN fits
  nsims <- 1e2
  tmaxes <- c(10, 20, 50)
  nfold <- min(tmaxes)
  gen_types <- c("blin", "bilinear")
  burn <- max(tmaxes)
  R2target <- .75
  diagval <- 0
  # qs <- c(0, .5, .9)
  
  
  for(q in qs){ 
    
    # Find scaling values for A and B
    ks <- 2^seq(-2, 2, by=.05)
    temp <- cbind(ks, 0, 0)
    for(i in 1:length(ks)){
      
      k <- ks[i]
      
      seedA <- 1
      set.seed(seedA)
      A <- outer(rnorm(S), rnorm(S), "*") /3.5
      diag(A) <- NA
      A[which(abs(A) < quantile(abs(A)[!is.na(abs(A))], q))] <- 0
      diag(A) <- diagval
      # abs(eigen(A)$values)
      
      
      seedB <- 11
      set.seed(seedB)
      B <- outer(rnorm(L), rnorm(L), "*")/3.5
      diag(B) <- NA
      B[which(abs(B) < quantile(abs(B)[!is.na(abs(B))], q))] <- 0
      diag(B) <- diagval
      # abs(eigen(B)$values)
      
      theta_blin <- k*kronecker(B, diag(S)) + k*kronecker(diag(L), t(A))
      if(max(abs(eigen(theta_blin)$values)) < 1){
        vecSigma_blin <- solve(diag(S^2*L^2) - kronecker(theta_blin, theta_blin), c(diag(S*L)))
        s1blin <- sum(diag(crossprod(theta_blin) %*% matrix(vecSigma_blin, S*L, S*L)))
        
        temp[i,2] <- R2blin <- s1blin / (s1blin + S*L)
      } else {
        temp[i,2] <- NA
      }
      
      theta_bilinear <- kronecker(B, A)*k^2
      if(max(abs(eigen(theta_bilinear)$values)) < 1){
        
        vecSigma_bilinear <- solve(diag(S^2*L^2) - kronecker(theta_bilinear, theta_bilinear), c(diag(S*L)))
        s1bilinear <- sum(diag(crossprod(theta_bilinear) %*% matrix(vecSigma_bilinear, S*L, S*L)))
        temp[i,3] <- R2bilinear <- s1bilinear / (s1bilinear + S*L)
      } else {
        temp[i,3] <- NA
      }
      
      # cat("frac=", i/length(ks), "\n")
    }
    
    
    k1 <- ks[which.min(abs(temp[,2] - R2target))]
    Ablin <- k1*A
    Bblin <- k1*B
    
    k2 <- ks[which.min(abs(temp[,3] - R2target))]
    Abilinear <- k2*A
    Bbilinear <- k2*B
    
    colnames(temp) <- c("k", "R2_blin", "R2_bilinear")
    write.table(temp, file.path(outdir, paste0("ks_q", 100*q, ".txt")))
    # print(temp)
    
    
    R2s <- array(0, c(nsims, 4, length(tmaxes), 2))
    Ahats <- array(0, c(S, S, nsims, 4, length(tmaxes), 2, nfold))
    Bhats <- array(0, c(L, L, nsims, 4, length(tmaxes), 2, nfold))
    dimnames(R2s)[[4]] <- dimnames(Ahats)[[6]] <- dimnames(Bhats)[[6]] <- gen_types
    
    
    
    
    for(gen_type in gen_types){
      
      for(i in 1:nsims){
        set.seed(i)
        Y0 <- matrix(rnorm(S*L), S, L)*tau
        Yall <- array(0, c(S, L, max(tmaxes) + 1 + burn))
        
        if(gen_type == "bilinear"){
          A <- Abilinear
          B <- Bbilinear
          Yall[,,1] <- t(A) %*% Y0 %*% B + tau*matrix(rnorm(S*L), S, L)
          for(j in 1 + 1:(burn + max(tmaxes)) ){
            Yall[,,j] <- t(A) %*% Yall[,,j-1] %*% B + tau*matrix(rnorm(S*L), S, L)
          }
        } else if (gen_type == "blin"){
          A <- Ablin
          B <- Bblin
          Yall[,,1] <- t(A) %*% Y0 + Y0 %*% B + tau*matrix(rnorm(S*L), S, L)
          for(j in 1 + 1:(burn + max(tmaxes)) ){
            Yall[,,j] <- t(A) %*% Yall[,,j-1] + Yall[,,j-1] %*% B + tau*matrix(rnorm(S*L), S, L)
          }
        }
        
        for(k in 1:length(tmaxes)){
          tmax <- tmaxes[k]
          
          Y <- Yall[,,burn + 2:(tmax + 1)]
          D <- Yall[,,burn + 1:tmax]
          Yp <- array(0, c(dim(Y), 4))
          if(tmax > nfold){
            set.seed(1)
            icv <- sample((1:tmax) %% nfold + 1)
          } else {
            icv <- 1:nfold
          }
          
          
          for(j in 1:nfold){
            itest <- which(icv == j)
            itrain <- which(icv != j)
            
            Ytrain <- Y[,,itrain]
            Dtrain <- D[,,itrain]
            Dtest <- D[,,itest]
            
            res <- mlm.ALS(Ytrain, Dtrain, verbose=F, tol = 1e-10, imax=1000)
            Yp[,,itest,1] <- drop(tprod(Dtest, res[[1]]))
            Ahats[,,i,1,k,gen_type,j] <- t(res[[1]][[1]])
            Bhats[,,i,1,k,gen_type,j] <- t(res[[1]][[2]])
            
            r <- fit_MLE_array_additive(Ytrain, Dtrain, X =NA, verbose=FALSE, use_cov=FALSE, tol = 1e-8)
            Yp[,,itest,2] <- amprod(Dtest, r$A, 1) + amprod(Dtest, t(r$B), 2)
            Ahats[,,i,2,k,gen_type,j] <- r$A
            Bhats[,,i,2,k,gen_type,j] <- r$B
            
            r <- additive_regression(Ytrain, Dtrain, NA, test_rank = F, type="biten", use_cov=FALSE, penalty = 1, whichlambda="min")
            Yp[,,itest,3] <- amprod(Dtest, r$A, 1) + amprod(Dtest, t(r$B), 2)
            Ahats[,,i,3,k,gen_type,j] <- r$A
            Bhats[,,i,3,k,gen_type,j] <- r$B
            
            r <- fit_MLE_array(Ytrain, Dtrain, NA, 1, 1, verbose=F, init="ran", sigma_init=1, use_cov=FALSE)
            Yp[,,itest,4] <- amprod(Dtest, r$A, 1) + amprod(Dtest, t(r$B), 2)
            Ahats[,,i,4,k,gen_type,j] <- r$A
            Bhats[,,i,4,k,gen_type,j] <- r$B
          }
          
          y0 <- sum(Y^2)
          R2s[i,,k, gen_type] <- sapply(1:4, function(z) 1 - sum((Yp[,,,z] - Y)^2) / y0)
        }
        
        cat("sim ", i, "\n")
        if(as.logical(writeout)){
          save(R2s, file=file.path(outdir, paste0("R2_q", q*100, ".rda")))
        }
        if(as.logical(writeout)){
          save(Ahats, Bhats, file=file.path(outdir_coef, paste0("coefs_q", q*100, ".rda")))
        }
      }
      
    }
    
    
    apply(R2s[,,,1,drop=FALSE], 2:3, mean)
    apply(R2s[,,,2,drop=FALSE], 2:3, mean)
    
    apply(R2s[,,,1,drop=FALSE], 2:3, sd)
    apply(R2s[,,,2,drop=FALSE], 2:3, sd)
    
    
    
  }
  
  
  
}
####


#### Plots for simulation study 
plot_cv <- function(qs=c(0,.5,.9),                         # fractions of zero entries
                    plotdir=file.path(getwd(), "plots"),    # where to write
                    readdir=file.path(getwd(), "results")) # where to read
                    
{
  
  require("RColorBrewer")
  dir.create(plotdir, showWarnings = F)
  colors <- c(brewer.pal(11, "Spectral")[c(1,8:10)], "gray60")
  fittypes <- c("Bilinear", "BLIN", "BLIN Sparse", "BLIN Reduced" )
  tmaxes <- c(10, 20, 50)
  
  setwd("readdir")
  for(q in qs){
    
    ps <- paste0("_q", q*100)
    load(paste0("R2", ps, ".rda"))
    
    dimnames(R2s)[[2]] <- fittypes   #<- c("Bilinear", "BLIN", "BLIN Sparse", "BLIN Reduced" )
    dimnames(R2s)[[3]] <- tmaxes    # <- c(10, 20, 50)
    gentypes <- dimnames(R2s)[[4]]
    
    
    blw <- 2.5
    count <- 0
    for(g in gentypes){
      count <- count + 1
      
      ncol <- length(tmaxes)
      boxtoplot <- NULL
      for(i in 1:ncol){
        temp <- apply(R2s[,,i,g], 2, quantile, c(.00, .1, .5, .9, 1), na.rm=TRUE)
        boxtoplot <- cbind(boxtoplot, temp)
        if(i < ncol){
          boxtoplot <- cbind(boxtoplot, NA)
        }
      }
      
      
      png(file.path(plotdir, paste0("box_R2out_gen_", g, ps,".png")), width=7, height=5, res=200, units="in")
      par(mar = c(3.2,3.2,.5,3.2), mgp=c(2.2,.6,0), cex.lab=1.5, cex.axis=1.5) 
      boxplot(boxtoplot,  
              xlab="T", 
              ylim=c(-1,1),
              xaxt="n",
              border=colors,
              range=0,
              boxlwd=blw, whisklty=1, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw, medlwd=blw)
      n2 <- length(fittypes)+1
      axis(1, tmaxes, at=seq(n2/2, ncol*n2, by=n2))                              
      abline(h=c(0, .75), lty=2, lwd=2, col="gray60")
      
      if(count == 1 & q %in% c(0.50, 0.90)){
        legend("bottomright", paste0(fittypes), lwd=blw, col=colors, lt=1,  bty="n", cex=1.55, title="Estimated Model")
      }
      
      dev.off()
    }
  }
  
  
}
####





