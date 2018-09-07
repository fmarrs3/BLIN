

#### Fit BLIN model to SID data
fit_blin_sid <- function(outdir=file.path(getwd(), "results"),  # where to write
                         readdir=file.path(getwd()),      # where to read data from
                         verbose=F # write out into window?
)
{
  
  dir.create(outdir)

  use_cov <- F
  lag <- 2
  diff <- T
  
  #### data
  data <- read_state_interaction_data(use_cov=F, agg=1, remakeVAR=T, lag=lag, diff=diff)
  Y <- data$Y  ; D <- data$D  ;  X <- data$X
  resp_type <- dimnames(Y)[[3]]
  ####
  

  #### Sparse version
  whichlambda <- "min"  # choose "min" or "1se"  lamabda from lasso fit
  P2 <- matrix(NA, 1, 4)


  #### cvfitting
  YP<-Y*0


  for(j in 1:4)
  {
    out <- additive_regression(Y[,,j,], D[,,j,], X[,,j,,], type="biten", use_cov=use_cov, test_rank=F, penalty=1, whichlambda=whichlambda)   # lasso fit
    if(use_cov){
      Xb <- drop(amprod(X[,,j,,], out$beta, 4))
      cov_name <- "_withcov"
    } else {
      Xb <- 0
      cov_name <- ""
    }
    YP[,,j,]<- amprod(D[,,j,], out$A, 1) + amprod(D[,,j,], t(out$B), 2) + Xb
    write.table(out$A, file=file.path(outdir, paste0("A_fullfit_biten_sparse_FM_", resp_type[j], "_lag", lag, "_diff", as.numeric(diff), "_L", whichlambda, cov_name,".txt")))
    write.table(out$B, file=file.path(outdir, paste0("B_fullfit_biten_sparse_FM_", resp_type[j], "_lag", lag, "_diff", as.numeric(diff), cov_name, "_L", whichlambda, ".txt")))
    write.table(out$beta, file=file.path(outdir, paste0("beta_fullfit_biten_sparse_FM_", resp_type[j], "_lag", lag, "_diff", as.numeric(diff), cov_name, "_L", whichlambda, ".txt")))

    #### evaluation
    P2[1,] <- 1-apply((YP-Y)^2,3,mean,na.rm=TRUE)/
      apply(Y^2,3,mean,na.rm=TRUE)
    ####

    write.table(P2, file=file.path(outdir, paste0("P2_fullfit_biten_sparse_FM", "_lag", lag, "_diff", as.numeric(diff),"_L", whichlambda,  cov_name,".txt")))
    if(verbose){
      cat("Response ", j, "\n")
    }
  }
  # cat(cv,"\n")
  ####

  if(verbose){
    cat("\nDone with sparse BLIN fit\n")
  }
  ####
  
  
  
}
####


#### Plot results from SID fits
plot_blin_sid <- function(outdir=file.path(getwd(), "plots"),  # where to write
                          readdir=file.path(getwd(), "results"), 
                          countriesdir=getwd(), 
                          fittype="biten_sparse")
{
  plotdir <- outdir
  dir.create(outdir)
  
  
  diff <- TRUE
  lag <- 2
  resp <- c("mn", "mp", "vn", "vp")
    
  spacings <- cbind(c(-8,-9.5,-8,-8), c(-9,-8.5,-9,-8)) - 1.5
  
  rownames(spacings)[1:length(resp)] <- resp
  horizontal <- T
  
  # colors <- brewer.pal(8, "Set2")[-6]
  colors <- brewer.pal(11, "Spectral")[c(1,3,5,7,9,11)]
  codes <- read.table(file.path(countriesdir, "countries.txt"), header=T, stringsAsFactors = F)
  resp <-  c("mn", "mp", "vn", "vp")
  countries <- as.vector(codes$country)
  
  for(r in resp){
    A <- read.table(file.path(readdir, paste0("A_fullfit_", fittype, "_FM_", r, "_lag", lag,"_diff", as.numeric(diff),"_Lmin.txt")))
    B <- read.table(file.path(readdir, paste0("B_fullfit_", fittype, "_FM_", r, "_lag", lag, "_diff", as.numeric(diff), "_Lmin.txt")))
    
    A <- as.matrix(A)
    B <- as.matrix(B)
    diag(A) <- diag(B) <- NA
    colnames(A) <- colnames(B) <- rownames(A) <- rownames(B) <- countries
    
    
    maxval <- max(abs(c(c(A), c(B))), na.rm=T)
    maxdeg <- max( c(rowSums(A!=0, na.rm=T) + colSums(A!=0, na.rm=T),  rowSums(B!=0, na.rm=T) + colSums(B!=0, na.rm=T)) )
    
    mplot <- A*(A>0)
    groups <- codes$continent
    
    pdf(file.path(plotdir, paste0("Aall_", r,"_arc.pdf")), width=6, height=4)
    par(mar=c(0,0,0,0))  #, cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    arcplot_fm2(A, groups, colors, linespace = spacings[r, 1], maxwt=maxval, maxdeg=maxdeg, horizontal=horizontal)
    dev.off()
    
    pdf(file.path(plotdir, paste0("Ball_", r,"_arc.pdf")), width=6, height=4)
    par(mar=c(0,0,0,0))  #, cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
    arcplot_fm2(B, groups, colors, linespace = spacings[r, 2], maxwt=maxval, maxdeg=maxdeg, horizontal=horizontal)
    dev.off()
    
  } 
  
  
  
  # data <- read_state_interaction_data(use_cov=F, agg=1, remakeVAR=T, lag=lag, diff=diff)
  # Y <- data$Y  ; D <- data$D  ;  X <- data$X
  # resp_type <- dimnames(Y)[[3]]
  # S <- dim(Y)[1]   ;   L <- dim(Y)[2]   ;   tmax <- dim(Y)[4]
  # 
  # r <-  c("mn", "mp", "vn", "vp")
  # fittype <- "biten_sparse"
  # 
  # arccolors <- brewer.pal(11, "Spectral")[c(1,3,5,7,9,11)]
  # icolors <- match(codes$continent, sort(unique(codes$continent)))
  # 
  # for(m in 1:4){
  #   # outdir <- paste0(instub, "_j", j)
  #   r <- resp_type[m]
  #   if(lag == 2){ps <- paste0("_lag", lag, "_diff", as.numeric(diff), "_Lmin")} else {ps <- ""}
  #   
  #   # A <- as.matrix( read.table(list.files(pattern=paste0("A_fullfit_", fittype, "_FM_", r, ps, ".txt"))) )
  #   # B <- as.matrix( read.table(list.files(pattern=paste0("B_fullfit_", fittype, "_FM_", r, ps, ".txt"))) )
  #   A <- as.matrix( read.table(file.path(readdir, paste0("A_fullfit_", fittype, "_FM_", r, "_lag", lag,"_diff", as.numeric(diff), "_Lmin.txt"))) )
  #   B <- as.matrix( read.table(file.path(readdir, paste0("B_fullfit_", fittype, "_FM_", r, "_lag", lag, "_diff", as.numeric(diff), "_Lmin.txt"))) )
  #   
  #   
  #   Yhat <- amprod(D[,,m,], A, 1) + amprod(D[,,m,], t(B), 2)
  #   e <- Y[,,m,] - Yhat 
  #   
  #   
  #   # Matricizing resids
  #   mm <- mat(e,1)
  #   mm[is.na(mm)] <- 0
  #   ee <- tcrossprod(mm)
  #   pplot <- prcomp(ee, scale=T)$rotation
  #   pdf(file.path(plotdir,paste0("exporter_resid_PCA_", resp_type[m], ".pdf")), width=4, height=4)
  #   par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  #   plot(NA, NA, xlim=range(pplot[,1]), ylim=range(pplot[,2]), xlab="PCA 1", ylab="PCA 2")
  #   text(pplot[,1], pplot[,2], labels=rownames(pplot), offset=0, col=arccolors[icolors], cex=1.4) #rownames(pplot)) # text(rownames(pplot)))
  #   dev.off()
  #   
  #   
  #   mm <- mat(e,2)
  #   mm[is.na(mm)] <- 0
  #   ee <- tcrossprod(mm)
  #   pplot <- prcomp(ee, scale=T)$rotation
  #   pdf(file.path(plotdir,paste0("importer_resid_PCA_", resp_type[m], ".pdf")), width=4, height=4)
  #   par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  #   plot(NA, NA, xlim=range(pplot[,1]), ylim=range(pplot[,2]), xlab="PCA 1", ylab="PCA 2")
  #   text(pplot[,1], pplot[,2], labels=rownames(pplot), col=arccolors[icolors], offset=0, cex=1.4) #rownames(pplot)) # text(rownames(pplot)))
  #   dev.off()
  #   
  #   
  #   mm <- mat(e,3)
  #   mm[is.na(mm)] <- 0
  #   ee <- tcrossprod(mm)
  #   eecor <- cov2cor(ee)
  #   rows <- matrix(1:tmax, tmax, tmax)
  #   cols <- matrix(1:tmax, tmax, tmax, byrow = TRUE)
  #   imat <- abs(rows - cols)
  #   ncor <- 10
  #   blw <- 2
  #   quants <- c(.025,.1,.5,.9,.975)
  #   boxtoplot <- matrix(0, 5, ncor)
  #   for(j in 1:ncor){
  #     boxtoplot[,j] <- quantile(eecor[imat == j], quants, na.rm=TRUE)
  #   }
  #   
  #   pdf(file.path(plotdir,paste0("year_resid_boxplot_", resp_type[m], ".pdf")), width=4, height=4)
  #   par(mar=c(3,3,.5,.5), cex.axis=.8, cex.lab=.8, mgp=c(1.7,.5,0))
  #   boxplot(boxtoplot, xlab="lag", ylab="Estimated correlation", range=0, border=arccolors[1],
  #           boxlwd=blw, whisklty=1, whisklwd=blw, staplelty=1, staplelwd=blw, outlwd=blw)
  #   abline(h=0, col="gray80")
  #   dev.off()
  #   
  # }
  
}
####

