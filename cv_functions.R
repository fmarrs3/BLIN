

#### Run simulation study
run_first_cv_sims <- function(nsims=1000, 
                              verbose=F,
                              outdir=file.path(getwd(), "results"))  # where to write
{
  use_cov <- gen_cov <- F
  k <- m <- 2   # ranks of BLIN fits
  S <- 10
  L <- 9
  tau <- 50
  dir.create(outdir, showWarnings = F)
  
    for(tmax in c(5,10,20)){
    
      for(gen_type in c("biten_full", "bilinear", "biten")){
        
        # Make error file
        if(use_cov){ post <- "_withcov" } else {post <- ""}
        error_file <- paste0("errors_gen_", gen_type, "_S", S, "_tau", tau, "_t", tmax, post, ".txt")
        out <- data.frame(matrix(NA, nrow=0, ncol=29))
        names(out) <- c("gen_type", "fit_type", "S", "L", "tmax", "tau", "sim", "seed", "rankA", "rankB", "init", "CV", "rmse_out", "made_out", "rmse_in", "made_in", "R2_out", "R2_in", 
                        "fittime", "withcov", "obs_sparsityA", "obs_sparsityB", "precA", "recallA", "precB", "recallB", "MSEA", "MSEB", "MSEdiag")
        write.table(out, file.path(outdir, error_file), quote=F, sep = "\t", row.names = F)
        for(sim in 1:nsims){
          
          # Generate data
          seed <- sim
          data <- generate_data(S, L, tmax, tau, 
                                rankA <- S*(gen_type != "biten") + 2*(gen_type == "biten"), 
                                rankB <- L*(gen_type != "biten") + 2*(gen_type == "biten"), 
                                muAB=0, sigmaAB=1, sigmaD=1, genAB=T, use_cov=gen_cov, seed=seed, gen_type=gen_type)
          Atrue <- data$A
          Btrue <- data$B
          Y <- data$Y
          D <- data$D
          X <- data$X
          
          for(type in c("biten", "bilinear", "biten_full", "biten_sparse")){
            results <- misspec_model_fit(type, Y, X, D, k, m, use_cov, gen_type, S, L, tmax, tau, sim, seed, rankA, rankB, Atrue, Btrue)
            
            write.table(results$out, file.path(outdir, error_file), quote=F, sep = "\t", row.names = F, col.names=F, append=T)
          }
          
        }
        if(verbose){
          cat("tmax: \t", tmax, "\t data generation type: \t", gen_type, "\n")
        }
        
        
        
        if(gen_type == "biten_full"){   # sparse version
          # Make error file
          if(use_cov){ post <- "_withcov" } else {post <- ""}
          error_file <- paste0("errors_gen_sparse_", gen_type, "_S", S, "_tau", tau, "_t", tmax, post, ".txt")
          out <- data.frame(matrix(NA, nrow=0, ncol=29))
          names(out) <- c("gen_type", "fit_type", "S", "L", "tmax", "tau", "sim", "seed", "rankA", "rankB", "init", "CV", "rmse_out", "made_out", "rmse_in", "made_in", "R2_out", "R2_in", 
                          "fittime", "withcov", "obs_sparsityA", "obs_sparsityB", "precA", "recallA", "precB", "recallB", "MSEA", "MSEB", "MSEdiag")
          write.table(out, file.path(outdir, error_file), quote=F, sep = "\t", row.names = F)
          for(sim in 1:nsims){
            
            # Generate data
            seed <- sim
            data <- generate_data(S, L, tmax, tau, 
                                  rankA <- S*(gen_type != "biten") + 2*(gen_type == "biten"), 
                                  rankB <- L*(gen_type != "biten") + 2*(gen_type == "biten"), 
                                  muAB=0, sigmaAB=1, sigmaD=1, genAB=T, use_cov=gen_cov, seed=seed, gen_type=gen_type, sparse=2/S)
            Atrue <- data$A
            Btrue <- data$B
            Y <- data$Y
            D <- data$D
            X <- data$X
            
            for(type in c("biten", "bilinear", "biten_full", "biten_sparse")){
              results <- misspec_model_fit(type, Y, X, D, k, m, use_cov, gen_type, S, L, tmax, tau, sim, seed, rankA, rankB, Atrue, Btrue)
              
              write.table(results$out, file.path(outdir, error_file), quote=F, sep = "\t", row.names = F, col.names=F, append=T)
            }
            
          }
          
          
          if(verbose){
            cat("tmax: \t", tmax, "\t data generation type: \t", gen_type, "sparse \n")
          }
          
        }
        
        
        
      }
    }
  
}
####


#### Plots for simulation study
plot_first_cv <- function(plotdir=file.path(getwd(), "plots"), 
                          readdir=file.path(getwd(), "results"))
{
  
  require("RColorBrewer")
  dir.create(plotdir, showWarnings = F)

  
  post <- ".txt" # "t10_withcov.txt"
  colors <- c(brewer.pal(8, "Set3")[c(1,6,5,4)], "gray60") #[2:4], "black", "gray60")
  colorsnew <- c(brewer.pal(11, "Spectral")[c(8:10,1)], NA)
  
  # inpath <- "./results_rerun1000_sparse_par"   # ./results_rerun1000_par
  for(inpath in readdir){
    
    
    # Read in data
    d <- NULL
    count <- 0
    files <- list.files(path = inpath, pattern=paste0("errors_+.*(",post,")+"))
    for(f in files){
      count <- count + 1
      data <- read.table(file.path(inpath, f), header=T)
      if(length(grep("sparse", f))   > 0){
        data$gen_type <- "biten_sparse"
      }
      
      # if(length(data$fittime) == 0){ data$fittime <- NA}  # missing column in OLS code
      if(count == 1){ d <- data} else {
        d <- rbind(d, data)
      }
    }
    
    # unique(d$gen_type)
    # d$gen_type[d$gen_type == 1] <- "bilinear"
    # d$gen_type[d$gen_type == 2] <- "biten_full"
    # d$gen_type[d$gen_type == 3] <- "biten"
    # 
    plot_types <- c("biten_full", "biten", "biten_sparse", "bilinear") 
      # c("bilinear", "biten_full", "biten", "biten_sparse")   #, "mu")  mean model at zero for both
    gen_types <- unique(d$gen_type)  #  c("bilinear", "biten_full", "biten")
    

    
    
    pull_data_Sfixed_time <- function(gen_type, tau, S, plot_types, plot_var, trange, qin, meanin=T)
    {
      box <- box_avg <- NULL
      
      if(min(trange) < 10){   # do quantiles only when tmax < 5
        box_avg_q <- NULL
        
        for(tmax in trange){
          # plot_var <- "rmse_out"
          for(type in plot_types){
            rows <- d$gen_type == gen_type & d$fit_type == type & d$tau == tau & d$S == S & d$tmax == tmax   # rows of interest
            # unavg <- d[rows, plot_var]
            unavg <- as.numeric(sapply(unique(d[rows, "sim"]), function(x) mean(d[rows, plot_var][d[rows, "sim"] == x])))
            unavg2 <- quantile(unavg, c(.025,.1,.5,.9,.975), na.rm=T)
            unavg2[3] <- mean(unavg, na.rm=T)
            box_avg_q <- cbind(box_avg_q, unavg2)   # quantiles
          }
          box_avg_q <- cbind(box_avg_q, matrix(NA, nrow=nrow(box_avg_q), ncol=1))   # columns of NAs for space between boxplots
        }
        
        box_avg_q <- box_avg_q[,-ncol(box_avg_q)]
        
        
        return(list(box=NA, box_avg=NA, box_avg_q=box_avg_q))
      } else {
        for(tmax in trange){
          # plot_var <- "rmse_out"
          for(type in plot_types){
            rows <- d$gen_type == gen_type & d$fit_type == type & d$tau == tau & d$S == S & d$tmax == tmax   # rows of interest
            unavg <- d[rows, plot_var]
            box <- cbind(box, as.numeric(unavg))   # not averaged over sims, i.e. not avg CV
            
            box_avg <- cbind(box_avg, as.numeric(sapply(unique(d[rows, "sim"]), function(x) mean(d[rows, plot_var][d[rows, "sim"] == x]))))  # average over simulation numbers
          }
          box <- cbind(box, matrix(NA, nrow=nrow(box), ncol=1))   # columns of NAs for space between boxplots
          box_avg <- cbind(box_avg, matrix(NA, nrow=nrow(box_avg), ncol=1))   # columns of NAs for space between boxplots
        }
        
        box <- box[,-ncol(box)]   # remove last column
        box_avg <- box_avg[,-ncol(box_avg)]
        
        # Quantiles of average box
        box_avg_q <- apply(box_avg, 2, quantile, qin, na.rm=T)
        if(meanin){
          box_avg_q[3,] <- apply(box_avg, 2, mean, na.rm=T)
        }
        
        return(list(box=box, box_avg=box_avg, box_avg_q=box_avg_q))
      }
    }
    
    
    
    
    # Time on x-axis, quantiles
    taurange <- 50
    trange <- c(5,10,20)
    blw <- 2.5
    qin <- c(0,.1,.5,.9,1)
    
    for(S in c(10)){
      for(tau in taurange){
        for(gen_type in gen_types){
          
          # Plot MSPE (out of sample)
          plot_var <- "R2_out"
          boxout <- pull_data_Sfixed_time(gen_type, tau, S, plot_types, plot_var, trange, qin, meanin=F)$box_avg_q
          plot_var <- "R2_in"
          boxin <- pull_data_Sfixed_time(gen_type, tau, S, plot_types, plot_var, trange, meanin=F)$box_avg_q
          
          yrange <- c(-1,1)
          if(gen_type == "bilinear"){ 
            yax = "s"
            margs <- c(3,2,.5,.5)
          } else {
            yax = "n"
            margs <- c(3,.5,.5,.5)
          }
          
          pdf(file.path(plotdir, paste0("line_gen_", gen_type, "_S", S,  "_tau", tau, ".pdf")), width=4, height=4.4)
          par(mar=margs, cex.axis=1.3, cex.lab=1.3, mgp=c(1.7,.5,0))
          plot(NA, NA, xlab=expression(T), xaxt="n",ylim=yrange, yaxt=yax, xlim=c(1,ncol(boxout)), ylab="")
          jitter = .17
          for(j in 1:4){
            # j <- 1
            xseq <- 1:3
            yseq <- (xseq-1)*5 + j
            boxoutplot <- boxout[3,yseq]
            boxoutplot[boxoutplot < -1] <- NA
            points(yseq, boxoutplot, pch=1, type="o", lwd=1.6, col=colorsnew[j], cex=1.2)
            boxarrows <- boxout[, yseq]
            # boxarrows[, is.na(boxoutplot)] <- NA
            arrows(x0=yseq, x1=yseq, y0=boxarrows[1,], y1=boxarrows[5,], code=3, angle=90, length=.06, lwd=1.6, col=colorsnew[j])
            
            points(yseq+jitter, boxin[3,yseq], pch=4, type="o", lwd=1.6, col=colorsnew[j], cex=1.2, lty=2)
            arrows(x0=yseq+jitter, x1=yseq+jitter, y0=boxin[1,yseq], y1=boxin[5,yseq], code=3, angle=90, length=0.06, lwd=1.6, col=colorsnew[j], lty=1)
          }
          
          
          atrange <- seq((length(plot_types)+1)/2, length(trange)*(length(plot_types)+1), by=(length(plot_types)+1))
          axis(1, trange, at=atrange)
          abline(v=atrange[-length(atrange)] + length(plot_types)/2 + 1/2, lty=1, lwd=.7, col="gray80")
          abline(h=0, col="black", lwd=2, lty=3)
          if(gen_type == "biten"){
            legend("bottomleft", legend=c("Bilinear", "BLIN", "BLIN, rank 2", "BLIN, sparse",  'in-sample'), col=c(colorsnew[c(4,1:3)], "gray50"), lwd=.9*1.6, bty = "n", cex = 1.2,
                   lty=c(rep(1,4),2), pch=c(rep(1,4),4) )
          }
          dev.off()
          
          ### Old plots
          # # Plot MSPE (out of sample)
          # plot_var <- "R2_out"
          # box <- pull_data_Sfixed_time(gen_type, tau, S, plot_types, plot_var, trange, qin, meanin=F)$box_avg_q
          # # if(gen_type == "bilinear"){
          # #   yrange <- range(box, na.rm=T)
          # # } else { yrange <- c(0,1)}
          # yrange <- c(-1,1)
          # if(gen_type == "bilinear"){ yax = "s"} else {yax = "n"}
          # 
          # pdf(file.path(plotdir, paste0("box_R2_out_gen_", gen_type, "_S", S,  "_tau", tau, ".pdf")), width=4, height=4)
          # par(mar=c(3,2,.5,2), cex.axis=1.3, cex.lab=1.3, mgp=c(1.7,.5,0))
          # boxplot(box, border=colors, 
          #         # ylab=expression("CV-average Out-of-sample R"^"2"), 
          #         xlab=expression(T), 
          #         xaxt="n", 
          #         ylim=yrange, range=0, whisklty=1, boxlwd=blw, stapleltwd=blw, whisklwd=blw, medlwd=blw, staplelwd=blw,
          #         yaxt=yax)
          # atrange <- seq((length(plot_types)+1)/2, length(trange)*(length(plot_types)+1), by=(length(plot_types)+1))
          # axis(1, trange, at=atrange)
          # abline(v=atrange[-length(atrange)] + length(plot_types)/2 + 1/2, lty=1, lwd=.7, col="gray80")
          # abline(h=0, col="black", lwd=2, lty=3)
          # # if(gen_type == "bilinear"){
          # #   legend("topleft", legend=c("Bilinear", "BiTEN", "BiTEN, rank 2", "BiTEN, sparse"), col=colors, lwd=.9*blw, bty = "n", cex = .9)
          # # }
          # dev.off()
          # 
          # 
          # # Plot MSPE in-sample
          # plot_var <- "R2_in"
          # box <- pull_data_Sfixed_time(gen_type, tau, S, plot_types, plot_var, trange, meanin=F)$box_avg_q
          # # if(gen_type == "bilinear"){
          # #   yrange <- range(box, na.rm=T)
          # # } else { yrange <- c(0,1)}
          # 
          # pdf(file.path(plotdir, paste0("box_R2_in_gen_", gen_type, "_S", S,  "_tau", tau, ".pdf")), width=4, height=3.5)
          # par(mar=c(.5,2,.5,2), cex.axis=1.3, cex.lab=1.3, mgp=c(1.7,.5,0))
          # boxplot(box, border=colors, 
          #         # ylab=expression("CV-average In-sample R"^"2"), 
          #         xlab=expression(T), xaxt="n", ylim=c(-1,1), 
          #         range=0, whisklty=1, boxlwd=blw, stapleltwd=blw, whisklwd=blw, medlwd=blw, staplelwd=blw,
          #         yaxt=yax)
          # atrange <- seq((length(plot_types)+1)/2, length(trange)*(length(plot_types)+1), by=(length(plot_types)+1))
          # # axis(1, trange, at=atrange)
          # abline(v=atrange[-length(atrange)] + length(plot_types)/2 + 1/2, lty=1, lwd=.7, col="gray80")
          # abline(h=0, col="black", lwd=2, lty=3)
          # if(gen_type == "bilinear"){
          #   legend("bottomleft", legend=c("Bilinear", "BiTEN", "BiTEN, rank 2", "BiTEN, sparse"), col=colors, lwd=1.3*blw, bty = "n", cex = 1.3)
          # }
          # dev.off()
          
        }
      }}
    
  }
}
####





