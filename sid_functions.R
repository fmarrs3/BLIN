

#### Fit BLIN model to SID data
fit_blin_sid <- function(outdir=file.path(getwd(), "results"),  # where to write
                         readdir=file.path(getwd()),      # where to read 
                         verbose=FALSE) # write out into window?
{
  diff <- TRUE
  use_cov <- FALSE
  ps <- "_multi2"
  lagrange <- 5
  ntail <- 543 -  lagrange
  whichlambda <- "min"
  
  P2 <- kb <- AICs <- array(NA, c(2, lagrange, lagrange, lagrange))
  if(as.logical(use_cov)){ stop("Only coded for without covariates")}
  dir.create(outdir)

  lagvec <- c(5,3,1)   # lag for each mode
      
  #### data
  data <- read_state_interaction_data_multi(lagvec=lagvec, diff=diff, readdir=getwd())
  Y <- data$Y  ; D <- data$D  ;  X <- data$X
  ####
  
  #### fitting
  YP<-Y*0
  out <- additive_regression_multi(Y, D, X, whichlambda=whichlambda)   # lasso fit

  # Write out coefficients
  for(j in 1:length(out$B)){
    write.table(out$B[[j]], file=file.path(outdir, paste0("B", j, "_fullfit_biten_sparse_FM_", "_lag", paste0(c(lagvec), collapse=""), "_diff", as.numeric(diff),  ps, ".txt")))
  }
  ####

}
####
  
  


#### Plot results from SID fits
plot_blin_sid <- function(outdir=file.path(getwd(), "plots"),  # where to write
                          readdir=file.path(getwd(), "results"),    # where to read
                          countriesdir=getwd(),    # where to read countries from
                          fittype="biten_sparse")   # 
{
  
  plotdir <- outdir
  dir.create(plotdir)
  require(igraph)
  require(RColorBrewer)
  resp <- c("mn", "mp", "vn", "vp")
  diff <- TRUE
  lagvec <- c(5,3,1)
  ps <- "_multi2"
  
  # Read in networks and country abbreviations
  codes <- read.table(file.path(countriesdir, "countries.txt"), header=T, stringsAsFactors = F)
  countries <- as.vector(codes$country)
  
  

  j <- 1
  filein <- list.files(path=readdir, pattern=paste0("B", j, "_fullfit_biten_sparse_FM_", "_lag", paste0(c(lagvec), collapse=""), "_diff", as.numeric(diff),  ps, ".txt"))
  if(length(filein) == 1){
    B <- vector("list", 3)
    for(j in 1:length(B)){
      B[[j]] <- as.matrix( read.table(file=file.path(readdir, paste0("B", j, "_fullfit_biten_sparse_FM_", "_lag", paste0(c(lagvec), collapse=""), "_diff", as.numeric(diff),  ps, ".txt"))) )
    }
  }
  length(B)
  
  # outer(diag(B[[1]]), diag(B[[2]]), "+")
  
  
  
  lb <- -.0037
  ub <- .0121
  
  
  A <- B
  for(j in 1:length(B)){
    A[[j]] <- B[[j]]
    if(j < 3){
      colnames(A[[j]]) <- rownames(A[[j]]) <- countries
    } else {
      colnames(A[[j]]) <- rownames(A[[j]]) <- resp
    }
    # A[[j]][A[[j]] < ub & A[[j]] > lb] <- 0
    diag(A[[j]]) <- 0
  }
  
  
  
  A <- B
  for(j in 1:length(B)){
    A[[j]] <- B[[j]]
    if(j < 3){
      colnames(A[[j]]) <- rownames(A[[j]]) <- countries
    } else {
      colnames(A[[j]]) <- rownames(A[[j]]) <- resp
    }
    # A[[j]][A[[j]] < ub & A[[j]] > lb] <- 0
    diag(A[[j]]) <- 0
  }
  
  # Plot params
  frac <- 1.5   # dialation of C relative to A and B
  vlc <- .9    # vertex label cex
  av <- 6.5     # vertex size intercept
  bv <- 165     # vertex size slope
  eaw <- 2.    # edge arrow width
  ae <- 6       # edge size intercept
  be <- 200     # edge size slsope
  ecurve <- -.15
  vcol="gray25"
  alpha_node <- .45
  vdist <- 1.42
  
  continents <- unique( codes$continent[order(codes$continent)] )
  circle_colors <- brewer.pal(11, "Spectral")[c(1,5,11,3,9,7)]
  # circle_colors <- brewer.pal(6, "Spectral")
  # circle_colors_trans <- paste0(circle_colors, 95)
  library("igraph")
  
  j <- 1
  C <- as.matrix(A[[j]])
  colnames(C) <- rownames(C) <- countries
  i1 <- which(C != 0, arr.ind=TRUE)
  EL <- cbind(countries[i1[,1]], countries[i1[,2]])
  
  net <- graph_from_edgelist(EL, directed=TRUE)
  E(net)$weight <- C[i1]                           
  V(net)$label <- names(V(net))
  V(net)$continent <- codes$continent[match(names(V(net)), codes$country)]
  V(net)$color <- alpha(circle_colors[ match(V(net)$continent, continents) ], alpha_node)
  ordering <- V(net)$label[ order(V(net)$continent) ] #[keep]
  L <- layout_in_circle(net, order=ordering)
  # V(net)$label <- codes$country
  # V(net)$color <- circle_colors_trans[match(codes$continent[ match(V(net)$label, codes$country) ], continents)]
  
  # Colors
  temp0 <- circle_colors[ match(codes$continent[match(as_edgelist(net)[,1], countries)], continents) ]
  temp <- cbind(temp0, abs(E(net)$weight) / max(abs( E(net)$weight )) )
  ismall <- which(E(net)$weight <= ub & E(net)$weight >= lb)
  ibig <- setdiff(1:length(E(net)$weight), ismall)
  temp[ismall,1] <- "gray85"
  # E(net)$weight[ismall] <- 0
  temp2 <- cbind(temp[,1], abs(E(net)$weight) / max(abs( E(net)$weight )) )
  E(net)$color <- apply(temp2, 1, function(z) alpha(z[1], as.numeric(z[2])))
  
  E(net)$lty <- 3 - 2*(E(net)$weight >0)
  
  V(net)$size <- av + bv*(rowSums(abs(C)) )[V(net)$label]
  E(net)$width <- ae + be*abs(E(net)$weight)
  
  
  png(file.path(plotdir, paste0("B_circle_", j, "_v3.png")), width=5, height=5, res=200, units="in")
  par(mar = c(0, 0, 0, 0), mgp=c(1.7,.5,0), cex.lab=1.1, cex.axis=1.1) 
  plot(net, 
       layout=L, 
       edge.curved=ecurve,
       edge.arrow.width	=eaw,
       vertex.label.font=2,
       vertex.label.family="Helvetica",
       vertex.label.color=vcol, 
       vertex.label.cex=vlc, 
       vertex.label.degree=pi/2,
       vertex.label.dist=vdist,
       vertex.frame.color=NA)
  
  tnet <- net
  temp2[ismall, 2] <- 0
  temp2[ibig, 2] <- 1
  E(tnet)$color <- apply(temp2, 1, function(z) alpha(z[1], as.numeric(z[2])))
  
  plot(tnet, add=TRUE, 
       layout=L, 
       edge.curved=ecurve,
       edge.arrow.width	=eaw,
       vertex.label.font=2,
       vertex.label.family	= "Helvetica",
       vertex.label.color=vcol, 
       vertex.label.cex=vlc, 
       vertex.label.degree=pi/2,
       vertex.label.dist=vdist,
       vertex.frame.color=NA)
  
  dev.off()
  
  
  
  j <- 2
  C <- as.matrix(A[[j]])
  colnames(C) <- rownames(C) <- countries
  i1 <- which(C != 0, arr.ind=TRUE)
  EL <- cbind(countries[i1[,1]], countries[i1[,2]])
  
  net <- graph_from_edgelist(EL, directed=TRUE)
  E(net)$weight <- C[i1]                           
  V(net)$label <- names(V(net))
  V(net)$continent <- codes$continent[match(names(V(net)), codes$country)]
  V(net)$color <- alpha(circle_colors[ match(V(net)$continent, continents) ], alpha_node)
  ordering <- V(net)$label[ order(V(net)$continent) ] #[keep]
  L <- layout_in_circle(net, order=ordering)
  # V(net)$label <- codes$country
  # V(net)$color <- circle_colors_trans[match(codes$continent[ match(V(net)$label, codes$country) ], continents)]
  
  # Colors
  temp0 <- circle_colors[ match(codes$continent[match(as_edgelist(net)[,1], countries)], continents) ]
  temp <- cbind(temp0, abs(E(net)$weight) / max(abs( E(net)$weight )) )
  ismall <- which(E(net)$weight <= ub & E(net)$weight >= lb)
  ibig <- setdiff(1:length(E(net)$weight), ismall)
  temp[ismall,1] <- "gray85"
  # E(net)$weight[ismall] <- 0
  temp2 <- cbind(temp[,1], abs(E(net)$weight) / max(abs( E(net)$weight )) )
  E(net)$color <- apply(temp2, 1, function(z) alpha(z[1], as.numeric(z[2])))
  
  E(net)$lty <- 3 - 2*(E(net)$weight >0)
  
  V(net)$size <- (av + bv*(rowSums(abs(C))  )[V(net)$label])   #+ colSums(abs(C))
  E(net)$width <- (ae + be*abs(E(net)$weight))
  
  
  png(file.path(plotdir, paste0("B_circle_", j, "_v3.png")), width=5, height=5, res=200, units="in")
  par(mar = c(0, 0, 0, 0), mgp=c(1.7,.5,0), cex.lab=1.1, cex.axis=1.1) 
  plot(net, 
       layout=L, 
       edge.curved=ecurve,
       edge.arrow.width	=eaw,
       vertex.label.font=2,
       vertex.label.family	= "Helvetica",
       vertex.label.color=vcol, 
       vertex.label.cex=vlc, 
       vertex.label.degree=pi/2,
       vertex.label.dist=vdist,
       vertex.frame.color=NA)
  
  tnet <- net
  temp2[ismall, 2] <- 0
  temp2[ibig, 2] <- 1
  E(tnet)$color <- apply(temp2, 1, function(z) alpha(z[1], as.numeric(z[2])))
  
  plot(tnet, add=TRUE, 
       layout=L, 
       edge.curved=ecurve,
       edge.arrow.width	=eaw,
       vertex.label.font=2,
       vertex.label.family	= "Helvetica",
       vertex.label.color=vcol, 
       vertex.label.cex=vlc, 
       vertex.label.degree=pi/2,
       vertex.label.dist=vdist,
       vertex.frame.color=NA)
  
  dev.off()
  # edge.lty=E(net)$edge.lty, 
  
  # edge.width=E(net)$edge.width,
  # vertex.color=V(net)$vertex.color, 
  
  
  # vertex.size = V(net)$vertex.size, 
  
  
  
  
  circle_colors <- brewer.pal(11, "Spectral")[c(5,1,9,11)]  
  circle_colors_trans <- paste0(circle_colors, 95)
  
  j <- 3
  C <- as.matrix(A[[j]])
  colnames(C) <- rownames(C) <- resp
  i1 <- which(C != 0, arr.ind=TRUE)
  EL <- cbind(resp[i1[,1]], resp[i1[,2]])
  
  
  net <- graph_from_edgelist(EL, directed=TRUE)
  E(net)$weight <- C[i1]                           
  V(net)$label <- names(V(net))
  # V(net)$continent <- codes$continent[match(names(V(net)), codes$country)]
  V(net)$color <- alpha(circle_colors, alpha_node)
  ordering <- order(V(net)$label)
  L <- layout_in_circle(net, order=c("vn", "mn", "vp", "mp"))
  # V(net)$label <- codes$country
  # V(net)$color <- circle_colors_trans[match(codes$continent[ match(V(net)$label, codes$country) ], continents)]
  
  # Colors
  temp0 <- circle_colors[match(as_edgelist(net)[,1], V(net)$label)]
  temp <- cbind(temp0, abs(E(net)$weight) / max(abs( E(net)$weight )) )
  ismall <- which(E(net)$weight <= ub & E(net)$weight >= lb)
  ibig <- setdiff(1:length(E(net)$weight), ismall)
  temp[ismall,1] <- "gray85"
  # E(net)$weight[ismall] <- 0
  temp2 <- cbind(temp[,1], abs(E(net)$weight) / max(abs( E(net)$weight )) )
  E(net)$color <- apply(temp2, 1, function(z) alpha(z[1], as.numeric(z[2])))
  
  
  
  
  E(net)$lty <- 3 - 2*(E(net)$weight >0)
  
  V(net)$size <- (av + bv*(rowSums(abs(C))  )[V(net)$label])*frac   # + colSums(abs(C))
  E(net)$width <- (ae + be*abs(E(net)$weight))*frac
  
  
  png(file.path(plotdir, paste0("B_circle_", j, "_v3.png")), width=5/frac, height=5/frac, res=200, units="in")
  par(mar = c(1, 1, 1, 1), mgp=c(1.7,.5,0), cex.lab=1.1, cex.axis=1.1) 
  plot(net, 
       layout=L, 
       edge.curved=ecurve,
       edge.arrow.width	=eaw*frac,
       vertex.label.font=2,
       vertex.label.family	= "Helvetica",
       vertex.label.color=vcol, 
       vertex.label.cex=vlc*frac+.3, 
       vertex.label.degree=pi/2,
       vertex.label.dist=vdist*frac+c(0, 2.5, 0, 1),
       vertex.frame.color=NA)
  
  tnet <- net
  temp2[ismall, 2] <- 0
  temp2[ibig, 2] <- 1
  E(tnet)$color <- apply(temp2, 1, function(z) alpha(z[1], as.numeric(z[2])))
  
  
  plot(tnet, add=TRUE, 
       layout=L, 
       edge.curved=ecurve,
       edge.arrow.width	=eaw*frac,
       vertex.label.font=2,
       vertex.label.family	= "Helvetica",
       vertex.label.color=vcol, 
       vertex.label.cex=vlc*frac+.3, 
       vertex.label.degree=pi/2,
       vertex.label.dist=vdist*frac+c(0, 2.5, 0, 1),
       vertex.frame.color=NA)
  
  dev.off()
  
  
  # dev.off()
  
  ####
  
  
  
  
  
  
  # 
  # # Plot params
  # frac <- .65   # dialation of B relative to A and C
  # vlc <- 1.4    # vertex label cex
  # av <- 20     # vertex size intercept
  # bv <- 400     # vertex size slope
  # eaw <- 2.3    # edge arrow width
  # ae <- 6       # edge size intercept
  # be <- 200     # edge size slsope
  # 
  # continents <- unique( codes$continent[order(codes$continent)] )
  # circle_colors <- brewer.pal(11, "Spectral")[c(1,5,11,3,9,7)]
  # # circle_colors <- brewer.pal(6, "Spectral")
  # circle_colors_trans <- paste0(circle_colors, 95)
  # # library("igraph")
  # 
  # j <- 1
  # C <- as.matrix(A[[j]])
  # colnames(C) <- rownames(C) <- countries
  # i1 <- which(C != 0, arr.ind=TRUE)
  # EL <- cbind(countries[i1[,1]], countries[i1[,2]])
  # 
  # net <- graph_from_edgelist(EL, directed=TRUE)
  # E(net)$weight <- C[i1]                           
  # V(net)$label <- names(V(net))
  # V(net)$continent <- codes$continent[match(names(V(net)), codes$country)]
  # V(net)$color <- circle_colors_trans[ match(V(net)$continent, continents) ]
  # ordering <- V(net)$label[ order(V(net)$continent) ] #[keep]
  # L <- layout_in_circle(net, order=ordering)
  # # V(net)$label <- codes$country
  # # V(net)$color <- circle_colors_trans[match(codes$continent[ match(V(net)$label, codes$country) ], continents)]
  # E(net)$color <- circle_colors[ match(codes$continent[match(as_edgelist(net)[,1], countries)], continents) ]
  # E(net)$lty <- 3 - 2*(E(net)$weight >0)
  # 
  # V(net)$size <- av + bv*(rowSums(abs(C)) + colSums(abs(C)) )[V(net)$label]
  # E(net)$width <- ae + be*E(net)$weight
  # 
  # 
  # png(file.path(plotdir, paste0("B_circle_", j, ".png")), width=5, height=5, res=200, units="in")
  # par(mar = c(0, 0, 0, 0), mgp=c(1.7,.5,0), cex.lab=1.1, cex.axis=1.1) 
  # plot(net, 
  #      layout=L, 
  #      edge.curved=.25,
  #      edge.arrow.width	=eaw,
  #      vertex.label.font=2,
  #      vertex.label.family	= "Helvetica",
  #      vertex.label.color="gray30", 
  #      vertex.label.cex=vlc, 
  #      # vertex.label.degree=pi/2,
  #      # vertex.label.dist=2.1,
  #      vertex.frame.color=NA)
  # dev.off()
  # 
  # 
  # 
  # j <- 2
  # C <- as.matrix(A[[j]])
  # colnames(C) <- rownames(C) <- countries
  # i1 <- which(C != 0, arr.ind=TRUE)
  # EL <- cbind(countries[i1[,1]], countries[i1[,2]])
  # 
  # net <- graph_from_edgelist(EL, directed=TRUE)
  # E(net)$weight <- C[i1]                           
  # V(net)$label <- names(V(net))
  # V(net)$continent <- codes$continent[match(names(V(net)), codes$country)]
  # V(net)$color <- circle_colors_trans[ match(V(net)$continent, continents) ]
  # ordering <- V(net)$label[ order(V(net)$continent) ] #[keep]
  # L <- layout_in_circle(net, order=ordering)
  # # V(net)$label <- codes$country
  # # V(net)$color <- circle_colors_trans[match(codes$continent[ match(V(net)$label, codes$country) ], continents)]
  # 
  # E(net)$color <- circle_colors[ match(codes$continent[match(as_edgelist(net)[,1], countries)], continents) ]
  # E(net)$lty <- 3 - 2*(E(net)$weight >0)
  # 
  # V(net)$size <- (av + bv*(rowSums(abs(C)) + colSums(abs(C)) )[V(net)$label])*frac
  # E(net)$width <- (ae + be*E(net)$weight)*frac
  # 
  # 
  # png(file.path(plotdir, paste0("B_circle_", j, ".png")), width=5/frac, height=5/frac, res=200, units="in")
  # par(mar = c(0, 0, 0, 0), mgp=c(1.7,.5,0), cex.lab=1.1, cex.axis=1.1) 
  # plot(net, 
  #      layout=L, 
  #      edge.curved=.3,
  #      edge.arrow.width	=eaw*frac,
  #      vertex.label.font=2,
  #      vertex.label.family	= "Helvetica",
  #      vertex.label.color="gray30", 
  #      vertex.label.cex=vlc*frac*1.4,
  #      # vertex.label.degree=pi/2,
  #      # vertex.label.dist=1.2,
  #      vertex.frame.color=NA)
  # dev.off()
  # # edge.lty=E(net)$edge.lty, 
  # 
  # # edge.width=E(net)$edge.width,
  # # vertex.color=V(net)$vertex.color, 
  # 
  # 
  # # vertex.size = V(net)$vertex.size, 
  # 
  # 
  # 
  # 
  # circle_colors <- brewer.pal(11, "Spectral")[c(5,1,9,11)]  
  # circle_colors_trans <- paste0(circle_colors, 95)
  # 
  # j <- 3
  # C <- as.matrix(A[[j]])
  # colnames(C) <- rownames(C) <- resp
  # i1 <- which(C != 0, arr.ind=TRUE)
  # EL <- cbind(resp[i1[,1]], resp[i1[,2]])
  # 
  # 
  # net <- graph_from_edgelist(EL, directed=TRUE)
  # E(net)$weight <- C[i1]                           
  # V(net)$label <- names(V(net))
  # # V(net)$continent <- codes$continent[match(names(V(net)), codes$country)]
  # V(net)$color <- circle_colors_trans
  # ordering <- order(V(net)$label)
  # L <- layout_in_circle(net, order=c("vn", "mn", "vp", "mp"))
  # # V(net)$label <- codes$country
  # # V(net)$color <- circle_colors_trans[match(codes$continent[ match(V(net)$label, codes$country) ], continents)]
  # E(net)$color <- circle_colors[match(as_edgelist(net)[,1], V(net)$label)]
  # E(net)$lty <- 3 - 2*(E(net)$weight >0)
  # 
  # V(net)$size <- av + bv*(rowSums(abs(C)) + colSums(abs(C)) )[V(net)$label]
  # E(net)$width <- ae + be*E(net)$weight
  # 
  # 
  # png(file.path(plotdir, paste0("B_circle_", j, ".png")), width=5, height=5, res=200, units="in")
  # par(mar = c(0, 0, 0, 0), mgp=c(1.7,.5,0), cex.lab=1.1, cex.axis=1.1) 
  # plot(net, 
  #      layout=L, 
  #      edge.curved=.25,
  #      edge.arrow.width	=eaw,
  #      vertex.label.font=2,
  #      vertex.label.family	= "Helvetica",
  #      vertex.label.color="gray30", 
  #      vertex.label.cex=vlc, 
  #      # vertex.label.degree=pi/2,
  #      # vertex.label.dist=2.1,2
  #      vertex.frame.color=NA)
  # dev.off()
  
  ####


}
####

