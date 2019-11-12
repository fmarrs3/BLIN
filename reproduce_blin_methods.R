# Reproduce results from 
# Marrs, F.W., Campbell, B.C., Fosdick, B.K., Cranmer, S.J., and B{\"o}hmelt, T. (2018) 
# ``Inferring influence networks from longitudinal bipartite relational data.'' 
# Submitted August 2018.



# set working directory to download location, e.g.
setwd("~/Desktop/BLIN_rep")

# call functions and load libraries, 
# downloads and installs `arcdiagram` from github and function files from Hoff (2015).
source("all_blin_functions.R")   



# Run simulations for cross-validation study
run_cv(nsims=11, qs=c(0.9))
# runs for nsims = # of error realizations
# simulations in the paper go to 100, but may take a lot of time to run
# qs is a vector of fractions of empty cells in coefficients. In paper, we ran for qs=c(0, 0.5, 0.9)


# Make plots for cross-validation study
plot_cv()



# Estimate BLIN model for the ICEWS data set
fit_blin_sid()

# Make plots of estimated networks
plot_blin_sid()


