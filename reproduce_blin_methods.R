# Reproduce results from 
# Marrs, F.W., Campbell, B.C., Fosdick, B.K., Cranmer, S.J., and B{\"o}hmelt, T. (2018) 
# ``Inferring influence networks from longitudinal bipartite relational data.'' 
# Submitted August 2018.



# set working directory to download location, e.g.
setwd("~/Downloads")

# call functions and load libraries, 
# downloads and installs `arcdiagram` from github and function files from Hoff (2015).
source("all_blin_functions.R")   



# Run simulations for first cross-validation study
run_first_cv_sims(nsims=20, verbose=TRUE)
# runs for nsims = # of error realizations
# simulations in the paper go to 1,000, but may take a lot of time to run

# Make plots for first cross-validation study
plot_first_cv()



# Estimate BLIN model for the ICEWS data set
fit_blin_sid(verbose=TRUE)

# Make plots of estimated networks
plot_blin_sid()


