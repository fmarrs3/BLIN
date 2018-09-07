# BLIN
This directory contains code to reproduce the results from "Inferring influence networks from longitudinal bipartite relational data" (2018) by Frank W. Marrs, Benjamin C. Campbell, Bailey K. Fosdick, Skyler J. Cranmer, and Tobias Böhmelt.


## Original Code
`reproduce_blin_methods.R` -- Master file. Run this to reproduce the results from Marrs et. al. (2018). 

`all_blin_functions.R` -- Load supporting functions and libraries. Requires internet connection to download data and existing code.

`MLE_functions.R` -- Functions to estimate BLIN model. 

`sid_functions.R` -- Functions to estimate BLIN model for ICEWS data set.

`cv_functions.R` -- Cross-validation functions. 

`misspec_run_cv.R` -- Cross-validation script. 


## Supporting code and data
To analyze the ICEWS data, the data must be downloaded. In addition, to estimate the bilinear model, we require the functions from "Multilinear tensor regression for longitudinal relational data" (2015) by Peter D. Hoff. To create the influence network plots, the R package `arcdiagram` is required, which is hosted on github (rather than CRAN). Running the master file will automatically download these files (by sourcing `all_blin_functions.R`). Lastly, some R packages may be required; all packages are listed  in `all_blin_functions.R`.



