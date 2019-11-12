# Call all function files

# Libraries
require("RColorBrewer")
require("Matrix")  
require("countrycode")   
require("fields")  
require("glmnet")  
require("MASS")  
require("abind")  
require("foreach")
require("doMC")
require("foreign")  
require("utils")
require("devtools")
require("igraph")



# Download bilinear code and data
download.file("https://www.stat.washington.edu/people/pdhoff/Code/hoff_2014/YX", file.path(getwd(), "YX"))
download.file("https://www.stat.washington.edu/people/pdhoff/Code/hoff_2014/functions_als.r", file.path(getwd(), "functions_als.r"))
download.file("https://www.stat.washington.edu/people/pdhoff/Code/hoff_2014/tfunctions.r", file.path(getwd(), "tfunctions.r"))


# Functions
source("functions_als.r")
source("MLE_functions.R")
source("cv_functions.R")
source("sid_functions.R")
