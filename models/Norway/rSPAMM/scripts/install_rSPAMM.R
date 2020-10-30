#!/usr/bin/env Rscript

# Creates documentation and installs the R-package


# Load library required for automatic documentation
library(devtools)
library(roxygen2)

# Generate documentation
document()

# Install library 
setwd("../")
install("rSPAMM")

# Get back to where you started
setwd("rSPAMM/")


