<img src="rSPAMM/instructions/figures/Logos.png" align="right" height="50px"/>

# The rSPAMM package

#### T. A. Øigård and M. Biuw

# Overview
The *rSPAMM* package is a R-package for performing assessment of harp and hooded seals in the West Ice (along the coast of Greenland) and the East Ice. Assessment of other populations is possible as long as similar input data is available. The package fits available pup production data and models the total abundance of the population of interest. Various catch options and the impact of the future abundance can be explored and graphical tools to visualize the results is available. 



# Installation

Make sure that the *rSPAMM* folder is set to be the working directory.

There are two ways to install the package:

1. In *RStudio* you select *Build -> Install and Restart* or press down *ctrl+Shift+B*. 

2. Source the *install_rSPAMM.R* file in the *scripts* folder, e.g.,

```{r}
source("scripts/install_rSPAMM.R")
```


In order to load the package type:
```{r}
library(rSPAMM)
```

Instructions on how to use *rSPAMM* is found in the *howToUse_rSPAMM.pdf* file in the */instructions* folder.

