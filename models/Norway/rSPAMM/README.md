<img src="vignettes/figures/Logos.png" align="right" height="50px"/>

# The rSPAMM package

#### T. A. Øigård and M. Biuw

# Overview
The *rSPAMM* package is a R-package for performing assessment of harp and hooded seals in the West Ice (along the coast of Greenland) and the East Ice. Assessment of other populations is possible as long as similar input data is available. The package fits available pup production data and models the total abundance of the population of interest. Various catch options and the impact of the future abundance can be explored and graphical tools to visualize the results is available. 


Instructions on how to use *rSPAMM* is found in the vignette which will be available when installing the R package. 


# Installation

To run the model and the various functions for assessment it is recommended to install the *rSPAMM* package as a R package. It is also possible to run the model by sourcing the relevant files in the /R folder.


To install the R-package directly from the repository use the following command (note: the R-package devtools has to be installed first)
``` r
devtools::install_github("https://github.com/NorskRegnesentral/rSPAMM.git", build_vignettes = TRUE)
``` 

In order to load the package type:
```{r}
library(rSPAMM)
```

Instructions on how to use *rSPAMM* is found in the vignette.

To load the vignette type:
```{r}
vignette("howToUse_rSPAMM",package = "rSPAMM")
```

The vignette will be opened in a system specified viewer. For a HTML version of the vignette (recommended) use

```{r}
browseVignettes()
```
and scroll down to the `rSPAMM` package.