---
title: "Data structure and data files for the rSPAMM package"
author: "Tor Arne Øigård and Martin Biuw"
date: "10/31/2020"
output: html_document
---


# Data sources used by the population dynamics model

To estimate the population size trajectory, the population dynamics model implemented in the *rSPAMM* package uses the following data sources:

1. Historical catch records 

2. Fecundity rates

3. Age specific proportions of mature females

4. Estimates of pup production

Two types of reproductive data are used in the model: information on the proportion of females that are mature at a given age (i.e., maturity ogive) and the proportion of mature females that are pregnant at a given year (i.e. fecundity rate). 

Will describe what type of data is used, which data files are needed, and the format the various data are stored.

All data described below will be loaded and used by the model to fit the population dynamics. 

# Data files
The following data files are needed. The format of each file is specified in the next section.

1. catch_data.dat

2. pup_production.dat

3. fecundity.dat

4. wgharp.ogi

5. wgharp.ogp

6. priors.dat

## Catch data
Catch data are found in the file *catch_data.dat* and the content of the content of the file looks like this:

    ##     V1     V2    V3
    ## 1 1946 105631 64742
    ## 2 1947 105631 64742
    ## 3 1948 105631 64742
    ## 4 1949 105631 64742
    ## 5 1950 105631 64742
    ## 6 1951  67793 76448

The first column contains the year, the seccond column contains the catch of pups, and the third column contains the catch of the 1+ population.

## Pup production estimates
Estimated pup production is found in the file *pup_production.dat*. The content of the file looks like this: 

    ##      V1     V2    V3
    ## 1  1998 286260 0.150
    ## 2  2000 322474 0.098
    ## 3  2000 339710 0.105
    ## 4  2002 330000 0.103
    ## 5  2003 328000 0.181
    ## 6  2004 231811 0.190
    ## 7  2004 234000 0.205
    ## 8  2005 122658 0.162
    ## 9  2008 123104 0.199
    ## 10 2009 157000 0.108
    ## 11 2010 163032 0.198
    ## 12 2013 128786 0.237

In the first column you see the survey year, in the seccond column you see the estimated pup production, and in the third the estimated CV of that survey is given.

## Fecundity data
Fecundity data is found in the *fecundity.dat* file. The content of the file looks like this:

    ##     V1   V2   V3
    ## 1 1990 0.84 0.06
    ## 2 1991 0.84 0.06
    ## 3 1992 0.84 0.06
    ## 4 1993 0.84 0.06
    ## 5 2006 0.68 0.06
    ## 6 2011 0.84 0.06
    ## 7 2018 0.91 0.03

In the first column you have the year for when the fecundity is estimated, in the seccond column you have the estimated fecundity, and in the thirst column you have the estimated standard error of the fecundity estimate.

## Birth ogive
In the files *wgharp.ogi* and *wgharp.ogp* you find the birth ogive data, and the from which time periods the birth ogive curves applies for. The content of the *wgharp.ogi* looks like

    ##    Age   P1   P2   P3   P4   P5
    ## 1    1 0.00 0.00 0.00 0.00 0.00
    ## 2    2 0.00 0.00 0.00 0.01 0.00
    ## 3    3 0.01 0.00 0.00 0.02 0.00
    ## 4    4 0.17 0.00 0.02 0.05 0.00
    ## 5    5 0.64 0.24 0.08 0.11 0.00
    ## 6    6 0.90 0.62 0.21 0.25 0.52
    ## 7    7 0.98 0.81 0.40 0.55 0.77
    ## 8    8 0.99 0.91 0.59 0.90 0.89
    ## 9    9 1.00 0.95 0.75 0.99 0.95
    ## 10  10 1.00 0.98 0.85 1.00 0.97
    ## 11  11 1.00 0.99 0.91 1.00 0.99
    ## 12  12 1.00 0.99 0.95 1.00 0.99
    ## 13  13 1.00 1.00 0.97 1.00 0.99
    ## 14  14 1.00 1.00 0.98 1.00 1.00
    ## 15  15 1.00 1.00 0.99 1.00 1.00
    ## 16  16 1.00 1.00 1.00 1.00 1.00
    ## 17  17 1.00 1.00 1.00 1.00 1.00
    ## 18  18 1.00 1.00 1.00 1.00 1.00
    ## 19  19 1.00 1.00 1.00 1.00 1.00
    ## 20  20 1.00 1.00 1.00 1.00 1.00

In the first column the age group is listed and in the following columns a birth ogive curve is listed for a given period. Every time a new birth ogive estimated is available, a new column is added. The content of the *wgharp.ogp* looks like 

    ##   Pstart Pstop
    ## 1   1962  1972
    ## 2   1976  1985
    ## 3   1988  1993
    ## 4   2006  2006
    ## 5   2018  2018

The first column denotes the start of the time period a given birth ogive curve applies for, and the second column denotes the end of the time period the birth ogive curve is valid. The reason for this is that usually data for several years are pooled together when estimating the birth ogive curve in order to ensure enough data to obtain a reasonable estimated of the birth ogive. Note that the number of rows in the *wgharp.ogp* file corresponds to the number of columns in the *wgharp.ogi* file (minus 1 since one column represents the age group).

## Priors used
Priors used for the various parameters are specified in the *priors.dat* file. The content of the file looks like this:

    ##        V1    V2
    ## 1 1.0e+06 2e+05
    ## 2 9.0e-02 5e-02
    ## 3 2.7e-01 5e-02

Each row specify mean and the standard deviation of a Gaussian prior used for a given parameter. In the first colum the mean value is found and in the second column the standard deviation is found. In the first row the mean and the standard deviation of the Normal prior distribution of the initial population size $K$ is specified, in the seccond row the mean and the standard deviation of the Normal prior distribution of the 1+ mortality $M$ is specified, and in the third row the mean and the standard deviation of the Normalprior distribution of the pup mortality $M_0$ is specified.

# How to load the data files

When cloning the repository and installing the R package *rSPAMM*, a demo data set is installed. In addition a full data set is available in the *wk_WKSEALS-2020/data/Norway/* folder. 

## The demo data set

The demo data is reproductive data, catch data, pup production estimates and priors used for population dynamics modelling of the harp seal population in the East Ice (White Sea).

To load the demo data:
```{r eval = TRUE}
data("harpeastDemo")
```

The demo data is a list called `harpeast` containing two lists called `data` and `parameters`. The `data` list contains the data needed to fit the population dynamics model and the `parameters` list contains the parameters estimated by the model along with initial values of them.

```{r eval = TRUE}
names(harpeast$data)

names(harpeast$parameters)
```

To use the demo data set in the following examples it would be easiest to split the list in two separate lists `data` and `parameters`, i.e.,

```{r eval = TRUE}
data = harpeast$data
parameters = harpeast$parameters
```

## The full data set
The full data set is stored in the repositorys data folder (in the Norway folder).

It does not matter what is set to working directory, but when loading these data to use with the *rSPAMM* you might want to adjust the *dataPath* variable defined below. In this example it is assumed that the working directory is set to the *rSPAMM* model folder, i.e., that the working directory is the root folder of the *rSPAMM* package.

In order to load the data you have to specify which population you want to load the data for. The various alternatives are `harpeast`, `harpwest`, and `hooded` (which is the hooded seal population in the West Ice - the ice along the east coast of Greenland). In the examples that follows we will use the `harpeast` population. 

To load the data used you run the following function:
```{r eval=FALSE}
dataPath = "../../../../data/Norway/"
population = "harpeast"
dataFiles = paste0(dataPath,population)
data <- load.data(population = dataFiles)
```

```{r echo=FALSE}
data = harpeast$data
parameters = harpeast$parameters
```

To explore the data object you can
```{r}
names(data)
```

The data object is a list and you can further explore the actual values used by e.g.,
```{r}
data$Amax

data$pupProductionData

```

# Parameters to be estimed and how to load them
As briefly mentioned earlier the population dynamics model is described by three parameters. The initial population size $K$, the mortality of the 1+ population $M$, and the pup mortality $M_0$. The initial population size is the population size for the year the model is fit from, and this is determined by the availability of the catch data. The earliest catch data is from 1946, so the model is fit from 1946 and up to present. It also predicts the population dynamics into the future, and this is default set to 15 years. 

To load the initial values of the parameters to be estimated run:
```{r eval=FALSE}
paramPath = "../../../../data/Norway/"
population = "harpeast"
priorFile = paste0(paramPath,population)
parameters = load.initial.values(population = priorFile)
```
The parameters object is also a list and the content looks like this for the `harpeast` population
```{r}
parameters
```
Here `logK` is the log transformed initial population size, `M0tilde` is the logit transformed pup mortality, and `Mtilde` is the logit transformed 1+ mortality (seals of age 1 and greater). The reason that transformed parameters are estimated instead of the non transformed parameters is that the log transformation of the initial population ensures that the model provides a stricktly positive estimate of the initial population. The logit tranformation of the mortalities ensures that the estimates are bounded between 0 and 1. 

