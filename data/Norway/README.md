---
title: "Data structure and data files for the rSPAMM package"
author: "Tor Arne Oigard and Martin Biuw"
date: "10/31/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

## Data files
All data that is described below will be loaded and used by the model to fit the population dynamics. 

### Catch data
Catch data are found in the file *catch_data.dat* and the content of the file looks like this:

    ##     V1     V2    V3
    ## 1 1946 105631 64742
    ## 2 1947 105631 64742
    ## 3 1948 105631 64742
    ## 4 1949 105631 64742
    ## 5 1950 105631 64742
    ## 6 1951  67793 76448

The first column contains the year, the seccond column contains the catch of pups, and the third column contains the catch of the 1+ population.

### Pup production estimates
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

### Fecundity data
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

### Birth ogive
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

### Priors used
Priors used for the various parameters are specified in the *priors.dat* file. The content of the file looks like this:

    ##        V1    V2
    ## 1 1.0e+06 2e+05
    ## 2 9.0e-02 5e-02
    ## 3 2.7e-01 5e-02

Each row specify mean and the standard deviation of a Gaussian prior used for a given parameter. In the first colum the mean value is found and in the second column the standard deviation is found. In the first row the prior distribution of the initial population size $K$ is specified, in the seccond row the prior distribution of the 1+ mortality $M$ is specified, and in the third row the prior distribution of the pup mortality $M_0$ is specified.

