#################################################################################
#########			CanHarp.Rdata				#########
#########	   Northwest Atlantic Harp seal population data		#########
#########	considered into the Canadian population dynamic model	#########
#################################################################################

The CanHarp object is a list with the following structure
List of 5
 $ removal   :'data.frame':	68 obs. of  7 variables:
  ..$ year       : int [1:68] 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 ...
  ..$ pupInclSL  : num [1:68] 219537 219448 199520 273297 357915 ...
  ..$ adultInclSL: num [1:68] 198636 141746 168120 151305 92466 ...
  ..$ totalInclSL: num [1:68] 418173 361194 367640 424602 450381 ...
  ..$ pupNoSL    : num [1:68] 207800 207712 186393 261523 347932 ...
  ..$ adultNoSL  : num [1:68] 117492 83358 98957 89164 54235 ...
  ..$ totalNoSL  : num [1:68] 325292 291070 285350 350687 402167 ...
 $ pup_prod  :'data.frame':	67 obs. of  4 variables:
  ..$ year    : int [1:67] 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 ...
  ..$ pup_prod: int [1:67] NA NA NA NA NA NA NA 235000 NA NA ...
  ..$ SE      : int [1:67] NA NA NA NA NA NA NA 120428 NA NA ...
  ..$ CV      : num [1:67] NA NA NA NA NA ...
 $ pregnancy :'data.frame':	250 obs. of  4 variables:
  ..$ years: int [1:250] 1952 1953 1954 1964 1965 1966 1967 1968 1969 1970 ...
  ..$ age  : int [1:250] 4 4 4 4 4 4 4 4 4 4 ...
  ..$ total: int [1:250] 0 0 4 11 30 7 10 27 25 13 ...
  ..$ npreg: int [1:250] NA NA 0 0 1 0 0 0 1 0 ...
 $ fecundity : logi NA
 $ parameters:List of 4
  ..$ InitPop   :'data.frame':	1 obs. of  3 variables:
  .. ..$ pup    : num 80000
  .. ..$ onePlus: num 8605516
  .. ..$ year   : num 1952
  ..$ PupM      :'data.frame':	1 obs. of  2 variables:
  .. ..$ mean: num 0.316
  .. ..$ sd  : num 0.079
  ..$ AdultM    :'data.frame':	1 obs. of  2 variables:
  .. ..$ mean: num 0.03
  .. ..$ sd  : num 0.02
  ..$ maturCurve: logi NA
  
#######  Details ####################################  

# removal #
Catch level from 1952 to 2019 provided as pup, adult (1+) and total catch including Struck and Lost (i.e. animal stuck but not recovered; Indicated with the suffix "InclSL"), or not including Struck and Lost (indicated with the suffix "NoSL")

# pup_prod #
Number of pups estimated from aerial surveys counts and Capture-Mark-Recapture (1960;1978; 1979; 1980; 1983; 1990; 1994; 1999; 2004; 2008; 2012; 2017)
Standard error and Confidence intervals around the estimates are also provided

# pregnancy #
Number of female pregnant per age class (4 to 7 years old, then 8+).
Data are from a sampling program conducted from 1952 to 2019. Sampling size is indicated in the column "total"

# fecundity #
Used in the Norwegian model but not available here - need to discuss how to convert pregnancy rate to fecundity or if it is even needed to be able to run these data into the Norwegian model.

# parameters #
	# InitPop #
	Initial population size given to the model in 1952 separated in pup and onePlus (note: in the canadian model, this is provided as numbers in 26 age classes)

	# PupM #
	Pup mortality (and standard deviation) estimated by the last run of the canadian population dynamics model
	
	# AdultM #
	Adult mortality value plus a standard deviation around this value. Note: In the last version of the canadian model, this value is fixed so the sd provided here is just for testing purpose and will need justification.

	# maturCurve #
	Used in the Norwegian model but not available here
