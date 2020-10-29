#' Inverse logit transformation - type 2
#'
#' Performs an inverse logit transforms the input number
#' @param x Number to be inverse logit transformed
#' @return ilogit_x Inverse logit transformed value of x
#' @keywords logit 
#' @export
#' @examples
#' ilogit1(x) 

ilogit1 <- function(x){
  return(1/(1+exp(-x)))
}

#' Logit transformation
#'
#' Logit transforms the input number
#' @param x Number to be logit transformed
#' @return logit_x Logit transformed value of x
#' @keywords logit 
#' @export
#' @examples
#' logit(x) 

logit <- function(x){
  return(log(x/(1-x)))
}

#' Inverse logit transformation - type 1
#'
#' Performs an inverse logit transforms the input number
#' @param x Number to be inverse logit transformed
#' @return ilogit_x Inverse logit transformed value of x
#' @keywords logit 
#' @export
#' @examples
#' ilogit(x)
ilogit <- function(x){
  return(exp(x)/(1+exp(x)))
}

