#' One hot encoder for natural numbers without the 0.
#' 
#' Given a natural number, return the natural number equivalent to its
#' one-hot encoding. Examples: 3 -> 100 -> 4, 5 -> 10000 -> 16
#' 
#' @param nat the natural number to convert
#' @return the converted number
one_hot <- function(nat){
  return(2^(nat-1))
}

#' Geometric distribution sampler truncated to a maximum
#' 
#' A geometric distribution sampler with probability 'p' restricted to values
#' inside [1, max]. Because of this restriction, very low values of 'p' 
#' coupled with low 'max' return increasingly uniform populations in 
#' the interval [1, max].
#' 
#' @param p the parameter of the geometric distribution
#' @param max the maximum value allowed to be sampled
#' @return the sampled value
#' @importFrom stats runif
trunc_geom <- function(p, max){
  return(floor(log(1 - runif(1)*(1 - (1 - p)^max)) / log(1 - p)))
}