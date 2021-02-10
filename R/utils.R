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

########### ICO-Merge: Delete the experimental functions

#' Experimental function that translates a natPosition vector into a DBN network.
#' 
#' This function will be used in the experiments, but should not appear in the 
#' final version of the package. Uses a position vector and transforms it into 
#' a DBN.
#' 
#' @param cl a position vector of natural numbers
#' @param ordering_raw the ordering of the variables
#' @param n_arcs the total number of arcs 
#' @param nodes the name of all the nodes in the network
#' @return a bn object
#' @export
bn_translate_exp = function(cl, ordering_raw, n_arcs, nodes){
  arc_mat <- cl_to_arc_matrix_cpp(cl, ordering_raw, private$n_arcs)
  
  net <- bnlearn::empty.graph(nodes)
  bnlearn::arcs(net) <- arc_mat
  
  return(net)
}

