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
#' @param ps a position vector of natural numbers
#' @param ordering_raw the ordering of the variables
#' @param n_arcs the total number of arcs 
#' @param nodes the name of all the nodes in the network
#' @return a bn object
bn_translate_exp = function(ps, ordering_raw, n_arcs, nodes){
  arc_mat <- cl_to_arc_matrix_cpp(ps, ordering_raw, n_arcs)
  
  net <- bnlearn::empty.graph(nodes)
  bnlearn::arcs(net) <- arc_mat
  
  return(net)
}

#' Experimental function that recounts the number of arcs in the position
#' @param ps a position vector of natural numbers
#' @return the number of arcs
recount_arcs_exp = function(ps){
  n_arcs <- 0
  
  for(i in 1:length(ps))
    n_arcs <- n_arcs + bitcount(ps[i])
  
  return(n_arcs)
}

# Generates the names of the nodes in t_0 and in all the network
nodes_gen_exp <- function(ordering, size){
  res <- list(ordering_t_0 = NULL, nodes = NULL)
  
  res$ordering_t_0 <- sapply(ordering, function(x){paste0(x, "_t_0")}, USE.NAMES = F)
  tmp <- matrix(nrow = size, ncol = length(ordering))
  tmp <- as.data.table(tmp)
  names(tmp) <- ordering
  tmp <- dbnR::fold_dt(tmp, size)
  res$nodes <- names(tmp)
  
  return(res)
}

# Generates the names of n variables.
ordering_gen_exp <- function(n){
  res <- rep("", n)
  for(i in 1:n)
    res[i] <- paste0("X", i-1)
  
  return(res)
}

#' Experimental function that generates a random DBN and samples a dataset that defines it
#' 
#' @param n_vars number of desired variables per time-slice
#' @param size desired size of the networks
#' @param min_mu minimum mean allowed for the variables
#' @param max_mu maximum mean allowed for the variables
#' @param min_sd minimum standard deviation allowed for the variables
#' @param max_sd maximum standard deviation allowed for the variables
#' @param min_coef minimum coefficient allowed for the parent nodes
#' @param max_coef maximum coefficient allowed for the parent nodes
#' @param seed the seed of the experiment
#' @return the number of arcs
#' @import data.table
#' @export
generate_random_network_exp <- function(n_vars, size, min_mu, max_mu,
                                        min_sd, max_sd, min_coef, max_coef, seed = NULL){
  res <- list(net = NULL, f_dt = NULL)
  set.seed(seed)
  
  # First. we generate a random position and translate it into a DBN structure
  
  # Generate the names of the variables in the network
  ordering_raw <- ordering_gen_exp(n_vars)
  nodes_l <- nodes_gen_exp(ordering_raw, size)
  ps <- rep(0, n_vars * n_vars)
  
  # Generate a random position
  for(i in 1:length(ps))
    ps[i] <- floor(runif(1, 0, 2^(size-1)))
  
  n_arcs <- recount_arcs_exp(ps)
  res$net <- bn_translate_exp(ps, ordering_raw, n_arcs, nodes_l$nodes)
  
  # Second, we generate a dataset that represents the same relationships portrayed by the network
  dt <- as.data.table(matrix(nrow = 10000, ncol = length(nodes_l$nodes)))
  names(dt) <- nodes_l$nodes
  dt[, names(dt) := lapply(.SD, function(x){stats::rnorm(length(x),
                                                  runif(1, min_mu, max_mu), 
                                                  runif(1, min_sd, max_sd))})]
  
  # Apply the effects of the arcs in the dataset
  bn_arcs <- bnlearn::arcs(res$net)
  for(i in 1:n_arcs){
    from <- bn_arcs[i,1]
    to <- bn_arcs[i,2]
    dt[, (to) := get(to) + get(from) * runif(1, min_coef, max_coef)]
  }
  
  res$f_dt <- dt
  
  return(res)
}



