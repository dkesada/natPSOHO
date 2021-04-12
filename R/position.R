#' R6 class that defines DBNs as vectors of natural numbers
#' 
#' A natPosition represents a single HO-DBN structure with a vector. Its function
#' is to encode the solutions in the PSO framework. Each particle will have a 
#' position.
natPosition <- R6::R6Class("natPosition", 
  inherit = natCauslist,
  public = list(
    #' @description 
    #' Constructor of the 'natPosition' class
    #' @param nodes a vector with the names of the nodes
    #' @param ordering a vector with the names of the nodes in t_0
    #' @param ordering_raw a vector with the names of the nodes without the appended "_t_0"
    #' @param max_size Maximum number of timeslices of the DBN
    #' @param p the parameter of the sampling truncated geometric distribution
    #' If lesser or equal to 0, a uniform distribution will be used instead. 
    #' @return A new 'natPosition' object
    #' @importFrom dbnR fold_dt
    initialize = function(nodes, ordering, ordering_raw, max_size, p = 0.06){
      #initial_size_check(size) --ICO-Merge
      
      super$initialize(ordering, ordering_raw)
      private$nodes <- nodes
      private$max_size <- max_size
      private$cl <- private$generate_random_position(length(ordering), p)
      private$n_arcs <- private$recount_arcs()
      private$p <- p
    },
    
    get_n_arcs = function(){return(private$n_arcs)},
    
    #' @description 
    #' Translate the vector into a DBN network
    #' 
    #' Uses this object private cl and transforms it into a DBN.
    #' @return a dbn object
    bn_translate = function(){
      arc_mat <- cl_to_arc_matrix_cpp(private$cl, private$ordering_raw, private$n_arcs)
      net <- bnlearn::empty.graph(private$nodes, check.args = FALSE)
      bnlearn::arcs(net, check.cycles = FALSE, check.illegal = FALSE, check.bypass = TRUE) <- arc_mat
      
      return(net)
    },
    
    #' @description 
    #' Add a velocity to the position
    #' 
    #' Given a natVelocity object, add it to the current position.
    #' @param vl a natVelocity object
    add_velocity = function(vl){
      private$n_arcs <- nat_pos_plus_vel_cpp(private$cl, vl$get_cl(), vl$get_cl_neg(), private$n_arcs)
    }
  ),
  
  private = list(
    #' @field n_arcs Number of arcs in the network
    n_arcs = NULL,
    #' @field max_size Maximum number of timeslices of the DBN
    max_size = NULL,
    #' @field p Parameter of the sampling truncated geometric distribution
    p = NULL,
    #' @field nodes Names of the nodes in the network
    nodes = NULL,
    
    #' @description 
    #' Return the static node ordering
    #' 
    #' This function takes as input a dbn and return the node ordering of the
    #' variables inside a timeslice. This ordering is needed to understand a
    #' position vector.
    #' @param net a dbn or dbn.fit object
    #' @return the ordering of the nodes in t_0
    dbn_ordering = function(net){
      return(grep("t_0", names(net$nodes), value = TRUE))
    },
    
    #' @description 
    #' Translate a DBN into a position vector
    #' 
    #' This function takes as input a network from a DBN and transforms the 
    #' structure into a vector of natural numbers if it is a valid DBN. Valid 
    #' DBNs have only inter-timeslice edges and only allow variables in t_0 to 
    #' have parents.
    #' @param net a dbn object
    cl_translate = function(net){
      private$cl <- create_natcauslist_cpp(private$cl, net$nodes, private$ordering)
    },
    
    #' @description 
    #' Generates a random position
    #' 
    #' This function takes as input the number of variables, the maximum size
    #' and the parameter p and returns a random position with arcs 
    #' sampled either from the uniform distribution or from a truncated 
    #' geometric distribution. Much faster than the binary implementation with
    #' lists of lists and random bn generation into translation.
    #' @param n_vars the number of variables in t_0
    #' @param p the parameter of the truncated geometric sampler. If lesser or
    #' equal to 0, a uniform distribution will be used instead.
    #' @return a random position
    generate_random_position = function(n_vars, p){
      res <- init_cl_cpp(n_vars * n_vars)
      
      if(p <= 0){
        res <- floor(runif(n_vars * n_vars, 0, 2^(private$max_size - 1)))
      }
      
      else{
        for(i in 1:length(res))
          res[i] <- trunc_geom(p, 2^(private$max_size - 1))
      }
      
      return(res)
    },
    
    #' @description 
    #' Recount the number of arcs in the cl
    #' @return the number of arcs
    recount_arcs = function(){
      private$n_arcs <- 0
      for(i in 1:length(private$cl))
        private$n_arcs <- private$n_arcs + bitcount(private$cl[i])
      
      return(private$n_arcs)
    }
    
  )
)


