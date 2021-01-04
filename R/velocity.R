#' R6 class that defines velocities affecting causality lists in the PSO
#' 
#' The velocities will be defined as a causality list where each element in
#' a causal unit is a pair (v, node) with v being either 0, 1 or -1. 0 means 
#' that arc remained the same, 1 means that arc was added and -1 means that arc 
#' was deleted.
natVelocity <- R6::R6Class("natVelocity",
  inherit = natCauslist,
  public = list(
    #' @description 
    #' Constructor of the 'natVelocity' class. Only difference with the
    #' natCauslist one is that it has a negative cl attribute.
    #' @param ordering a vector with the names of the nodes in t_0
    #' @return A new 'natVelocity' object
    initialize = function(ordering, max_size){
      super$initialize(ordering)
      private$abs_op <- 0
      private$max_size <- max_size
      private$cl_neg <- private$cl
    },
    
    get_cl_neg = function(){return(private$cl_neg)},
    
    #' @description 
    #' Getter of the abs_op attribute.
    #' 
    #' return the number of operations that the velocity performs
    get_abs_op = function(){return(private$abs_op)},
    
    #' @description 
    #' Setter of the abs_op attribute. Intended for inside use only. 
    #' This should be a 'protected' function in Java-like OOP, but there's no 
    #' such thing in R6. This function should not be used from outside the
    #' package.
    #' 
    #' @param n the new number of operations that the velocity performs
    set_abs_op = function(n){private$abs_op = n},
    
    #' @description 
    #' Randomizes the Velocity's directions.
    #' 
    #' @param probs the weight of each value {-1,0,1}. They define the probability that each of them will be picked 
    #' @param p the parameter of the geometric distribution
    randomize_velocity = function(probs = c(10, 65, 25), p = 0.06){
      numeric_prob_vector_check(probs)
      
      for(i in 1:length(private$cl)){
        op <- rmultinom(n = 1, size = 1, prob = probs)
        if(op[3] == 1){
          private$cl[i] <- trunc_geom(p, 2^(private$max_size - 1))
          private$abs_op <- private$abs_op + count_bits(private$cl[i])
        }
        else if (op[1] == 1){
          private$cl_neg[i] <- trunc_geom(p, 2^(private$max_size - 1))
          private$abs_op <- private$abs_op + count_bits(private$cl_neg[i])
        }
      }
    },
    
    #' @description 
    #' Given a position, returns the velocity that gets this position to the
    #' other.
    #' 
    #' @param ps a Position object
    #' return the Velocity that gets this position to the new one
    subtract_positions = function(ps1, ps2){
      res <- pos_minus_pos_cpp(ps1$get_cl(), ps2$get_cl(), private$cl)
      
      private$cl <- res[[1]]
      private$abs_op <- res[[2]]
    },
    
    #' @description 
    #' Add both velocities directions
    #' 
    #' @param vl a Velocity object
    add_velocity = function(vl){
      
      res <- vel_plus_vel_cpp(private$cl, vl$get_cl(), private$abs_op)
      
      private$cl <- res[[1]]
      private$abs_op <- res[[2]]
    },
    
    #' @description 
    #' Multiply the Velocity by a constant real number
    #' 
    #' This function multiplies the Velocity by a constant real number. 
    #' It is non deterministic by definition. When calculating k*|V|, the 
    #' result will be floored and bounded to the set [-max_op, max_op], where max_op
    #' is the maximum number of arcs that can be present in the network.
    #' 
    #' @param k a real number
    cte_times_velocity = function(k){
      # initial_numeric_check(k) --ICO-Merge
      
      if(k == 0){
        private$cl <- initialize_cl_cpp(private$ordering, private$size)
        private$abs_op <- 0
      }
      
      else{
        max_op <- (private$size - 1) * length(private$ordering) * length(private$ordering)
        res = cte_times_vel_cpp(k, private$cl, private$abs_op, max_op)
        
        private$cl = res[[1]]
        private$abs_op = res[[2]]
      }
    }
  ),
  private = list(
    #' @field abs_op Total number of operations 1 or -1 in the velocity
    abs_op = NULL,
    #' @field max_size Maximum number of timeslices of the DBN
    max_size = NULL,
    #' @field cl_neg Negative part of the velocity
    cl_neg = NULL
  )
)