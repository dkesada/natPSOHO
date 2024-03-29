#' R6 class that defines velocities in the PSO
#' 
#' The velocities will be defined as two natural vectors where each element in
#' them represents the arcs from a temporal family of nodes to a receiving
#' node. 1-bits in the binary representation of this number represent arc 
#' additions/deletions 
natVelocity <- R6::R6Class("natVelocity",
  inherit = natCauslist,
  public = list(
    #' @description 
    #' Constructor of the 'natVelocity' class. Only difference with the
    #' natCauslist one is that it has a negative cl attribute.
    #' @param ordering a vector with the names of the nodes in t_0
    #' @param ordering_raw a vector with the names of the nodes without the appended "_t_0"
    #' @param max_size maximum number of timeslices of the DBN
    #' @return A new 'natVelocity' object
    initialize = function(ordering, ordering_raw, max_size){
      super$initialize(ordering, ordering_raw)
      private$abs_op <- 0
      private$max_size <- max_size
      private$cl_neg <- init_cl_cpp(length(private$cl))
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
          if(p <= 0)
            private$cl[i] <- floor(runif(1, 0, 2^(private$max_size - 1)))
          else
            private$cl[i] <- trunc_geom(p, 2^(private$max_size - 1))
          private$abs_op <- private$abs_op + bitcount(private$cl[i])
        }
        else if (op[1] == 1){
          if(p <= 0)
            private$cl_neg[i] <- floor(runif(1, 0, 2^(private$max_size - 1)))
          else
            private$cl_neg[i] <- trunc_geom(p, 2^(private$max_size - 1))
          private$abs_op <- private$abs_op + bitcount(private$cl_neg[i])
        }
      }
    },
    
    #' @description 
    #' Given two positions, returns the velocity that gets the first position to the
    #' other one.
    #' 
    #' @param ps1 the origin natPosition object
    #' @param ps2 the objective natPosition object
    #' @return the natVelocity that gets the ps1 to ps2
    subtract_positions = function(ps1, ps2){
      private$abs_op <- nat_pos_minus_pos_cpp(ps1$get_cl(), ps2$get_cl(), private$cl, private$cl_neg)
    },
    
    #' @description 
    #' Add both velocities directions
    #' 
    #' @param vl a Velocity object
    add_velocity = function(vl){
      private$abs_op <- nat_vel_plus_vel_cpp(private$cl, private$cl_neg, vl$get_cl(), vl$get_cl_neg(), private$abs_op, vl$get_abs_op())
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
      
      # If k < 0, invert the cl and the cl_neg
      if(k < 0){ 
        tmp <- private$cl
        private$cl <- private$cl_neg
        private$cl_neg <- tmp
        k <- abs(k)
      }
      
      if(k == 0){
        private$cl <- init_cl_cpp(length(private$cl))
        private$cl_neg <- init_cl_cpp(length(private$cl_neg))
        private$abs_op <- 0
      }
      
      else
        private$abs_op = nat_cte_times_vel_cpp(k, private$cl, private$cl_neg, private$abs_op, private$max_size)
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