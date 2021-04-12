#' R6 class that defines causal lists in the PSO
#' 
#' The causal lists will be the base of the positions and the velocities
#' in the pso part of the algorithm. They will not have the same structure
#' as their binary counterparts, but their class skeleton will serve as a
#' base.
natCauslist <- R6::R6Class("natCauslist",
  public = list(
    #' @description 
    #' Constructor of the 'natCauslist' class
    #' @param ordering a vector with the names of the nodes in t_0
    #' @param ordering_raw a vector with the names of the nodes without the appended "_t_0"
    #' @return A new 'natCauslist' object
    initialize = function(ordering, ordering_raw){
      #initial_size_check(size) --ICO-Merge
      
      private$ordering <- ordering
      private$ordering_raw <- ordering_raw
      private$cl <- init_cl_cpp(length(ordering) * length(ordering))
    },
    
    get_cl = function(){return(private$cl)},
    
    get_ordering = function(){return(private$ordering)}
    
  ),
  private = list(
    #' @field cl List of causal units
    cl = NULL,
    #' @field ordering String vector defining the order of the nodes in t_0
    ordering = NULL,
    #' @field ordering String vector defining the order of the nodes without the appended "_t_0"
    ordering_raw = NULL
    
  )
)
