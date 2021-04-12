#' R6 class that defines a Particle in the PSO algorithm
#' 
#' A particle has a Position, a Velocity and a local best
natParticle <- R6::R6Class("natParticle",
 public = list(
   #' @description 
   #' Constructor of the 'natParticle' class
   #' @param nodes a vector with the names of the nodes
   #' @param ordering a vector with the names of the nodes in t_0
   #' @param ordering_raw a vector with the names of the nodes without the appended "_t_0"
   #' @param max_size maximum number of timeslices of the DBN
   #' @param v_probs vector of probabilities for the velocity sampling
   #' @param p parameter of the truncated geometric distribution 
   #' @return A new 'natParticle' object
   initialize = function(nodes, ordering, ordering_raw, max_size, v_probs, p){
     #initial_size_check(size) --ICO-Merge
     
     private$ps <- natPosition$new(nodes, ordering, ordering_raw, max_size, p)
     private$vl <- natVelocity$new(ordering, ordering_raw, max_size)
     private$vl$randomize_velocity(v_probs, p)
     private$vl_gb <- natVelocity$new(ordering, ordering_raw, max_size)
     private$vl_lb <- natVelocity$new(ordering, ordering_raw, max_size)
     private$lb <- -Inf
   },
   
   #' @description 
   #' Evaluate the score of the particle's position
   #' 
   #' Evaluate the score of the particle's position.
   #' Updates the local best if the new one is better.
   #' @param dt dataset to evaluate the fitness of the particle
   #' @return The score of the current position
   eval_ps = function(dt){
     struct <- private$ps$bn_translate()
     score <- bnlearn::score(struct, dt, type = "bge", check.args = F, targets = private$ps$ordering) # For now, unoptimized bge. Any Gaussian score could be used
     
     if(score > private$lb){
        private$lb <- score 
        private$lb_ps <- private$ps
     }
     
     return(score)
   },
   
   #' @description 
   #' Update the position of the particle with the velocity
   #' 
   #' Update the position of the particle given the constants after calculating
   #' the new velocity
   #' @param in_cte parameter that varies the effect of the inertia
   #' @param gb_cte parameter that varies the effect of the global best
   #' @param gb_ps position of the global best
   #' @param lb_cte parameter that varies the effect of the local best
   #' @param r_probs vector that defines the range of random variation of gb_cte and lb_cte
   update_state = function(in_cte, gb_cte, gb_ps, lb_cte, r_probs){ # max_vl = 20
      # 1.- Inertia of previous velocity
      private$vl$cte_times_velocity(in_cte)
      # 2.- Velocity from global best
      op1 <- gb_cte * runif(1, r_probs[1], r_probs[2])
      private$vl_gb$subtract_positions(private$ps, gb_ps)
      private$vl_gb$cte_times_velocity(op1)
      # 3.- Velocity from local best
      op2 <- lb_cte * runif(1, r_probs[1], r_probs[2])
      private$vl_lb$subtract_positions(private$ps, private$lb_ps)
      private$vl_lb$cte_times_velocity(op2)
      # 4.- New velocity
      private$vl$add_velocity(private$vl_gb)
      private$vl$add_velocity(private$vl_lb)
      # 5.- Reduce velocity if higher than maximum. Awful results when the limit is low, so dropped for now.
      # if(private$vl$get_abs_op() > max_vl)
      #    private$vl$cte_times_velocity(max_vl / private$vl$get_abs_op())
      # 6.- New position
      private$ps$add_velocity(private$vl)
      # 7.- If a node has more parents than the maximum, reduce them (TODO)
   },
   
   get_ps = function(){return(private$ps)},
   
   get_vl = function(){return(private$vl)},
   
   get_lb = function(){return(private$lb)},
   
   get_lb_ps = function(){return(private$lb_ps)}
 ),
 
 private = list(
   #' @field ps position of the particle
   ps = NULL,
   #' @field cl velocity of the particle
   vl = NULL,
   #' @field velocity that takes the particle to the global best
   vl_gb = NULL, # Just to avoid instantiating thousands of velocities
   #' @field velocity that takes the particle to the local best
   vl_lb = NULL,
   #' @field lb local best score obtained
   lb = NULL,
   #' @field lb_ps local best position found
   lb_ps = NULL
 )
)