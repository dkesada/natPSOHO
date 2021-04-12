#' R6 class that defines the PSO controller
#' 
#' The controller will encapsulate the particles and run the algorithm
natPsoCtrl <- R6::R6Class("natPsoCtrl",
  public = list(
    #' @description 
    #' Constructor of the 'natPsoCtrl' class
    #' @param nodes a vector with the names of the nodes
    #' @param max_size maximum number of timeslices of the DBN
    #' @param n_inds number of particles that the algorithm will simultaneously process
    #' @param n_it maximum number of iterations of the pso algorithm
    #' @param in_cte parameter that varies the effect of the inertia
    #' @param gb_cte parameter that varies the effect of the global best
    #' @param lb_cte parameter that varies the effect of the local best
    #' @param v_probs vector that defines the random velocity initialization probabilities
    #' @param p parameter of the truncated geometric distribution for sampling edges
    #' @param r_probs vector that defines the range of random variation of gb_cte and lb_cte
    #' @param cte boolean that defines whether the parameters remain constant or vary as the execution progresses
    #' @return A new 'natPsoCtrl' object
    initialize = function(nodes, max_size, n_inds, n_it, in_cte, gb_cte, lb_cte,
                          v_probs, p, r_probs, cte){
      #initial_size_check(size) --ICO-Merge
      # Missing security checks --ICO-Merge
      
      ordering <- grep("_t_0", nodes, value = TRUE) 
      private$initialize_particles(nodes, ordering, max_size, n_inds, v_probs, p)
      private$gb_scr <- -Inf
      private$n_it <- n_it
      private$in_cte <- in_cte
      private$gb_cte <- gb_cte
      private$lb_cte <- lb_cte
      private$r_probs <- r_probs
      private$cte <- cte
      if(!cte){
        private$in_var <- in_cte / n_it # Decrease inertia
        private$gb_var <- (1-gb_cte) / n_it # Increase gb
        private$lb_var <- lb_cte / n_it # Decrease gb
      }
    },
    
    #' @description 
    #' Getter of the cluster attribute
    #' @return the cluster attribute
    get_cl = function(){return(private$cl)},
    
    #' @description 
    #' Transforms the best position found into a bn structure and returns it
    #' @return the size attribute
    get_best_network = function(){return(private$gb_ps$bn_translate())},
    
    #' @description 
    #' Main function of the pso algorithm.
    #' @param dt the dataset from which the structure will be learned
    run = function(dt){
      # Missing security checks --ICO-Merge
      private$evaluate_particles(dt)
      pb <- utils::txtProgressBar(min = 0, max = private$n_it, style = 3)
      # Main loop of the algorithm.
      for(i in 1:private$n_it){
        # Inside loop. Update each particle
        for(p in private$parts)
          p$update_state(private$in_cte, private$gb_cte, private$gb_ps, private$lb_cte, private$r_probs)
        
        if(!private$cte)
          private$adjust_pso_parameters()
        
        private$evaluate_particles(dt)
        utils::setTxtProgressBar(pb, i)
      }
      close(pb)
    }
  ),
  private = list(
    #' @field parts list with all the particles in the algorithm
    parts = NULL,
    #' @field cl cluster for the parallel computations
    cl = NULL,
    #' @field n_it maximum number of iterations of the pso algorithm
    n_it = NULL,
    #' @field in_cte parameter that varies the effect of the inertia
    in_cte = NULL,
    #' @field gb_cte parameter that varies the effect of the global best
    gb_cte = NULL,
    #' @field lb_cte parameter that varies the effect of the local best
    lb_cte = NULL,
    #' @field b_ps global best position found
    gb_ps = NULL,
    #' @field b_scr global best score obtained
    gb_scr = NULL,
    #' @field r_probs vector that defines the range of random variation of gb_cte and lb_cte
    r_probs = NULL,
    #' @field cte boolean that defines whether the parameters remain constant or vary as the execution progresses
    cte = NULL,
    #' @field in_var decrement of the inertia each iteration
    in_var = NULL,
    #' @field gb_var increment of the global best parameter each iteration
    gb_var = NULL,
    #' @field lb_var increment of the local best parameter each iteration
    lb_var = NULL,
    
    #' @description 
    #' If the names of the nodes have "_t_0" appended at the end, remove it
    #' @param ordering a vector with the names of the nodes in t_0
    #' @return the ordering with the names cropped
    crop_names = function(ordering){
      #sapply(ordering, function(x){gsub("_t_0", "", x)}, USE.NAMES = F)
      crop_names_cpp(ordering)
    },
    
    #' @description 
    #' Initialize the particles for the algorithm to random positions and velocities.
    #' @param nodes a vector with the names of the nodes
    #' @param ordering a vector with the names of the nodes in t_0
    #' @param max_size maximum number of timeslices of the DBN
    #' @param n_inds number of particles that the algorithm will simultaneously process
    #' @param v_probs vector that defines the random velocity initialization probabilities
    #' @param p parameter of the truncated geometric distribution for sampling edges
    initialize_particles = function(nodes, ordering, max_size, n_inds, v_probs, p){
      #private$parts <- parallel::parLapply(private$cl,1:n_inds, function(i){Particle$new(ordering, size)})
      ordering_raw <- private$crop_names(ordering)
      private$parts <- vector(mode = "list", length = n_inds)
      
      # private$parts <- init_list_cpp(natParticle$new, n_inds, nodes, ordering, ordering_raw, max_size, v_probs, p) # Slower than pure R
      
      for(i in 1:n_inds)
        private$parts[[i]] <- natParticle$new(nodes, ordering, ordering_raw, max_size, v_probs, p)
    },
    
    #' @description 
    #' Evaluate the particles and update the global best
    #' @param dt the dataset used to evaluate the position
    evaluate_particles = function(dt){
      for(p in private$parts){
        scr <- p$eval_ps(dt)
        if(scr > private$gb_scr){
          private$gb_scr <- scr
          private$gb_ps <- p$get_ps()
        }
      }
    },
    
    #' @description 
    #' Modify the PSO parameters after each iteration
    adjust_pso_parameters = function(){
      private$in_cte <- private$in_cte - private$in_var
      private$gb_cte <- private$gb_cte + private$gb_var
      private$lb_cte <- private$lb_cte - private$lb_var
    }
  )
)