#' R6 class that defines DBNs as causality lists
#' 
#' A causality list has a list with causal units, a size representing the
#' Markovian order of the network and a specific node ordering.
natPosition <- R6::R6Class("natPosition", 
  inherit = natCauslist,
  public = list(
    #' @description 
    #' Constructor of the 'natPosition' class
    #' @param net dbn or dbn.fit object defining the network
    #' @param max_size Maximum number of timeslices of the DBN
    #' @param nodes A list with the names of the nodes in t_0 in the network
    #' If its not null, a random position will be generated for those nodes.
    #' @param p the parameter of the sampling truncated geometric distribution
    #' If lesser or equal to 0, a uniform distribution will be used instead. 
    #' @return A new 'natPosition' object
    initialize = function(net, max_size, nodes = NULL, p = 0.06){
      #initial_size_check(size) --ICO-Merge
      
      if(!is.null(nodes)){
        #initial_nodes_check(nodes)
        super$initialize(nodes)
        tmp <- matrix(0, nrow = max_size, ncol = length(nodes))
        colnames(tmp) <- nodes
        private$nodes <- names(dbnR::fold_dt(as.data.frame(tmp), max_size))
        private$cl <- private$generate_random_position(length(nodes), max_size, p)
        private$n_arcs <- private$recount_arcs()
        
      }
      else{
        #initial_dbn_check(net) --ICO-Merge
        initial_dbn_to_causlist_check(net)
        super$initialize(private$dbn_ordering(net))
        private$nodes <- names(net$nodes)
        private$n_arcs <- dim(net$arcs)[1]
        private$cl_translate(net)
      }
      
      private$max_size <- max_size
      private$p <- p
    },
    
    get_n_arcs = function(){return(private$n_arcs)},
    
    #' @description 
    #' Translate the causality list into a DBN network
    #' 
    #' Uses this object private causality list and transforms it into a DBN.
    #' @return a dbn object
    bn_translate = function(){
      arc_mat <- cl_to_arc_matrix_cpp(private$cl, private$ordering_raw, private$n_arcs)
      
      net <- bnlearn::empty.graph(private$nodes)
      bnlearn::arcs(net) <- arc_mat
      
      return(net)
    },
    
    #' @description 
    #' Add a velocity to the position
    #' 
    #' Given a Velocity object, add it to the current position.
    #' @param vl a Velocity object
    add_velocity = function(vl){
      res = pos_plus_vel_cpp(private$cl, vl$get_cl(), private$n_arcs)
      private$cl = res[[1]]
      private$n_arcs = res[[2]]
    },
    
    #' @description 
    #' Given another position, returns the velocity that gets this position to the
    #' other.
    #' 
    #' @param ps a Position object
    #' return the Velocity that gets this position to the new one
    subtract_position = function(ps){
      res <- Velocity$new(private$ordering, private$size)
      res$subtract_positions(self, ps)
      
      return(res)
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
    #' causal list.
    #' @param net a dbn or dbn.fit object
    #' @return the ordering of the nodes in t_0
    dbn_ordering = function(net){
      return(grep("t_0", names(net$nodes), value = TRUE))
    },
    
    #' @description 
    #' Translate a DBN into a causality list
    #' 
    #' This function takes as input a network from a DBN and transforms the 
    #' structure into a causality list if it is a valid DBN. Valid DBNs have only
    #' inter-timeslice edges and only allow variables in t_0 to have parents.
    #' @param net a dbn object
    #' @return a causlist object
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
    #' @param max_size the maximum size of the DBN
    #' @param p the parameter of the truncated geometric sampler. If lesser or
    #' equal to 0, a uniform distribution will be used instead.
    #' @return a random position
    generate_random_position = function(n_vars, max_size, p){
      res <- rep(0, n_vars * n_vars)
      
      if(p <= 0){
        res <- floor(runif(n_vars * n_vars, 0, max_size + 1))
      }
      
      else{
        for(i in 1:length(res))
          res[i] <- trunc_geom(p, max_size)
      }
      
      return(res)
    },
    
    #' @description 
    #' Recount the number of arcs in the cl
    #' @return the number of arcs
    recount_arcs = function(){
      res <- 0
      for(i in 1:length(private$cl))
        res <- res + count_bits(private$cl[i])
      
      return(res)
    }
    
  )
)


