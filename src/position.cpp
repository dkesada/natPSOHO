#include "include/position.h"

//' Create a natural causal list from a DBN. This is the C++ backend of the function.
//' 
//' @param cl an initialized causality list
//' @param net a dbn object treated as a list of lists
//' @param ordering a vector with the names of the variables in order
//' @return the natCauslist equivalent to the DBN
// [[Rcpp::export]]
Rcpp::NumericVector create_natcauslist_cpp(Rcpp::NumericVector &cl, Rcpp::List &net, StringVector &ordering) {
  Rcpp::List aux;
  Rcpp::StringVector parents;
  std::string node;
    
  // Translation into natural causal list
  for(unsigned int i = 0; i < ordering.size(); i++){
    node = ordering[i];
    aux = net[node];
    parents = aux["parents"];

    for(unsigned int j = 0; j < parents.size(); j++){
      node = parents[j];
      insert_node_natcl(cl, ordering, node, i);
    }
  }
  
  return cl;
}

//' Create a matrix with the arcs defined in a causlist object
//' 
//' @param cl a causal list
//' @param ordering a list with the order of the variables in t_0
//' @param rows number of arcs in the network
//' @return a StringMatrix with the parent nodes and the children nodes
// [[Rcpp::export]]
Rcpp::CharacterMatrix cl_to_arc_matrix_cpp(const Rcpp::NumericVector &cl, Rcpp::CharacterVector &ordering,
                                           unsigned int rows){
  Rcpp::StringMatrix res (rows, 2);
  int slice, j, k;
  k = 0;
  
  for(unsigned int i = 0; i < cl.size(); i++){
    slice = cl[i];
    j = 1;

    while(slice > 0){
      if(slice % 2 == 1){
        include_arc(res, ordering, i, j, k);
      }
      
      slice = slice >> 1;
      j++;
    }
  }
  
  return res;
}

//' Add a velocity to a position
//' 
//' @param cl the position's causal list
//' @param vl the velocity's causal list
//' @param n_arcs number of arcs present in the position
//' @return a list with the modified position and the new number of arcs
// [[Rcpp::export]]
Rcpp::List pos_plus_vel_cpp(Rcpp::List &cl, Rcpp::List &vl, int n_arcs){
  Rcpp::List slice_cl, slice_vl, cu_cl, cu_vl, pair_cl, pair_vl;
  Rcpp::NumericVector dirs_cl, dirs_vl;
  Rcpp::List res (2);
  
  for(unsigned int i = 0; i < cl.size(); i++){
    slice_cl = cl[i];
    slice_vl = vl[i];
    
    for(unsigned int j = 0; j < slice_cl.size(); j++){
      pair_cl = slice_cl[j];
      pair_vl = slice_vl[j];
      dirs_cl = pair_cl[1];
      dirs_vl = pair_vl[1];
      dirs_cl = add_dirs_vec(dirs_cl, dirs_vl, n_arcs);
      
      pair_cl[1] = dirs_cl;
      slice_cl[j] = pair_cl;
    }
    
    cl[i] = slice_cl;
  }
  
  res[0] = cl;
  res[1] = n_arcs;
  
  return res;
}

//' Initialize the particles
//' 
//' @param nodes the names of the nodes
//' @param size the size of the DBN
//' @param n_inds the number of particles
//' @return a list with the randomly initialized particles
// [[Rcpp::export]]
Rcpp::List init_list_cpp(Rcpp::StringVector nodes, unsigned int size, unsigned int n_inds){
  Rcpp::List res (n_inds);
  Environment psoho("package:natPsoho");
  Environment env = psoho["natPosition"];
  Function new_ps = env["new"];
  
  for(unsigned int i = 0; i < n_inds; i++){
    Environment ps;
    ps = new_ps(NULL, size, nodes);
    res[i] = ps;
  }
  
  return res;
}