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
  for(int i = 0; i < ordering.size(); i++){
    node = ordering[i];
    aux = net[node];
    parents = aux["parents"];

    for(int j = 0; j < parents.size(); j++){
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
  
  for(int i = 0; i < cl.size(); i++){
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
//' @param vl the velocity's positive causal list
//' @param vl_neg velocity's negative causal list
//' @param n_arcs number of arcs present in the position. Remainder: can't return integers by reference, they get casted to 1 sized vectors
//' @return the new position by reference and the new number of arcs by return
// [[Rcpp::export]]
int nat_pos_plus_vel_cpp(Rcpp::NumericVector &cl, const Rcpp::NumericVector &vl, const Rcpp::NumericVector &vl_neg, int n_arcs){
  int pos, new_pos, vl_i, vl_neg_i, n_prev, n_post;
  
  for(int i = 0; i < cl.size(); i++){
    pos = cl[i];
    vl_i = vl[i];
    vl_neg_i = vl_neg[i];
    new_pos = pos | vl_i;
    new_pos = bitwise_sub(new_pos, vl_neg_i);
    n_prev = bitcount(pos);
    n_post = bitcount(new_pos);
    
    cl[i] = new_pos;
    n_arcs += n_post - n_prev;
  }
  
  return n_arcs;
}