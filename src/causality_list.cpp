#include "include/causality_list.h"

// Insert an arc in the correspondent temporal family. 
// 
// @param cl a causality list
// @param ordering a list with the order of the variables in t_0
// @param node the node to insert
// @param i the causal unit in which to insert.
void insert_node_natcl(Rcpp::NumericVector &cl, const StringVector &ordering, std::string node, unsigned int i){
  Rcpp::StringVector tuple = find_name_and_index(node);
  std::string tmp;
  tmp = tuple[1];
  int idx = std::stoi(tmp);
  tmp = tuple[0];
  int ordering_idx = find_index(ordering, tmp);
  int arcs = cl[i * 3 + ordering_idx];
  
  idx = one_hot_cpp(idx);
  arcs = arcs | idx;
  cl[i * 3 + ordering_idx] = arcs;
}


