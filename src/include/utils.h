#ifndef Rcpp_head
#define Rcpp_head
#include <Rcpp.h>
using namespace Rcpp;
#endif

#ifndef utils_op
#define utils_op

#include <regex>
#include <random>

int find_index(std::string node);
Rcpp::StringVector rename_nodes_cpp(Rcpp::StringVector &nodes, unsigned int size);
Rcpp::StringVector rename_slices(Rcpp::StringVector nodes, unsigned int size);
Rcpp::List random_directions(Rcpp::NumericVector probs, unsigned int size);
#endif