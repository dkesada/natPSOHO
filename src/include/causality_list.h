#ifndef Rcpp_head
#define Rcpp_head
#include <Rcpp.h>
using namespace Rcpp;
#endif

#include "utils.h"

#ifndef nat_cl_op
#define nat_cl_op
void insert_node_natcl(Rcpp::NumericVector &cl, const StringVector &ordering, std::string node, unsigned int i);
#endif