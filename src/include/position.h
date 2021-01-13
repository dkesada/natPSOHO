#ifndef Rcpp_head
#define Rcpp_head
#include <Rcpp.h>
using namespace Rcpp;
#endif

#include "utils.h"
#include "causality_list.h"

#ifndef natps_op
#define natps_op
Rcpp::NumericVector create_natcauslist_cpp(Rcpp::NumericVector &cl, Rcpp::List &net, StringVector &ordering);
Rcpp::CharacterMatrix cl_to_arc_matrix_cpp(const Rcpp::NumericVector &cl, Rcpp::CharacterVector &ordering, unsigned int rows);
int nat_pos_plus_vel_cpp(Rcpp::NumericVector &cl, const Rcpp::NumericVector &vl, const Rcpp::NumericVector &vl_neg, int n_arcs);
#endif

