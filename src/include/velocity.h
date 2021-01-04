#ifndef Rcpp_head
#define Rcpp_head
#include <Rcpp.h>
using namespace Rcpp;
#endif

#include "utils.h"
#include "causality_list.h"

#ifndef vl_op
#define vl_op
Rcpp::List randomize_vl_cpp(Rcpp::List &vl, NumericVector &probs, int seed);
Rcpp::List pos_minus_pos_cpp(Rcpp::List &cl, Rcpp::List &ps, Rcpp::List &vl);
void nat_vel_plus_vel_cpp(Rcpp::NumericVector &vl1, Rcpp::NumericVector &vl1_neg, Rcpp::NumericVector &vl2, Rcpp::NumericVector &vl2_neg, int &abs_op, int abs_op2);
Rcpp::List cte_times_vel_cpp(float k, Rcpp::List vl, unsigned int abs_op, unsigned int max_op);
#endif