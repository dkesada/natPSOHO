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
int nat_vel_plus_vel_cpp(Rcpp::NumericVector &vl1, Rcpp::NumericVector &vl1_neg, 
                          Rcpp::NumericVector &vl2, Rcpp::NumericVector &vl2_neg, 
                          int abs_op1, int abs_op2);
void add_nat_vel(int &num1, int num2, int &abs_op);
int nat_cte_times_vel_cpp(float k, Rcpp::NumericVector &vl, int abs_op, int max_size);
#endif