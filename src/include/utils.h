#ifndef Rcpp_head
#define Rcpp_head
#include <Rcpp.h>
using namespace Rcpp;
#endif

#ifndef natutils_op
#define natutils_op

#include <regex>
#include <random>
#include <vector>

int one_hot_cpp(int nat);
int bitcount(unsigned x);
Rcpp::StringVector find_name_and_index(std::string node);
void include_arc(Rcpp::StringMatrix &res, const Rcpp::StringVector &ordering, int i, int j, int &k);
int find_index(const Rcpp::StringVector &ordering, std::string node);
std::vector<int> find_open_positions(const Rcpp::NumericVector &cl, const Rcpp::NumericVector &cl_neg, int max_int);
Rcpp::NumericVector find_open_bits(int x, bool remove, int max_int);
int bitwise_sub(int pos, int vel);
int debug_cpp();
#endif