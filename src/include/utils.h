#ifndef Rcpp_head
#define Rcpp_head
#include <Rcpp.h>
using namespace Rcpp;
#endif

#ifndef nat_utils_op
#define nat_utils_op

#include <regex>
#include <random>
#include <vector>
#include <string>

static const std::vector<unsigned int> MASKS = { // --ICO-Merge: delete
  0x1,
  0x3,
  0xF,
  0xFF,
  0xFFFF};

int one_hot_cpp(int nat);
int bitcount(unsigned x);
Rcpp::StringVector find_name_and_index(std::string node);
void include_arc(Rcpp::StringMatrix &res, const Rcpp::StringVector &ordering, int i, int j, int &k);
int find_index(const Rcpp::StringVector &ordering, std::string node);
std::vector<int> find_open_positions(const Rcpp::NumericVector &cl, const Rcpp::NumericVector &cl_neg, int max_int);
Rcpp::NumericVector find_open_bits(int x, bool remove, int max_int);
int bitwise_sub(int x1, int x2);
Rcpp::List init_list_cpp(const Rcpp::Function &new_part, int n_inds, const Rcpp::StringVector &nodes, const Rcpp::StringVector &ordering, const Rcpp::StringVector &ordering_raw, int max_size, const Rcpp::NumericVector &v_probs, float p);
std::vector<int> find_open_bits_log(int x, bool remove, int max_int); // --ICO-Merge: delete
std::vector<int> find_open_bits_log_rec(int x, int idx, int depth); // --ICO-Merge: delete
Rcpp::NumericVector init_cl_cpp(int n_nodes);
Rcpp::StringVector crop_names_cpp(Rcpp::StringVector names);
int debug_cpp(int x, bool op, bool remove, int max_int);

#endif