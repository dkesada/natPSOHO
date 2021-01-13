#include "include/utils.h"

//' One-hot encoder for natural numbers without the 0
//' 
//' Given a natural number, return the natural number equivalent to its
//' one-hot encoding. Instead of pow, the '<<' operator will be used.
//' Examples: 3 -> 100 -> 4, 5 -> 10000 -> 16
//' @param nat the natural number to convert
//' @return the converted number
// [[Rcpp::export]]
int one_hot_cpp(int nat){
  return(1 << (nat - 1));
}

// Bitcount implementation from the book 'Hacker's Delight'
// Basically, a divide and conquer algorithm that sums the number of bits
// in two halves. It is not done recursively because the size of integers is
// fixed to 2^5, and so only 5 "mask and add" steps are needed.
// A less efficient but more readable version would be:
//
// int bitcount(unsigned int x)
// {
//   x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
//   x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
//   x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
//   x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
//   x = (x & 0x0000FFFF) + ((x >> 16)& 0x0000FFFF);
//   return x;
// }
int bitcount(unsigned x){
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x + (x >> 4)) & 0x0F0F0F0F;
  x = x + (x >> 8);
  x = x + (x >> 16);
  return x & 0x0000003F;
}

// Return the name of the node and its time slice
// 
// @param node a string with the name of the node
// @return a list with the name of the node and an integer with the time slice that the node belongs to
Rcpp::StringVector find_name_and_index(std::string node){
  Rcpp::StringVector res (2);
  std::string delim = "_t_";
  size_t pos;
  
  pos = node.find(delim);
  res[0] = node.substr(0, pos);
  res[1] = node.substr(pos + delim.length());
  
  return res;
}

void include_arc(Rcpp::StringMatrix &res, const Rcpp::StringVector &ordering, int i, int j, int &k){
  std::string from, to;
  int from_idx, to_idx;
  
  to_idx = i / ordering.size();
  to = ordering[to_idx];
  to += "_t_0";
  from_idx = i % ordering.size();
  from = ordering[from_idx];
  from += "_t_" + std::to_string(j);
  
  res(k, 0) = from;
  res(k, 1) = to;
  k++;
}

// Find the position of the node in the ordering. The node should be findable
// in the ordering. If not, an out of bounds index is returned
// 
// @param ordering a list with the order of the variables in t_0
// @param node the desired node
// @return the position of the node in the ordering
int find_index(const Rcpp::StringVector &ordering, std::string node){
  int i = 0;
  bool found = false;
  std::string name;
  
  while(i < ordering.size() && !found){
    name = ordering[i];
    if(name.find(node) != std::string::npos)
      found = true;
    else
      i++;
  }
  
  return i;
}

// Add two directions whose value has to be in the set {-1,0,1}
// 
// @param d1 first direction
// @param d2 second direction
// @param n_arcs the number of arcs present in the resulting causal list
// @return the result of adding them
int add_dirs(int d1, int d2, int &n_arcs){
  int res = d1 + d2;
  
  if(res < 0)
    res = 0;
  else if(res > 1)
    res = 1;
  
  if(res > d1)
    n_arcs++;
  else if(res < d1)
    n_arcs--;
  
  return res;
}

// Add two directions vectors whose value has to be in the set {-1,0,1}
// 
// @param d1 first NumericVector direction
// @param d2 second NumericVector direction
// @param n_arcs the number of arcs present in the resulting causal list
// @return the result of adding them
Rcpp::NumericVector add_dirs_vec(const NumericVector &d1, const NumericVector &d2, int &n_arcs){
  Rcpp::NumericVector res (d1.size());
  
  for(int i = 0; i < d1.size(); i++){
    res[i] = add_dirs(d1[i], d2[i], n_arcs);
  }
  
  return res;
}

// Subtract two directions whose value has to be in the set {-1,0,1}
// 
// @param d1 first direction
// @param d2 second direction
// @param n_arcs the number of arcs operations present in the resulting Velocity
// @return the result of subtracting them
int subtract_dirs(int d1, int d2, int &n_abs){
  int res = d1 - d2;
  
  if(d1 != d2)
    n_abs++;
  
  return res;
}

// Subtract two directions vectors whose value has to be in the set {-1,0,1}
// 
// @param d1 first NumericVector direction
// @param d2 second NumericVector direction
// @param n_arcs the number of arcs operations present in the resulting Velocity
// @return the result of subtracting them
Rcpp::NumericVector subtract_dirs_vec(const NumericVector &d1, const NumericVector &d2, int &n_abs){
  Rcpp::NumericVector res (d1.size());
  
  for(int i = 0; i < d1.size(); i++){
    res[i] = subtract_dirs(d1[i], d2[i], n_abs);
  }
  
  return res;
}

// Add two velocity directions whose value has to be in the set {-1,0,1}
// 
// @param d1 first direction
// @param d2 second direction
// @param abs_op the number of {1,-1} operations present in the resulting Velocity
// @return the result of adding them
int add_vel_dirs(int d1, int d2, int &abs_op){
  int res = d1 + d2;
  
  if(res < -1)
    res = -1;
  else if(res > 1)
    res = 1;
  
  if(res > d1 && res == 1)
    abs_op++;
  else if(res > d1 && res == 0)
    abs_op--;
  else if(res < d1 && res == 0)
    abs_op--;
  else if(res < d1 && res == -1)
    abs_op++;
  
  return res;
}

// Subtract two directions vectors whose value has to be in the set {-1,0,1}
// 
// @param d1 first NumericVector direction
// @param d2 second NumericVector direction
// @param abs_op the number of {1,-1} operations present in the resulting Velocity
// @return the result of adding them
Rcpp::NumericVector add_vel_dirs_vec(const NumericVector &d1, const NumericVector &d2, int &abs_op){
  Rcpp::NumericVector res (d1.size());
  
  for(int i = 0; i < d1.size(); i++){
    res[i] = add_vel_dirs(d1[i], d2[i], abs_op);
  }
  
  return res;
}

// Find the positions that are available to operate in a velocity
// 
// When adding or removing arcs via cte * vel, this finds positions where
// bits can be added or removed
//
// @param cl the positive causal list
// @param cl_neg the negative causal list
// @param max_int if 0 will search for integers greater than 0. Will search for < max_int otherwise.
// @return a vector with the open positions
std::vector<int> find_open_positions(const Rcpp::NumericVector &cl, const Rcpp::NumericVector &cl_neg, int max_int){
  std::vector<int> res(cl.size());
  bool add = max_int;
  int pos, pos_neg, j = 0;
  
  for(int i = 0; i < cl.size(); i++){
    pos = cl[i];
    pos_neg = cl_neg[i];
    pos |= pos_neg; // Remainder: arcs are shared between the positive and negative cl, so the number is obtained as a combination of both
    
    if((add && pos < max_int) || (!add && pos > 0)){
      res[j] = i;
      j++;
    }
  }
  
  res.resize(j);
  
  return res;
}

// Find the bits that are set to 0 or 1 in an integer
// 
// This can also be done recursively by masking, but in most cases the size
// of the network shouldn't exceed size 8 or so, which means only 8 iterations
// at worst. By using divide and conquer, I could get there in O(log(n)), but
// it would be more of a hassle than an improvement on the average case.
//
// @param x the integer to process
// @param remove if true, will search for bits set to 1. Will search for 0s otherwise.
// @param max_int the maximum integer allowed, i.e., the one with all possible bits set to 1
// @return a NumericVector with the open bits
Rcpp::NumericVector find_open_bits(int x, bool remove, int max_int){
  if(!remove)
    x ^= max_int;
  Rcpp::NumericVector res(bitcount(x));
  int i = 0, pos = 1;
  
  while(x != 0){
    if(x % 2){
      res[i] = pos;
      i++;
    }
    x >>= 1;
    pos++;
  }
  
  return res;
}

// Bitwise operator to remove 1s in the position when a 1 is in that same bit
// in the velocity. The truth table is:
//        vel
//        0 1
//  pos 0 0 0
//      1 1 0
int bitwise_sub_vel(int pos, int vel){
  return pos & (~vel);
}

// Just a debug function to try out stuff
// [[Rcpp::export]]
int debug_cpp(){
  std::vector<Rcpp::NumericVector> res(3);
  NumericVector tmp, pos_idx;
  int pos;
  
  tmp = seq(0, 10);
  pos_idx = sample(tmp, 1, false);
  pos = pos_idx[0];
  
  return pos;
}