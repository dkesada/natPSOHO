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
// [[Rcpp::export]]
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
// at worst. By using divide and conquer, I could get there in O(log(n)). The
// implementation I did is slower than this, probably because of appending
// the vectors of results in the recursion and because it ends up doing more
// operations by masking and shifting.
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

// Binary bitwise operator that returns 1 only when the first bit is 1 and the
// second one is 0. Kinda like subtracting 1 bit only when there is a 1 in the 
// other bit. It helps removing 1s in the position when a 1 is in that same bit
// in the velocity, or finding the velocity that takes you from a position to another.
// The truth table is:
//        x2
//       0 1
//  x1 0 0 0
//     1 1 0
int bitwise_sub(int x1, int x2){
  return x1 & (~x2);
}

//' Initialize the particles
//' 
//' Initialize the list with particles in C++. It is equivalent to initializing
//' them in R, so this will be dropped. --ICO-Merge: delete if obsolete
//' @param new_part function that creates a new particle
//' @param n_inds number of particles that the algorithm will simultaneously process
//' @param nodes a vector with the names of the nodes
//' @param ordering the names of the nodes
//' @param ordering_raw 
//' @param max_size the maximum size of the DBN
//' @param v_probs vector that defines the random velocity initialization probabilities
//' @param p parameter of the truncated geometric distribution for sampling edges
//' @return a list with the randomly initialized particles
// [[Rcpp::export]]
Rcpp::List init_list_cpp(const Rcpp::Function &new_part, int n_inds, const Rcpp::StringVector &nodes, const Rcpp::StringVector &ordering, 
                         const Rcpp::StringVector &ordering_raw, int max_size, const Rcpp::NumericVector &v_probs, 
                         float p){
  Rcpp::List res (n_inds);
  
  for(int i = 0; i < n_inds; i++){
    Environment part;
    part = new_part(nodes, ordering, ordering_raw, max_size, v_probs, p);
    res[i] = part;
  }
  
  return res;
}

// Find the bits that are set to 0 or 1 in an integer using divide and conquer
// 
// Recursive algorithm that finds the position of 1s or 0s by masking and 
// shifting bits. In most cases, it has a O(log(n)) cost, in the worst case 
// scenario of 0xFFFFFFFF, it has O(n) cost. Cool idea in theory, but it ends
// up being slower than the trivial way. Too many vectors appended, probably.
// --ICO-Merge: worse than the trivial approach. Delete 
//
// @param x the integer to process
// @param remove if true, will search for bits set to 1. Will search for 0s otherwise.
// @param max_int the maximum integer allowed, i.e., the one with all possible bits set to 1
// @return a NumericVector with the open bits
std::vector<int> find_open_bits_log(int x, bool remove, int max_int){
  if(!remove)
    x ^= max_int;
  return find_open_bits_log_rec(x, 1, 4);
}

// Recursive part of the algorithm that finds the position of 1s or 0s. Given a 
// number, if 
//
// --ICO-Merge: delete
// @param x the integer to process
// @param remove if true, will search for bits set to 1. Will search for 0s otherwise.
// @param max_int the maximum integer allowed, i.e., the one with all possible bits set to 1
// @return a NumericVector with the open bits
std::vector<int> find_open_bits_log_rec(int x, int idx, int depth){
  std::vector<int> res;
  
  if(x == 1)
    res.push_back(idx);
  
  else if (x > 1){
    int x1, x2, shift;
    shift = 1 << depth;
    x1 = x & MASKS[depth];
    x2 = x >> shift;
    
    res = find_open_bits_log_rec(x1, idx, depth-1);
    std::vector<int> res2 = find_open_bits_log_rec(x2, idx + shift, depth-1);
    if(res.size() > res2.size())
      res.insert(res.end(), res2.begin(), res2.end());
    else{
      res2.insert(res2.end(), res.begin(), res.end());
      res = res2;
    }
  }
  
  return res;
}

//' Initialize the nodes vector
//' 
//' Initialize the vector in C++
//' @param n_nodes number of receiving nodes
//' @return a list with the randomly initialized particles
// [[Rcpp::export]]
Rcpp::NumericVector init_cl_cpp(int n_nodes){
  Rcpp::NumericVector res (n_nodes);
  
  return res;
}

//' If the names of the nodes have "_t_0" appended at the end, remove it
//' @param names a vector with the names of the nodes in t_0
//' @return the vector with the names cropped
// [[Rcpp::export]]
Rcpp::StringVector crop_names_cpp(const Rcpp::StringVector &names){
  StringVector res = Rcpp::clone(names);
  std::string tmp;
  std::size_t pos;
  
  for(int i = 0; i < res.size(); i++){
    tmp = res[i];
    pos = tmp.find("_t_0");
    tmp = tmp.substr(0, pos);
    res[i] = tmp;
  }
  
  return(res);
}

// Just a debug function to try out stuff
// [[Rcpp::export]]
int debug_cpp(int x, bool op, bool remove, int max_int){
  if(op)
    find_open_bits(x, remove, max_int);
  else
    find_open_bits_log(x, remove, max_int);
  
  return 0;
}