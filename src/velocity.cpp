#include "include/velocity.h"

//' Substracts two Positions to obtain the Velocity that transforms one into the other
//' 
//' @param cl the first position's causal list
//' @param ps the second position's causal list
//' @param vl the Velocity's causal list
//' @return a list with the Velocity's causal list and the number of operations
// [[Rcpp::export]]
Rcpp::List pos_minus_pos_cpp(Rcpp::List &cl, Rcpp::List &ps, Rcpp::List &vl){
  Rcpp::List slice_cl, slice_ps, slice_vl, cu_cl, cu_ps, cu_vl, pair_cl, pair_ps, pair_vl;
  Rcpp::NumericVector dirs_cl, dirs_ps, dirs_vl;
  int n_abs = 0;
  Rcpp::List res (2);
  
  for(unsigned int i = 0; i < cl.size(); i++){
    slice_cl = cl[i];
    slice_ps = ps[i];
    slice_vl = vl[i];
    
    for(unsigned int j = 0; j < slice_cl.size(); j++){
      pair_cl = slice_cl[j];
      pair_ps = slice_ps[j];
      pair_vl = slice_vl[j];
      dirs_cl = pair_cl[1];
      dirs_ps = pair_ps[1];
      dirs_vl = subtract_dirs_vec(dirs_cl, dirs_ps, n_abs);
      
      pair_vl[1] = dirs_vl;
      slice_vl[j] = pair_vl;
    }
    
    vl[i] = slice_vl;
  }
  
  res[0] = vl;
  res[1] = n_abs;
  
  return res;
}

//' Adds two natVelocities 
//' 
//' Adds two natVelocities represented as two numeric vectors: one with the
//' positive part and one with the negative part. Adding them is a process that
//' does a bitwise 'or' with both the positive and negative parts of the two
//' velocities, adjusts the new abs_op, removes duplicated arcs in the final
//' velocity by using a bitwise 'xor' with both parts and adjusts the final abs_op.
//' The results are returned via modifying the original vl1 and vl1_neg by
//' reference and returning the final abs_op normally. I can't have an integer
//' edited by reference because it automatically gets casted and cannot be used
//' to return values. 
//' 
//' @param vl1 the first Velocity's positive part 
//' @param vl1_neg the first Velocity's negative part 
//' @param vl2 the second Velocity's positive part
//' @param vl2_neg the first Velocity's negative part 
//' @param abs_op1 the number of {1,-1} operations in the first velocity
//' @param abs_op2 the number of {1,-1} operations in the second velocity
//' @return the total number of resulting operations
// [[Rcpp::export]]
int nat_vel_plus_vel_cpp(Rcpp::NumericVector &vl1, Rcpp::NumericVector &vl1_neg,
                          Rcpp::NumericVector &vl2, Rcpp::NumericVector &vl2_neg, 
                          int abs_op1, int abs_op2){
  int pos1, pos2, neg1, neg2, mask, res;
  
  res = abs_op1 + abs_op2;
  for(unsigned int i = 0; i < vl1.size(); i++){
    pos1 = vl1[i];
    pos2 = vl2[i];
    neg1 = vl1_neg[i];
    neg2 = vl2_neg[i];
    
    add_nat_vel(pos1, pos2, res);
    add_nat_vel(neg1, neg2, res);
    mask = pos1 & neg1;
    
    if(mask){
      pos1 ^= mask;
      neg1 ^= mask;
      res -= 2 * bitcount(mask);
    }
    
    vl1[i] = pos1;
    vl1_neg[i] = neg1;
  }
  
  return res;
}

// Auxiliary function to make the main 'for' less verbose. It adds a natural
// number in two velocities by using an 'or' and removes duplicated operations 
// from abs_op.
void add_nat_vel(int &num1, int num2, int &abs_op){
  int mask = num1 & num2;
  if(mask)
    abs_op -= bitcount(mask);
  num1 |= num2;
}

//' Multiply a Velocity by a constant real number
//' 
//' @param k the constant real number
//' @param vl the Velocity's causal list
//' @param abs_op the final number of {1,-1} operations
//' @param max_op the maximum number of directions in the causal list
//' @return the new total number of operations 
// [[Rcpp::export]]
int cte_times_vel_cpp(float k, Rcpp::NumericVector &vl, int abs_op, int max_op){
  int res, n_op, idx, cmp;
  int l_pool = abs_op;
  bool invert = false;
  NumericVector tmp;
  
  // Process the k and max_op
  if(k < 0){
    k = fabs(k);
    invert = true;
  }
  
  n_op = floor(k * abs_op);
  res[1] = n_op;
    
  if(n_op > max_op){
    n_op = max_op;
    res = max_op;
  }
  
  n_op = abs_op - n_op;
  
  if(n_op < 0){ // Convert {0} into {1,-1}
    l_pool = max_op - abs_op; // Number of 0's remaining
    n_op = std::abs(n_op);
    pool = Rcpp::List(l_pool);
    cmp = 0;
  } 
  
  else{ // Convert {1,-1} into {0}
    n_op = std::abs(n_op);
    pool = Rcpp::List(l_pool);
    cmp = 1;
  }
  
  if(n_op > 0){
    // Loop through the cl to store the position and sign invert the 0's or the 1's depending on k greater or lesser than 0
    
    
    // Sample the position vector to position 0's or 1's in some or all of those positions
    pos = seq(0, (pool.size() - 1));
    pos = sample(pos, n_op, false);
    n_pool = Rcpp::List(n_op);
    for(unsigned int i = 0; i < pos.size(); i++){
      idx = pos[i];
      tmp = pool[idx];
      n_pool[i] = tmp;
    }
    
    // Operate the selected directions
  }
  
  res[0] = vl;
  
  return res;
}

