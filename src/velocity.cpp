#include "include/velocity.h"

//' Substracts two natPositions to obtain the natVelocity that transforms ps1 into ps2
//' 
//' @param ps1 the first position's causal list
//' @param ps2 the second position's causal list
//' @param vl the natVelocity's positive causal list
//' @param vl_neg the natVelocity's negative causal list
//' @return the velocity's causal lists by reference and the number of operations by return
// [[Rcpp::export]]
int nat_pos_minus_pos_cpp(const Rcpp::NumericVector &ps1, const Rcpp::NumericVector &ps2, Rcpp::NumericVector &vl, Rcpp::NumericVector &vl_neg){
  int ps1_i, ps2_i, vl_i, vl_neg_i;
  int n_abs = 0;
  
  for(int i = 0; i < ps1.size(); i++){
    ps1_i = ps1[i];
    ps2_i = ps2[i];
    
    vl_i = bitwise_sub(ps2_i, ps1_i);
    vl_neg_i = bitwise_sub(ps1_i, ps2_i);
    
    vl[i] = vl_i;
    vl_neg[i] = vl_neg_i;
    n_abs += bitcount(vl_i) + bitcount(vl_neg_i);
  }
  
  return n_abs;
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
                          const Rcpp::NumericVector &vl2, const Rcpp::NumericVector &vl2_neg, 
                          int abs_op1, int abs_op2){
  int pos1, pos2, neg1, neg2, mask, res;
  
  res = abs_op1 + abs_op2;
  for(int i = 0; i < vl1.size(); i++){
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
//' @param vl the Velocity's positive causal list
//' @param vl_neg the Velocity's negative causal list
//' @param abs_op the final number of {1,-1} operations
//' @param max_size the maximum size of the network
//' @return the new total number of operations 
// [[Rcpp::export]]
int nat_cte_times_vel_cpp(float k, Rcpp::NumericVector &vl, Rcpp::NumericVector &vl_neg, int abs_op, int max_size){
  int res, max_op, n_op, pos, pos_neg, pos_mix, pool_idx, pos_idx, bit_idx, bit_dest, max_int;
  bool remove;
  std::vector<int> pool;
  Rcpp::NumericVector pool_samp, bit_pool, bit_samp;
  
  max_int = one_hot_cpp(max_size) - 1;
  max_op = (max_size - 1) * vl.size();
  
  n_op = floor(k * abs_op);
  if(n_op > max_op)
    n_op = max_op;
  res = n_op;
  
  n_op = abs_op - n_op;
  remove = n_op > 0; // Whether to add or remove arcs
  n_op = std::abs(n_op);
  
  // Find a pool of possible integers in the cl and cl_neg to operate
  if(remove)
    pool = find_open_positions(vl, vl_neg, 0);
  else
    pool = find_open_positions(vl, vl_neg, max_int);
  
  for(int i = 0; i < n_op; i++){
    // Sample a position from the pool
    pool_samp = seq(0, pool.size() - 1);
    pool_samp = sample(pool_samp, 1, false);
    pool_idx = pool_samp[0];
    pos_idx = pool[pool_idx];
    
    // Find the pool of bits in that position
    pos = vl[pos_idx];
    pos_neg = vl_neg[pos_idx];
    pos_mix = pos | pos_neg;
    bit_pool = find_open_bits(pos_mix, remove, max_int);
    
    // Sample a bit and add it or remove it
    bit_samp = seq(0, bit_pool.size() - 1); // Sample the selected bit
    bit_samp = sample(bit_samp, 1, false);
    bit_idx = bit_samp[0];
    bit_idx = bit_pool[bit_idx];
    
    if(remove){
      if(pos & one_hot_cpp(bit_idx))
        pos ^= one_hot_cpp(bit_idx);
      else
        pos_neg ^= one_hot_cpp(bit_idx);
      pos_mix = pos | pos_neg;
      if(pos_mix == 0)
        pool.erase(pool.begin() + pool_idx);
    }
    
    else{
      bit_samp = seq(0, 1); // Sample whether to add the bit in the positive or negative cl
      bit_samp = sample(bit_samp, 1, false);
      bit_dest = bit_samp[0];
      if(bit_dest)
        pos_neg |= one_hot_cpp(bit_idx);
      else
        pos |= one_hot_cpp(bit_idx);
      pos_mix = pos | pos_neg;
      if(pos_mix == max_int)
        pool.erase(pool.begin() + pool_idx);
    }
    
    // Save the modification of the velocity
    vl[pos_idx] = pos;
    vl_neg[pos_idx] = pos_neg;
  }
  
  return res;
}

