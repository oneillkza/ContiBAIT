//Greedy ordering based on contig quality.
#include <R_ext/Print.h>
#include <Rcpp.h>
#include <boost/dynamic_bitset.hpp>
#include <stdio.h>
#include <stdlib.h>
#define  STD_COUT
#include "orderContigsGreedy.h" 
#include <map>
#include <vector>
#include <math.h>
using namespace Rcpp ;
using namespace boost;




int nRows;
int nCols;
IntegerMatrix strandMatrix;


float score_order(std::vector<int> order_vector)
{

  
  int numSCE = 0;
  int prev_states[nCols];
  for (int i = 0; i < nCols; i++){
    prev_states[i] = 0;
  } 
  for (int i = 0; i < (int)order_vector.size(); i++)
  {
    IntegerVector states = strandMatrix.row(order_vector[i]);
    for (int j = 0; j < nCols; j++)
    {
      if (states[j] == prev_states[j] || states[j] == 0)
      {
        continue;
      }
      if ((prev_states[j] == 1 && states[j] == 3) || (prev_states[j] == 3 && states[j] == 1 )){
        numSCE += 5;
      }
      if (prev_states[j] != 0)
      {
        numSCE++;
      }
      prev_states[j] = states[j];        
    }
   }
  return int ((nRows*nCols) - numSCE);
};


  
/*****************************************************************************
* Function to greedily order contigs based on their strand state using contig quality 
* 
****************************************************************************/

//@Param SStrandMatrix  Strand state for all contigs across all libraries
RcppExport  SEXP orderContigsGreedy(
                SEXP SStrandMatrix                
                )
{
  IntegerMatrix sMatrix(SStrandMatrix);
  strandMatrix = sMatrix; 
  nRows = strandMatrix.nrow();
  nCols = strandMatrix.ncol();
  std::vector<int> best_order;
  best_order.push_back(0);
  best_order.push_back(1);
  for(int i = 2; i < nRows; i++){
    int best_index = 0;
    int best_index_score = 0;
    for(int j = 0; j <= (int)best_order.size(); j++){
      std::vector<int> order;
      for(int k = 0; k < (int)best_order.size(); k++){
        order.push_back(best_order[k]);
      }
    if (j == (int)order.size()){
      order.push_back(j);
    } else{
      std::vector<int>::iterator it2 = order.begin();
      order.insert(it2 + j, i);
    }
      int temp_score = score_order(order);
      if (temp_score > best_index_score){
        best_index = j;
        best_index_score = temp_score;
      }      
    }
    if (best_index == (int)best_order.size()){
      best_order.push_back(i);
    }else{
    std::vector<int>::iterator it = best_order.begin();
    best_order.insert(it + best_index, i);
    }
  }
  IntegerVector result(nRows);
  for (int i = 0; i < nRows; i++){
    result[i] = best_order[i] + 1;
  }
 return Rcpp::List::create(Rcpp::Named("order") = result,
                       Rcpp::Named("score") = score_order(best_order));
    
}

