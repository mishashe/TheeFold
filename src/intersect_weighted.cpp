#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <random>
#include <chrono>
#include <math.h>
#include <stdint.h>
#define ARMA_WARN_LEVEL 0
using namespace Rcpp;
using namespace arma;
using namespace std;

//' function to get ChiSquare
//'
//' @param par vector of breaks
//' @param muts mutations vector
//' @param L genome length
//' @return ChiSquare
//' @export
// [[Rcpp::export]]
arma::mat intersect_weighted(arma::mat intervals1, arma::mat intervals2)
{
  std::vector<long long int> starts;
  std::vector<long long int> ends;
  std::vector<long long int> weights;
  
  
  for (long long int i1 = 0; i1 < intervals1.n_rows; i1++)
  {
    long long int x1 = intervals1(i1,0);
    long long int y1 = intervals1(i1,1);
    for (long long int i2 = 0; i2 < intervals2.n_rows; i2++)
    {
      long long int x2 = intervals2(i2,0);
      long long int y2 = intervals2(i2,1);
      long long int Min = std::max(x1, x2);
      long long int Max = std::min(y1, y2);
      
      if (Min < Max) {
        starts.push_back(Min);
        ends.push_back(Max);
        weights.push_back(intervals1(i1,2)*intervals2(i2,2));
      }
    }
  }
  long long int L = weights.size();

  arma::mat out(L, 3);
  for(int i = 0; i < L; i++){
    out(i,0) = starts[i];
    out(i,1) = ends[i];
    out(i,2) = weights[i];
  }
  return(out);
}





