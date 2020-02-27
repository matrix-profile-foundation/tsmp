#ifndef __MATH__
#define __MATH__

#include <Rcpp.h>

using namespace Rcpp;

#define MIN(y,x) ((x)<(y) && (x)==(x) ? (x) : (y))
#define MAX(y,x) ((x)>(y) && (x)==(x) ? (x) : (y))

double        std_rcpp(const NumericVector data, const bool na_rm);
NumericMatrix list_to_matrix(const List x); // unnused?
NumericVector diff_lag(const NumericVector x, const uint32_t lag); // unnused?
NumericVector diff2_lag(const NumericVector x, const uint32_t lag, const double v);
NumericVector fast_movsd_rcpp(const NumericVector data, const uint32_t window_size);
List          fast_avg_sd_rcpp(const NumericVector data, const uint32_t window_size);
int32_t       mode_rcpp(const NumericVector x);
NumericVector znorm_rcpp(const NumericVector data);
NumericVector binary_split_rcpp(const uint32_t n);
double        inner_product(NumericVector a, NumericVector b);
double        sum_of_squares(NumericVector a);
NumericVector sum2s_rcpp(NumericVector a, uint32_t w);
List          muinvn_rcpp(NumericVector a, uint32_t w);

#endif // __MATH__
