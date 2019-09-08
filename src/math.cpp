#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

//[[Rcpp::export]]
double std_rcpp(NumericVector data, bool na_rm = false) {

  // if there is NaN in vector the result will be NaN
  if(any(is_na(data))) {
    if(na_rm) {
      data = na_omit(data);
    } else {
      return NA_REAL;
    }
  }

  double result = sqrt(sum(pow((data - mean(data)), 2)) / data.length());

  return(result);
}

//[[Rcpp::export]]
int32_t mode_rcpp(NumericVector x) {

  // is slower than R implementation...
  NumericVector ux = unique(x);
  int32_t y = ux[which_max(table(match(x, ux)))];
  return y;
}

//[[Rcpp::export]]
NumericVector znorm_rcpp(NumericVector data) {
  double data_mean = mean(data);
  double data_dev = sqrt(sum(pow((data - data_mean), 2)) / data.length());

  if (data_dev == NA_REAL || data_dev <= 0.01) {
    return(data - data_mean);
  }
  else {
    return((data - data_mean) / (data_dev));
  }
}
