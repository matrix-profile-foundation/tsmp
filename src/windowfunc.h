#ifndef __WINDOWFUNC__
#define __WINDOWFUNC__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector movmin(NumericVector data, uint32_t window_size);
NumericVector movmax(NumericVector data, uint32_t window_size);

#endif // __WINDOWFUNC__
