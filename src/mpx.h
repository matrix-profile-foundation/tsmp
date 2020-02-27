#ifndef __MPX__
#define __MPX__

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

List mpx_rcpp(NumericVector a, uint16_t w, uint16_t minlag, bool idxs, bool euclidean);
List mpxab_rcpp(NumericVector a, NumericVector b, uint16_t w, bool idxs, bool euclidean);
List mpx_rcpp_parallel(NumericVector a, uint16_t w, uint16_t minlag, bool idxs, bool euclidean);
List mpxab_rcpp_parallel(NumericVector a, uint16_t w, uint16_t minlag, bool idxs, bool euclidean);

#endif // __MPX__
