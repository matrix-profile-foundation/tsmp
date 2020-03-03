/*===========================================================================*/
/* Original source from caTools library:                                     */
/* runfunc - running window functions                                        */
/* Copyright (C) 2005 Jarek Tuszynski                                        */
/* Distributed under GNU General Public License version 3                    */
/* Converted to Rcpp by Francisco Bischoff                                   */
/*===========================================================================*/

// Supports:
//               NA/NaN  -Inf/Inf  Edge
// movmin          Unk     Unk      No
// movmax          Unk     Unk      No

#include "math.h" // math first to fix OSX error
#include "windowfunc.h"

//[[Rcpp::export]]
NumericVector movmin(NumericVector data, uint32_t window_size) {
  uint32_t data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  uint32_t i, j, d, k;
  double min_res, res_out;
  NumericVector out(data_size - window_size + 1);

  k = 0;
  d = 0;

  min_res = res_out = R_PosInf;
  for (i = window_size - 1; i < data_size; i++) {
    // if point comining out of the window was window's min than we need to recalculate 'min_res'
    if (res_out == min_res) {
      min_res = R_PosInf;
      // find minimum over a window of length 'window_size'
      for (j = 0; j < window_size; j++) {
        min_res = MIN(min_res, data[d + j]);
      }
    } else {
      // if point comining out of the window was NOT window min than min of window's first
      //  'window_size - 1' points is still 'min_res', so we have to add a single point
      min_res = MIN(min_res, data[d + window_size - 1]);
    }

    res_out = data[d++]; // store point comming out of the window for future use and move window
    out[k++] = (min_res == R_PosInf ? R_NaReal : min_res); // save 'min_res' and move window
  }

  return out;
}

//[[Rcpp::export]]
NumericVector movmax(NumericVector data, uint32_t window_size) {
  uint32_t data_size = data.length();

  if (window_size <= 1) {
    return data;
  }

  if (window_size > data_size) {
    window_size = data_size;
  }

  uint32_t i, j, d, k;
  double max_res, res_out;
  NumericVector out(data_size - window_size + 1);

  k = 0;
  d = 0;

  max_res = res_out = R_NegInf;
  for (i = window_size - 1; i < data_size; i++) {
    // if point comining out of the window was window's max than we need to recalculate 'max_res'
    if (res_out == max_res) {
      max_res = R_NegInf;
      // find maximum over a window of length 'window_size'
      for (j = 0; j < window_size; j++) {
        max_res = MAX(max_res, data[d + j]);
      }
    } else {
      // if point comining out of the window was NOT window max than max of window's first
      //  'window_size - 1' points is still 'max_res', so we have to add a single point
      max_res = MAX(max_res, data[d + window_size - 1]);
    }

    res_out = data[d++]; // store point comming out of the window for future use and move window
    out[k++] = (max_res == R_NegInf ? R_NaReal : max_res); // save 'max_res' and move window
  }

  return out;
}
