#include <Rcpp.h>
//#include <emmintrin.h>
using namespace Rcpp;

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
NumericMatrix list_to_matrix(List x){
  int32_t nlines = x.size();
  uint32_t colmax = 0;

  for(int32_t i = 0; i < nlines; i++) {
    uint32_t currsize = as<NumericVector>(x[i]).size();
    if(colmax < currsize) {
      colmax = currsize;
    }
  }

  NumericMatrix m(nlines, colmax);

  for(int32_t i = 0; i < nlines; i++) {
    int32_t line = nlines - i - 1;
    uint32_t currsize = as<NumericVector>(x[i]).size();
    NumericMatrix::Row row = m(i, _);
    row = as<NumericVector>(x[i]);

    for(uint32_t j = currsize; j < colmax; j++) {
      row[j] = 0;
    }
  }

  return(m);
}

//[[Rcpp::export]]
NumericVector diff_lag(NumericVector x, uint32_t lag = 1){
  uint32_t n = x.size();
  NumericVector out(n - lag);

  for(uint32_t i = 0; i < (n - lag); i++){
    out[i] = x[i + lag] - x[i];
  }
  return out;
}

//[[Rcpp::export]]
NumericVector diff2_lag(NumericVector x, uint32_t lag = 1, double v = 0.0){
  uint32_t n = x.size();
  NumericVector out(n - lag + 1);

  out[0] = v;

  for(uint32_t i = 0; i < (n - lag); i++){
    out[i + 1] = x[i + lag] - x[i];
  }
  return out;
}

//[[Rcpp::export]]
NumericVector fast_movsd_rcpp(NumericVector data, uint32_t window_size) {

  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  data = data - mean(data);

  NumericVector data_diff = diff_lag(data, window_size);
  data_diff.push_front(sum(data[Range(0,(window_size - 1))]));
  NumericVector data_sum = cumsum(data_diff);
  NumericVector data_mean = data_sum / window_size;

  NumericVector data2 = pow(data, 2);
  NumericVector data2_diff = diff_lag(data2, window_size);
  data2_diff.push_front(sum(data2[Range(0,(window_size - 1))]));
  NumericVector data2_sum = cumsum(data2_diff);
  NumericVector data_sd2 = (data2_sum / window_size) - pow(data_mean, 2); // variance
  NumericVector data_sd = sqrt(data_sd2);

  return (data_sd);
}

//[[Rcpp::export]]
NumericVector fast2_movsd_rcpp(NumericVector data, uint32_t window_size) {

  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  data = data - mean(data);

  NumericVector data_sum = cumsum(diff2_lag(data, window_size, sum(data[Range(0,(window_size - 1))])));
  NumericVector data_mean = data_sum / window_size;

  NumericVector data2 = pow(data, 2);
  NumericVector data2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0,(window_size - 1))])));
  NumericVector data_sd2 = (data2_sum / window_size) - pow(data_mean, 2); // variance
  NumericVector data_sd = sqrt(data_sd2);

  return (data_sd);
}

//[[Rcpp::export]]
List fast_avg_sd_rcpp(NumericVector data, uint32_t window_size) {

  NumericVector data_diff = diff_lag(data, window_size);
  data_diff.push_front(sum(data[Range(0,(window_size - 1))]));
  NumericVector mov_sum = cumsum(data_diff);
  NumericVector mov_mean = mov_sum / window_size;

  NumericVector data2 = pow(data, 2);
  NumericVector data2_diff = diff_lag(data2, window_size);
  data2_diff.push_front(sum(data2[Range(0,(window_size - 1))]));
  NumericVector mov2_sum = cumsum(data2_diff);


  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  NumericVector data_zeromean = data - mean(data);

  data_diff = diff_lag(data_zeromean, window_size);
  data_diff.push_front(sum(data_zeromean[Range(0,(window_size - 1))]));
  NumericVector data_sum = cumsum(data_diff);
  NumericVector data_mean = data_sum / window_size;

  data2 = pow(data_zeromean, 2);
  data2_diff = diff_lag(data2, window_size);
  data2_diff.push_front(sum(data2[Range(0,(window_size - 1))]));
  NumericVector data2_sum = cumsum(data2_diff);
  NumericVector data_sd2 = (data2_sum / window_size) - pow(data_mean, 2); // variance
  NumericVector data_sd = sqrt(data_sd2);
  NumericVector data_sig = sqrt(1/(data_sd2 * window_size));

  return (List::create(
      Rcpp::Named("avg") = mov_mean,
      Rcpp::Named("sd") = data_sd,
      Rcpp::Named("sig") = data_sig,
      Rcpp::Named("sum") = mov_sum,
      Rcpp::Named("sqrsum") = mov2_sum
  ));
}



// [[Rcpp::export]]
int vecmin(IntegerVector x) {
  // Rcpp supports STL-style iterators
  IntegerVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}


// [[Rcpp::export]]
int vecmax(IntegerVector x) {
  // Rcpp supports STL-style iterators
  IntegerVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
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

//[[Rcpp::export]]
NumericVector binary_split_rcpp(uint32_t n) {

  NumericVector idxs(n);

  idxs[0] = 1;// We always begin by explore the first integer
  // After exploring the first integer, we begin splitting the interval 2:n

  std::deque<uint32_t> lb_list;
  std::deque<uint32_t> ub_list;

  lb_list.push_back(2);
  ub_list.push_back(n);

  uint32_t lb;
  uint32_t ub;
  uint32_t mid;

  for(uint32_t i = 1; i < n; i++) {
    lb = lb_list.front();
    ub = ub_list.front();
    mid = (lb + ub) / 2; // integer division is automatically floor()
    lb_list.pop_front();
    ub_list.pop_front();

    idxs[i] = mid;

    if(lb == ub) {
      continue;
    } else {
      if(lb < mid) {
        lb_list.push_back(lb);
        ub_list.push_back(mid - 1);
      }

      if(ub > mid) {
        lb_list.push_back(mid + 1);
        ub_list.push_back(ub);
      }
    }
  }

  return(idxs);
}

//[[Rcpp::export]]
double sq2s_rcpp(NumericVector a) {
  NumericVector h = a * a;
  NumericVector c = (pow(2,27) + 1.0) * a; // <-- can be replaced with fma where available
  NumericVector a1 = (c - (c - a));
  NumericVector a2 = a - a1;
  NumericVector a3 = a1 * a2;
  NumericVector r = a2 * a2 - (((h - a1 * a1) - a3) - a3);

  double p = h[0];
  double s = r[0];

  for (uint32_t i = 1; i < a.length(); i++) {
    double x = p + h[i];
    double z = x - p;
    s = s + (((p - (x - z)) + (h[i] - z)) + r[i]);
    p = x;
  }
  double res = p + s;

  return(res);
}
