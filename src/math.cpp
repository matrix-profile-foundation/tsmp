#include "math.h"

//[[Rcpp::export]]
double std_rcpp(const NumericVector data, const bool na_rm = false) {

  NumericVector the_data = data;

  // if there is NaN in vector the result will be NaN
  if(any(is_na(data))) {
    if(na_rm) {
      the_data = na_omit(data);
    } else {
      return NA_REAL;
    }
  }

  double result = sqrt(sum(pow((the_data - mean(the_data)), 2)) / the_data.length());

  return(result);
}

//[[Rcpp::export]]
NumericMatrix list_to_matrix(const List x){
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
   // int32_t line = nlines - i - 1;
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
NumericVector diff_lag(const NumericVector x, const uint32_t lag = 1){
  uint32_t n = x.size();
  NumericVector out(n - lag);

  for(uint32_t i = 0; i < (n - lag); i++){
    out[i] = x[i + lag] - x[i];
  }
  return out;
}

//[[Rcpp::export]]
NumericVector diff2_lag(const NumericVector x, const uint32_t lag = 1, const double v = 0.0) {
  uint32_t n = x.size();
  NumericVector out(n - lag + 1);

  out[0] = v;

  for(uint32_t i = 0; i < (n - lag); i++){
    out[i + 1] = x[i + lag] - x[i];
  }
  return out;
}

//[[Rcpp::export]]
NumericVector fast_movsd_rcpp(const NumericVector data, const uint32_t window_size) {

  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  NumericVector data_zeromean = data - mean(data);

  NumericVector data_sum = cumsum(diff2_lag(data_zeromean, window_size, sum(data_zeromean[Range(0,(window_size - 1))])));
  NumericVector data_mean = data_sum / window_size;

  NumericVector data2 = pow(data_zeromean, 2);
  NumericVector data2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0,(window_size - 1))])));
  NumericVector data_sd2 = (data2_sum / window_size) - pow(data_mean, 2); // variance
  NumericVector data_sd = sqrt(data_sd2);

  return (data_sd);
}

//[[Rcpp::export]]
List fast_avg_sd_rcpp(const NumericVector data, const uint32_t window_size) {

  NumericVector mov_sum = cumsum(diff2_lag(data, window_size, sum(as<NumericVector>(data[Range(0,(window_size - 1))]))));
  NumericVector mov_mean = mov_sum / window_size;
  NumericVector data2 = pow(data, 2);
  NumericVector mov2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0,(window_size - 1))])));

  // Improve the numerical analysis by subtracting off the series mean
  // this has no effect on the standard deviation.
  NumericVector data_zeromean = data - mean(data);

  NumericVector data_sum = cumsum(diff2_lag(data_zeromean, window_size, sum(data_zeromean[Range(0,(window_size - 1))])));
  NumericVector data_mean = data_sum / window_size;

  data2 = pow(data_zeromean, 2);
  NumericVector data2_sum = cumsum(diff2_lag(data2, window_size, sum(data2[Range(0,(window_size - 1))])));
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

//[[Rcpp::export]]
int32_t mode_rcpp(const NumericVector x) {

  // is slower than R implementation...
  NumericVector ux = unique(x);
  int32_t y = ux[which_max(table(match(x, ux)))];
  return y;
}

//[[Rcpp::export]]
NumericVector znorm_rcpp(const NumericVector data) {
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
NumericVector binary_split_rcpp(const uint32_t n) {

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
double inner_product(NumericVector a, NumericVector b) {
  double res = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);

  return(res);
}

//[[Rcpp::export]]
double sum_of_squares(NumericVector a) {
  double res = std::inner_product(a.begin(), a.end(), a.begin(), 0.0);

  return(res);
}


NumericVector sum2s_rcpp(NumericVector a, uint32_t w) {
  NumericVector res(a.length() - w + 1, 0);
  double accum = a[0];
  double resid = 0.0;

  for (uint32_t i = 1; i < w; i++) {
    double m = a[i];
    double p = accum;
    accum = accum + m;
    double q = accum - p;
    resid = resid + ((p - (accum - q)) + (m - q));
  }

  res[0] = accum + resid;

  for (uint32_t i = w; i < a.length(); i++) {
    double m = a[i - w];
    double n = a[i];
    double p = accum - m;
    double q = p - accum;
    double r = resid + ((accum - (p - q)) - (m + q));
    accum = p + n;
    double t = accum - p;
    resid = r + ((p - (accum - t)) + (n - t));
    res[i - w + 1] = accum + resid;
  }

  return(res);
}

List muinvn_rcpp(NumericVector a, uint32_t w) {
  // Functions here are based on the work in
  // Ogita et al, Accurate Sum and Dot Product
  // results here are a moving average and stable inverse centered norm based
  // on Accurate Sum and Dot Product, Ogita et al

  NumericVector sig(a.length() - w + 1, 0);
  NumericVector mu = sum2s_rcpp(a, w) / w;

  for (uint32_t i = 0; i < mu.length(); i++) {
    sig[i] = sum_of_squares(a[Range(i, i + w - 1)] - mu[i]);
  }

  sig = 1 / sqrt(sig);

  return (List::create(
      Rcpp::Named("avg") = mu,
      Rcpp::Named("sig") = sig
  ));
}
