#include <Rcpp.h>
#include <RcppBlaze.h>
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

//[[Rcpp::export]]
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
// [[Rcpp::export]]
List muinvn_rcpp(NumericVector a, uint32_t w) {
  // Functions here are based on the work in
  // Ogita et al, Accurate Sum and Dot Product
  // results here are a moving average and stable inverse centered norm based
  // on Accurate Sum and Dot Product, Ogita et al

  NumericVector sig(a.length() - w + 1, 0);
  NumericVector mu = sum2s_rcpp(a, w) / w;

  for (uint32_t i = 0; i < mu.length(); i++) {
    sig[i] = sq2s_rcpp(a[Range(i, i + w - 1)] - mu[i]);
  }

  sig = 1 / sqrt(sig);

  return (List::create(
      Rcpp::Named("avg") = mu,
      Rcpp::Named("sig") = sig
  ));
}

//[[Rcpp::export]]
NumericVector sum2s2_rcpp(NumericVector a, uint32_t w) {
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

// [[Rcpp::export]]
List muinvn2_rcpp(NumericVector a, uint32_t w) {
  // Functions here are based on the work in
  // Ogita et al, Accurate Sum and Dot Product
  // results here are a moving average and stable inverse centered norm based
  // on Accurate Sum and Dot Product, Ogita et al

  NumericVector mu(a.length() - w + 1, 0);
  double accum = a[0];
  double resid = 0.0;

  for (uint32_t i = 1; i < w; i++) {
    double m = a[i];
    double p = accum;
    accum = accum + m;
    double q = accum - p;
    resid = resid + ((p - (accum - q)) + (m - q));
  }

  mu[0] = (accum + resid) / w;

  for (uint32_t i = w; i < a.length(); i++) {
    double m = a[i - w];
    double n = a[i];
    double p = accum - m;
    double q = p - accum;
    double r = resid + ((accum - (p - q)) - (m + q));
    accum = p + n;
    double t = accum - p;
    resid = r + ((p - (accum - t)) + (n - t));
    mu[i - w + 1] = (accum + resid) / w;
  }

  NumericVector sig(a.length() - w + 1, 0);
  NumericVector aa(w);
  NumericVector h(w);
  NumericVector c(w);
  NumericVector a1(w);
  NumericVector a2(w);
  NumericVector a3(w);
  NumericVector r(w);

  for (uint32_t i = 0; i < mu.length(); i++) {
    for(uint32_t j = 0; j < w; j++) {
      aa[j] = a[i+j] - mu[i];
      h[j] = aa[j] * aa[j];
      c[j] = (pow(2,27) + 1.0) * aa[j];
      a1[j] = (c[j] - (c[j] - aa[j]));
      a2[j] = aa[j] - a1[j];
      a3[j] = a1[j] * a2[j];
      r[j] = a2[j] * a2[j] - (((h[j] - a1[j] * a1[j]) - a3[j]) - a3[j]);
    }

    double p = h[0];
    double s = r[0];

    for (uint32_t j = 1; j < w; j++) {
      double x = p + h[j];
      double z = x - p;
      s = s + (((p - (x - z)) + (h[j] - z)) + r[j]);
      p = x;
    }

    sig[i] = p + s;
  }

  sig = 1 / sqrt(sig);

  return (List::create(
      Rcpp::Named("avg") = mu,
      Rcpp::Named("sig") = sig
  ));
}

// [[Rcpp::export]]
List mpx_rcpp(NumericVector a, uint32_t w, uint32_t minlag) {

  // matrix profile using cross correlation,
  uint32_t n = a.length();

  Environment pkg = Environment::namespace_env("tsmp");
  Function fast_avg_sd = pkg["fast_avg_sd"];

  List msd = fast_avg_sd(a, w);

  NumericVector mu = msd["avg"];
  NumericVector sig = msd["sig"];
  // differentials have 0 as their first entry. This simplifies index
  // calculations slightly and allows us to avoid special "first line"
  // handling.
  //
  uint32_t diagmax = n - w + 1;

  NumericVector df = 0.5 * (a[Range(w, n - 1)] - a[Range(0, n - w - 1)]);
  df.push_front(0);
  NumericVector dg = (a[Range(w, n - 1)] - mu[Range(1, diagmax - 1)]) + (a[Range(0, n - w - 1)] - mu[Range(0, n - w - 1)]);
  dg.push_front(0);
  NumericVector mp = rep(-1.0, diagmax);
  NumericVector mpi = rep(R_NaN, diagmax);

  for (uint32_t diag = minlag + 1; diag <= diagmax; diag++) {
    double c = sum((a[Range(diag - 1, (diag - 1 + w - 1))] - mu[diag - 1]) * (a[Range(0, w - 1)] - mu[0]));
    for (uint32_t offset = 1; offset <= (n - w - diag + 2); offset++) {
      c = c + df[offset - 1] * dg[offset - 1 + diag - 1] + df[offset - 1 + diag - 1] * dg[offset - 1];
      double c_cmp = c * sig[offset - 1] * sig[offset - 1 + diag - 1];
      if (c_cmp > mp[offset - 1]) {
        mp[offset - 1] = c_cmp;
        // mpi[offset - 1] = offset - 1 + diag - 1;
        mpi[offset - 1] = offset - 1 + diag;
      }
      if (c_cmp > mp[offset - 1 + diag - 1]) {
        mp[offset - 1 + diag - 1] = c_cmp;
        // mpi[offset - 1 + diag - 1] = offset - 1;
        mpi[offset - 1 + diag - 1] = offset;
      }
    }
  }
  // to do ed
  mp[mp > 1.0] = 1.0;
  // mp = sqrt(2 * w * (1 - mp));


  return (List::create(
      Rcpp::Named("mp") = mp,
      Rcpp::Named("mpi") = mpi
  ));
}
