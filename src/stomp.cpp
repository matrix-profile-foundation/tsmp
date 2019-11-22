// #include <Rcpp.h>
// using namespace Rcpp;
//
// double inner_product(NumericVector a, NumericVector b) {
//   double res = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
//
//   return(res);
// }
//
// double sum_of_squares(NumericVector a) {
//   double res = std::inner_product(a.begin(), a.end(), a.begin(), 0.0);
//
//   return(res);
// }
//
// NumericVector sum2s_rcpp(NumericVector a, uint32_t w) {
//   NumericVector res(a.length() - w + 1, 0);
//   double accum = a[0];
//   double resid = 0.0;
//
//   for (uint32_t i = 1; i < w; i++) {
//     double m = a[i];
//     double p = accum;
//     accum = accum + m;
//     double q = accum - p;
//     resid = resid + ((p - (accum - q)) + (m - q));
//   }
//
//   res[0] = accum + resid;
//
//   for (uint32_t i = w; i < a.length(); i++) {
//     double m = a[i - w];
//     double n = a[i];
//     double p = accum - m;
//     double q = p - accum;
//     double r = resid + ((accum - (p - q)) - (m + q));
//     accum = p + n;
//     double t = accum - p;
//     resid = r + ((p - (accum - t)) + (n - t));
//     res[i - w + 1] = accum + resid;
//   }
//
//   return(res);
// }
//
// List muinvn_rcpp(NumericVector a, uint32_t w) {
//   // Functions here are based on the work in
//   // Ogita et al, Accurate Sum and Dot Product
//   // results here are a moving average and stable inverse centered norm based
//   // on Accurate Sum and Dot Product, Ogita et al
//
//   NumericVector sig(a.length() - w + 1, 0);
//   NumericVector mu = sum2s_rcpp(a, w) / w;
//
//   for (uint32_t i = 0; i < mu.length(); i++) {
//     sig[i] = sum_of_squares(a[Range(i, i + w - 1)] - mu[i]);
//   }
//
//   sig = 1 / sqrt(sig);
//
//   return (List::create(
//       Rcpp::Named("avg") = mu,
//       Rcpp::Named("sig") = sig
//   ));
// }
//
// //function(..., window_size, exclusion_zone = 1 / 2, verbose = 2)
//
// // [[Rcpp::export]]
// List stomp_rcpp(NumericVector a, uint16_t w, uint16_t minlag, bool indexes = false, uint32_t s_len = 0) {
//
//   double c, c_cmp;
//   uint32_t off_max, off_diag, offset;
//   // matrix profile using cross correlation,
//   uint32_t n = a.length();
//
//   List msd = muinvn_rcpp(a, w);
//
//   NumericVector mmu = msd["avg"];
//   NumericVector ssig = msd["sig"];
//   double* mu = &mmu[0];
//   double* sig = &ssig[0];
//
//   uint32_t diagmax = n - w + 1;
//   IntegerVector seq_diag = Range(minlag, diagmax - 1);
//
//   if(s_len == 0 || s_len > seq_diag.length()) {
//     s_len = seq_diag.length();
//   }
//
//   seq_diag = sample(seq_diag, s_len);
//
//   NumericVector mmp(diagmax, -1.0);
//   IntegerVector mmpi(diagmax, R_NaN);
//
//   double* mp = &mmp[0];
//   int* mpi = &mmpi[0];
//
//   // differentials have 0 as their first entry. This simplifies index
//   // calculations slightly and allows us to avoid special "first line"
//   // handling.
//
//   NumericVector ddf = 0.5 * (a[Range(w, n - 1)] - a[Range(0, n - w - 1)]);
//   ddf.push_front(0);
//   NumericVector ddg = (a[Range(w, n - 1)] - mmu[Range(1, diagmax - 1)]) + (a[Range(0, n - w - 1)] - mmu[Range(0, n - w - 1)]);
//   ddg.push_front(0);
//
//   double* df = &ddf[0];
//   double* dg = &ddg[0];
//
//   NumericVector bb = (a[Range(0, w - 1)] - mmu[0]);
//
//   for (IntegerVector::iterator diag = seq_diag.begin(); diag != seq_diag.end(); ++diag) {
//     c = inner_product((a[Range(*diag, (*diag + w - 1))] - mu[*diag]), bb);
//     off_max = (n - w - *diag + 1);
//     for (offset = 0; offset < off_max; offset++) {
//       off_diag = offset + *diag;
//       c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
//       c_cmp = c * sig[offset] * sig[off_diag];
//       if (c_cmp > mp[offset]) {
//         mp[offset] = c_cmp;
//         if(indexes) {
//           mpi[offset] = off_diag + 1;
//         }
//       }
//       if (c_cmp > mp[off_diag]) {
//         mp[off_diag] = c_cmp;
//         if(indexes) {
//           mpi[off_diag] = offset + 1;
//         }
//       }
//     }
//   }
//
//   // to do ed
//   mmp[mmp > 1.0] = 1.0;
//   mmp = sqrt(2 * w * (1 - mmp));
//
//   return (List::create(
//       Rcpp::Named("mp") = mmp,
//       Rcpp::Named("mpi") = mmpi
//   ));
// }
