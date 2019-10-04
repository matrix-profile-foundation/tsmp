#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//


// NumericVector mpx_v_rcpp(NumericVector a, uint32_t w, uint32_t minlag) {
//
//   // matrix profile using cross correlation,
//   uint32_t n = a.length();
//
//   // msd = fast_avg_sd(a, w);
//   NumericVector mu; //= msd$avg;
//   NumericVector sig; //= msd$sig;
//
//   // differentials have 0 as their first entry. This simplifies index
//   // calculations slightly and allows us to avoid special "first line"
//   // handling.
//   uint32_t diagmax = n - w + 1;
//   NumericVector df = c(0, (1 / 2) * (a[(1 + w):n] - a[1:(n - w)]));
//   NumericVector dg = c(0, (a[(1 + w):n] - mu[2:(diagmax)]) + (a[1:(n - w)] - mu[1:(n - w)]));
// NumericVector mp = rep(-1, diagmax);
//   NumericVector mpi = rep(NA, diagmax);
//
//   NumericVector seq_diag = (minlag + 1):diagmax;
//   seq_diag = sample(seq_diag, size = seq_diag.length());
//
//   for (diag in seq_diag) {
//     double c = sum((a[diag:(diag + w - 1)] - mu[diag]) * (a[1:w] - mu[1]));
//
//     NumericVector offset = 1:(n - w - diag + 2);
//     NumericVector off_diag = (offset + diag - 1);
//     NumericVector d = df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
//
//     d[1] = d[1] + c;
//     d = cumsum(d);
//
//     NumericVector d_cmp = d * sig[offset] * sig[off_diag];
//
//     mask = d_cmp > mp[offset];
//     mp[c(mask, rep(FALSE, diag - 1))] = d_cmp[mask];
//     mpi[c(mask, rep(FALSE, diag - 1))] = off_diag[mask];
//
//     mask = d_cmp > mp[off_diag];
//     mp[c(rep(FALSE, diag - 1), mask)] = d_cmp[mask];
//     mpi[c(rep(FALSE, diag - 1), mask)] = offset[mask];
//
//   }
//
//   mp[mp > 1] = 1.0;
//   mp = sqrt(2 * w * (1 - mp));
//
//   return(mp);
// }
