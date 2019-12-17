#include "mpx.hpp"
#include "math.hpp"

// [[Rcpp::export]]
List mpx_rcpp(NumericVector a, uint16_t w, uint16_t minlag, bool idxs = false, bool euclidean = true) {

  try {
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    // matrix profile using cross correlation,
    uint32_t n = a.length();

    List msd = muinvn_rcpp(a, w);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double* mu = &mmu[0];
    double* sig = &ssig[0];

    uint32_t diagmax = n - w + 1;
    IntegerVector seq_diag = Range(minlag, diagmax - 1);

    NumericVector mmp(diagmax, -1.0);
    IntegerVector mmpi(diagmax, R_NaN);

    double* mp = &mmp[0];
    int* mpi = &mmpi[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf = 0.5 * (a[Range(w, n - 1)] - a[Range(0, n - w - 1)]);
    ddf.push_front(0);
    NumericVector ddg = (a[Range(w, n - 1)] - mmu[Range(1, diagmax - 1)]) + (a[Range(0, n - w - 1)] - mmu[Range(0, n - w - 1)]);
    ddg.push_front(0);

    double* df = &ddf[0];
    double* dg = &ddg[0];

    NumericVector ww = (a[Range(0, w - 1)] - mmu[0]);

    try
    {
      for (IntegerVector::iterator diag = seq_diag.begin(); diag != seq_diag.end(); ++diag) {

        if(*diag % 1000) {
          Rcpp::checkUserInterrupt();
        }

        c = inner_product((a[Range(*diag, (*diag + w - 1))] - mu[*diag]), ww);
        off_max = (n - w - *diag + 1);
        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + *diag;
          c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
          c_cmp = c * sig[offset] * sig[off_diag];
          if (c_cmp > mp[offset]) {
            mp[offset] = c_cmp;
            if(idxs) {
              mpi[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp[off_diag]) {
            mp[off_diag] = c_cmp;
            if(idxs) {
              mpi[off_diag] = offset + 1;
            }
          }
        }
      }
    } catch (internal::InterruptedException &ex)
    {
      Rcout << "Process terminated.\n";
    }

    // to do ed
    mmp[mmp > 1.0] = 1.0;

    if(euclidean) { // correlation to ed
      mmp = sqrt(2 * w * (1 - mmp));
    }

    if(idxs) {
      return (List::create(
          Rcpp::Named("mp") = mmp,
          Rcpp::Named("pi") = mmpi
      ));
    } else {
      return (List::create(
          Rcpp::Named("mp") = mmp
      ));
    }
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// [[Rcpp::export]]
List mpxab_rcpp(NumericVector a, NumericVector b, uint16_t w, bool idxs = false, bool euclidean = true) {

  try
  {
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    // matrix profile using cross correlation,
    uint32_t a_len = a.length();
    uint32_t b_len = b.length();

    List msd_a = muinvn_rcpp(a, w);
    List msd_b = muinvn_rcpp(b, w);

    NumericVector mmu_a = msd_a["avg"];
    NumericVector ssig_a = msd_a["sig"];
    NumericVector mmu_b = msd_b["avg"];
    NumericVector ssig_b = msd_b["sig"];
    double* mu_a = &mmu_a[0];
    double* sig_a = &ssig_a[0];
    double* mu_b = &mmu_b[0];
    double* sig_b = &ssig_b[0];

    uint32_t diagmax_a = a_len - w + 1;
    IntegerVector seq_diag = Range(0, diagmax_a - 1);

    NumericVector mmp_a(diagmax_a, -1.0);
    IntegerVector mmpi_a(diagmax_a, R_NaN);

    double* mp_a = &mmp_a[0];
    int* mpi_a = &mmpi_a[0];

    uint32_t diagmax_b = b_len - w + 1;

    NumericVector mmp_b(diagmax_b, -1.0);
    IntegerVector mmpi_b(diagmax_b, R_NaN);

    double* mp_b = &mmp_b[0];
    int* mpi_b = &mmpi_b[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf_a = 0.5 * (a[Range(w, a_len - 1)] - a[Range(0, a_len - w - 1)]);
    ddf_a.push_front(0);
    NumericVector ddg_a = (a[Range(w, a_len - 1)] - mmu_a[Range(1, diagmax_a - 1)]) + (a[Range(0, a_len - w - 1)] - mmu_a[Range(0, a_len - w - 1)]);
    ddg_a.push_front(0);
    NumericVector ddf_b = 0.5 * (b[Range(w, b_len - 1)] - b[Range(0, b_len - w - 1)]);
    ddf_b.push_front(0);
    NumericVector ddg_b = (b[Range(w, a_len - 1)] - mmu_b[Range(1, diagmax_b - 1)]) + (b[Range(0, b_len - w - 1)] - mmu_b[Range(0, b_len - w - 1)]);
    ddg_b.push_front(0);

    double* df_a = &ddf_a[0];
    double* dg_a = &ddg_a[0];
    double* df_b = &ddf_b[0];
    double* dg_b = &ddg_b[0];

    NumericVector ww = (b[Range(0, w - 1)] - mmu_b[0]);

    try
    {
      for (IntegerVector::iterator diag = seq_diag.begin(); diag != seq_diag.end(); ++diag) {
        if(*diag % 1000) {
          Rcpp::checkUserInterrupt();
        }

        off_max = MIN(a_len - w - *diag + 1, b_len - w + 1);
        c = inner_product((a[Range(*diag, (*diag + w - 1))] - mu_a[*diag]), ww);
        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + *diag;
          c = c + df_a[off_diag] * dg_b[offset] + dg_a[off_diag] * df_b[offset];
          c_cmp = c * sig_b[offset] * sig_a[off_diag];
          if (c_cmp > mp_b[offset]) {
            mp_b[offset] = c_cmp;
            if(idxs) {
              mpi_b[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp_a[off_diag]) {
            mp_a[off_diag] = c_cmp;
            if(idxs) {
              mpi_a[off_diag] = offset + 1;
            }
          }
        }
      }

      ww = (a[Range(0, w - 1)] - mmu_a[0]);
      seq_diag = Range(0, diagmax_b - 1);

      for (IntegerVector::iterator diag = seq_diag.begin(); diag != seq_diag.end(); ++diag) {
        if(*diag % 1000) {
          Rcpp::checkUserInterrupt();
        }

        off_max = MIN(b_len - w - *diag + 1, a_len - w + 1);
        c = inner_product((b[Range(*diag, (*diag + w - 1))] - mu_b[*diag]), ww);
        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + *diag;
          c = c + df_b[off_diag] * dg_a[offset] + dg_b[off_diag] * df_a[offset];
          c_cmp = c * sig_a[offset] * sig_b[off_diag];
          if (c_cmp > mp_a[offset]) {
            mp_a[offset] = c_cmp;
            if(idxs) {
              mpi_a[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp_b[off_diag]) {
            mp_b[off_diag] = c_cmp;
            if(idxs) {
              mpi_b[off_diag] = offset + 1;
            }
          }
        }
      }
    } catch (internal::InterruptedException &ex)
    {
      Rcout  << "Process terminated.\n";
    }

    // to do ed
    mmp_a[mmp_a > 1.0] = 1.0;
    mmp_b[mmp_b > 1.0] = 1.0;

    if(euclidean) { // correlation to ed
      mmp_a = sqrt(2 * w * (1 - mmp_a));
      mmp_b = sqrt(2 * w * (1 - mmp_b));
    }

    if(idxs) {
      return (List::create(
          Rcpp::Named("mp") = mmp_a,
          Rcpp::Named("pi") = mmpi_a,
          Rcpp::Named("mpb") = mmp_b,
          Rcpp::Named("pib") = mmpi_b
      ));
    } else {
      return (List::create(
          Rcpp::Named("mp") = mmp_a,
          Rcpp::Named("mpb") = mmp_b
      ));
    }
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}
