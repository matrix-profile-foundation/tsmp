#include "math.h" // math first to fix OSX error
#include "mpx.h"

// MPX
//
// @param a Time Series
// @return a List
// [[Rcpp::export]]

List mpx_rcpp(NumericVector a, uint16_t w, uint16_t minlag, bool idxs = false, bool euclidean = true) {

  try {
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    bool partial = false;
    // matrix profile using cross correlation,
    uint32_t n = a.length();

    List msd = muinvn_rcpp(a, w);

    NumericVector mmu = msd["avg"];
    NumericVector ssig = msd["sig"];
    double *mu = &mmu[0];
    double *sig = &ssig[0];

    uint32_t profile_len = n - w + 1;
    IntegerVector seq_diag = Range(minlag, profile_len - 1);

    NumericVector mmp(profile_len, -1.0);
    IntegerVector mmpi(profile_len, R_NaN);

    double *mp = &mmp[0];
    int *mpi = &mmpi[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf = 0.5 * (a[Range(w, n - 1)] - a[Range(0, n - w - 1)]);
    ddf.push_front(0);
    NumericVector ddg = (a[Range(w, n - 1)] - mmu[Range(1, profile_len - 1)]) + (a[Range(0, n - w - 1)] - mmu[Range(0, n - w - 1)]);
    ddg.push_front(0);

    double *df = &ddf[0];
    double *dg = &ddg[0];

    NumericVector ww = (a[Range(0, w - 1)] - mmu[0]);

    try {
      for (IntegerVector::iterator diag = seq_diag.begin(); diag != seq_diag.end(); ++diag) {

        if (*diag % 1000) {
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
            if (idxs) {
              mpi[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp[off_diag]) {
            mp[off_diag] = c_cmp;
            if (idxs) {
              mpi[off_diag] = offset + 1;
            }
          }
        }
      }
    } catch (internal::InterruptedException &ex) {
      partial = true;
      Rcout << "Process terminated.\n";
    }

    // to do ed
    mmp[mmp > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mmp = sqrt(2 * w * (1 - mmp));
    }

    if (idxs) {
      return (List::create(
          Rcpp::Named("mp") = mmp,
          Rcpp::Named("pi") = mmpi,
          Rcpp::Named("partial") = partial
      ));
    } else {
      return (List::create(
          Rcpp::Named("mp") = mmp,
          Rcpp::Named("partial") = partial
      ));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

// [[Rcpp::export]]
List mpxab_rcpp(NumericVector a, NumericVector b, uint16_t w, bool idxs = false, bool euclidean = true) {

  try {
    bool partial = false;
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
    double *mu_a = &mmu_a[0];
    double *sig_a = &ssig_a[0];
    double *mu_b = &mmu_b[0];
    double *sig_b = &ssig_b[0];

    uint32_t profile_len_a = a_len - w + 1;
    IntegerVector seq_diag = Range(0, profile_len_a - 1);

    NumericVector mmp_a(profile_len_a, -1.0);
    IntegerVector mmpi_a(profile_len_a, R_NaN);

    double *mp_a = &mmp_a[0];
    int *mpi_a = &mmpi_a[0];

    uint32_t profile_len_b = b_len - w + 1;

    NumericVector mmp_b(profile_len_b, -1.0);
    IntegerVector mmpi_b(profile_len_b, R_NaN);

    double *mp_b = &mmp_b[0];
    int *mpi_b = &mmpi_b[0];

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector ddf_a = 0.5 * (a[Range(w, a_len - 1)] - a[Range(0, a_len - w - 1)]);
    ddf_a.push_front(0);
    NumericVector ddg_a = (a[Range(w, a_len - 1)] - mmu_a[Range(1, profile_len_a - 1)]) + (a[Range(0, a_len - w - 1)] - mmu_a[Range(0, a_len - w - 1)]);
    ddg_a.push_front(0);
    NumericVector ddf_b = 0.5 * (b[Range(w, b_len - 1)] - b[Range(0, b_len - w - 1)]);
    ddf_b.push_front(0);
    NumericVector ddg_b = (b[Range(w, a_len - 1)] - mmu_b[Range(1, profile_len_b - 1)]) + (b[Range(0, b_len - w - 1)] - mmu_b[Range(0, b_len - w - 1)]);
    ddg_b.push_front(0);

    double *df_a = &ddf_a[0];
    double *dg_a = &ddg_a[0];
    double *df_b = &ddf_b[0];
    double *dg_b = &ddg_b[0];

    NumericVector ww = (b[Range(0, w - 1)] - mmu_b[0]);

    try {
      for (IntegerVector::iterator diag = seq_diag.begin(); diag != seq_diag.end(); ++diag) {
        if (*diag % 1000) {
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
            if (idxs) {
              mpi_b[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp_a[off_diag]) {
            mp_a[off_diag] = c_cmp;
            if (idxs) {
              mpi_a[off_diag] = offset + 1;
            }
          }
        }
      }

      ww = (a[Range(0, w - 1)] - mmu_a[0]);
      seq_diag = Range(0, profile_len_b - 1);

      for (IntegerVector::iterator diag = seq_diag.begin(); diag != seq_diag.end(); ++diag) {
        if (*diag % 1000) {
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
            if (idxs) {
              mpi_a[offset] = off_diag + 1;
            }
          }
          if (c_cmp > mp_b[off_diag]) {
            mp_b[off_diag] = c_cmp;
            if (idxs) {
              mpi_b[off_diag] = offset + 1;
            }
          }
        }
      }
    } catch (internal::InterruptedException &ex) {
      partial = true;
      Rcout  << "Process terminated.\n";
    }

    // to do ed
    mmp_a[mmp_a > 1.0] = 1.0;
    mmp_b[mmp_b > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mmp_a = sqrt(2 * w * (1 - mmp_a));
      mmp_b = sqrt(2 * w * (1 - mmp_b));
    }

    if (idxs) {
      return (List::create(
          Rcpp::Named("mp") = mmp_a,
          Rcpp::Named("pi") = mmpi_a,
          Rcpp::Named("mpb") = mmp_b,
          Rcpp::Named("pib") = mmpi_b,
          Rcpp::Named("partial") = partial
      ));
    } else {
      return (List::create(
          Rcpp::Named("mp") = mmp_a,
          Rcpp::Named("mpb") = mmp_b,
          Rcpp::Named("partial") = partial
      ));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

//##### Parallel version #####

struct MatrixProfileP : public Worker {
  // input
  const RVector<double> a;
  const uint16_t w;
  const RVector<double> df;
  const RVector<double> dg;
  const RVector<double> mu;
  const RVector<double> sig;
  const RVector<double> ww;

  // output
  RVector<double> mp;
  RVector<int> mpi;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  MatrixProfileP(const NumericVector a, const uint16_t w, const NumericVector df, const NumericVector dg,
                 const NumericVector mmu, const NumericVector sig, const NumericVector ww,
                 NumericVector mp, IntegerVector mpi) :
    a(a), w(w), df(df), dg(dg), mu(mmu), sig(sig), ww(ww), mp(mp), mpi(mpi) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) { // minlag:profile_len
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    uint32_t n = a.length();
    std::vector<double> aa(w);

    for (uint32_t diag = begin; diag < end; diag++) {

      for (int i = 0; i < w; i++) {
        aa[i] = a[diag + i] - mu[diag];
      }

      c = std::inner_product(aa.begin(), aa.end(), ww.begin(), 0.0);

      off_max = (n - w - diag + 1);

      for (offset = 0; offset < off_max; offset++) {
        off_diag = offset + diag;
        c = c + df[offset] * dg[off_diag] + df[off_diag] * dg[offset];
        c_cmp = c * sig[offset] * sig[off_diag];
        if (c_cmp > mp[offset]) {
          mp[offset] = c_cmp;
          mpi[offset] = off_diag + 1;
        }
        if (c_cmp > mp[off_diag]) {
          mp[off_diag] = c_cmp;
          mpi[off_diag] = offset + 1;
        }
      }
    }
  }
};

// [[Rcpp::export]]
List mpx_rcpp_parallel(NumericVector a, uint16_t w, uint16_t minlag, bool idxs = false, bool euclidean = true) {

  try {
    // matrix profile using cross correlation,
    bool partial = false;
    uint32_t n = a.length();

    List msd = muinvn_rcpp(a, w);

    NumericVector mu = msd["avg"];
    NumericVector sig = msd["sig"];

    uint32_t profile_len = n - w + 1;
    NumericVector mp(profile_len, -1.0);
    IntegerVector mpi(profile_len, R_NaN);

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.

    NumericVector df = 0.5 * (a[Range(w, n - 1)] - a[Range(0, n - w - 1)]);
    df.push_front(0);
    NumericVector dg = (a[Range(w, n - 1)] - mu[Range(1, profile_len - 1)]) + (a[Range(0, n - w - 1)] - mu[Range(0, n - w - 1)]);
    dg.push_front(0);

    NumericVector ww = (a[Range(0, w - 1)] - mu[0]);

    MatrixProfileP matrix_profile(a, w, df, dg, mu, sig, ww, mp, mpi);

    try {
      parallelFor(minlag, profile_len, matrix_profile, 100);
    } catch (internal::InterruptedException &ex) {
      partial = true;
      Rcout << "Process terminated.\n";
    }

    // to do ed
    mp[mp > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mp = sqrt(2 * w * (1 - mp));
    }

    if (idxs) {
      return (List::create(
          Rcpp::Named("mp") = mp,
          Rcpp::Named("pi") = mpi,
          Rcpp::Named("partial") = partial
      ));
    } else {
      return (List::create(
          Rcpp::Named("mp") = mp,
          Rcpp::Named("partial") = partial
      ));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}

struct MatrixProfilePAB : public Worker {
  // input
  const RVector<double> a;
  const RVector<double> b;
  const uint16_t w;
  const RVector<double> df_a;
  const RVector<double> df_b;
  const RVector<double> dg_a;
  const RVector<double> dg_b;
  const RVector<double> mu_a;
  const RVector<double> mu_b;
  const RVector<double> sig_a;
  const RVector<double> sig_b;
  const RVector<double> ww_a;
  const RVector<double> ww_b;

  // output
  RVector<double> mp_a;
  RVector<double> mp_b;
  RVector<int> mpi_a;
  RVector<int> mpi_b;

  // AB == 0, BA == 1
  uint8_t ab_ba;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  MatrixProfilePAB(const NumericVector a, const NumericVector b, const uint16_t w,
                   const NumericVector df_a, const NumericVector df_b,
                   const NumericVector dg_a, const NumericVector dg_b,
                   const NumericVector mu_a, const NumericVector mu_b,
                   const NumericVector sig_a, const NumericVector sig_b,
                   const NumericVector ww_a, const NumericVector ww_b,
                   NumericVector mp_a, NumericVector mp_b,
                   IntegerVector mpi_a, IntegerVector mpi_b) :
                    a(a), b(b), w(w),
                    df_a(df_a), df_b(df_b),
                    dg_a(dg_a), dg_b(dg_b),
                    mu_a(mu_a), mu_b(mu_b),
                    sig_a(sig_a), sig_b(sig_b),
                    ww_a(ww_a), ww_b(ww_b),
                    mp_a(mp_a), mp_b(mp_b),
                    mpi_a(mpi_a), mpi_b(mpi_b),
                    ab_ba(0) {}



  void set_ab() {
    this->ab_ba = 0;
  }

  void set_ba() {
    this->ab_ba = 1;
  }

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) { // minlag:profile_len
    double c, c_cmp;
    uint32_t off_max, off_diag, offset;
    uint32_t a_len = a.length();
    uint32_t b_len = b.length();
    std::vector<double> inn(w);

    if(ab_ba == 0) {

      for (uint32_t diag = begin; diag < end; diag++) {

        for (int i = 0; i < w; i++) {
          inn[i] = a[diag + i] - mu_a[diag];
        }

        off_max = MIN(a_len - w - diag + 1, b_len - w + 1);

        c = std::inner_product(inn.begin(), inn.end(), ww_b.begin(), 0.0);

        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          c = c + df_a[off_diag] * dg_b[offset] + dg_a[off_diag] * df_b[offset];
          c_cmp = c * sig_b[offset] * sig_a[off_diag];
          if (c_cmp > mp_b[offset]) {
            mp_b[offset] = c_cmp;
            mpi_b[offset] = off_diag + 1;
          }
          if (c_cmp > mp_a[off_diag]) {
            mp_a[off_diag] = c_cmp;
            mpi_a[off_diag] = offset + 1;
          }
        }
      }
    } else {

      for (uint32_t diag = begin; diag < end; diag++) {

        for (int i = 0; i < w; i++) {
          inn[i] = b[diag + i] - mu_b[diag];
        }

        off_max = MIN(b_len - w - diag + 1, a_len - w + 1);

        c = std::inner_product(inn.begin(), inn.end(), ww_a.begin(), 0.0);

        for (offset = 0; offset < off_max; offset++) {
          off_diag = offset + diag;
          c = c + df_b[off_diag] * dg_a[offset] + dg_b[off_diag] * df_a[offset];
          c_cmp = c * sig_a[offset] * sig_b[off_diag];
          if (c_cmp > mp_a[offset]) {
            mp_a[offset] = c_cmp;
            mpi_a[offset] = off_diag + 1;
          }
          if (c_cmp > mp_b[off_diag]) {
            mp_b[off_diag] = c_cmp;
            mpi_b[off_diag] = offset + 1;
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
List mpxab_rcpp_parallel(NumericVector a, NumericVector b, uint16_t w, bool idxs = false, bool euclidean = true) {

  try {
    // matrix profile using cross correlation,
    bool partial = false;
    uint32_t a_len = a.length();
    uint32_t b_len = b.length();

    List msd_a = muinvn_rcpp(a, w);
    List msd_b = muinvn_rcpp(b, w);

    NumericVector mu_a = msd_a["avg"];
    NumericVector sig_a = msd_a["sig"];
    NumericVector mu_b = msd_b["avg"];
    NumericVector sig_b = msd_b["sig"];

    uint32_t profile_len_a = a_len - w + 1;
    uint32_t profile_len_b = b_len - w + 1;

    NumericVector mp_a(profile_len_a, -1.0);
    NumericVector mp_b(profile_len_b, -1.0);

    IntegerVector mpi_a(profile_len_a, R_NaN);
    IntegerVector mpi_b(profile_len_b, R_NaN);

    // differentials have 0 as their first entry. This simplifies index
    // calculations slightly and allows us to avoid special "first line"
    // handling.
    NumericVector df_a = 0.5 * (a[Range(w, a_len - 1)] - a[Range(0, a_len - w - 1)]);
    df_a.push_front(0);
    NumericVector dg_a = (a[Range(w, a_len - 1)] - mu_a[Range(1, profile_len_a - 1)]) + (a[Range(0, a_len - w - 1)] - mu_a[Range(0, a_len - w - 1)]);
    dg_a.push_front(0);
    NumericVector df_b = 0.5 * (b[Range(w, b_len - 1)] - b[Range(0, b_len - w - 1)]);
    df_b.push_front(0);
    NumericVector dg_b = (b[Range(w, a_len - 1)] - mu_b[Range(1, profile_len_b - 1)]) + (b[Range(0, b_len - w - 1)] - mu_b[Range(0, b_len - w - 1)]);
    dg_b.push_front(0);

    NumericVector ww_a = (a[Range(0, w - 1)] - mu_a[0]);
    NumericVector ww_b = (b[Range(0, w - 1)] - mu_b[0]);

    MatrixProfilePAB matrix_profile(a, b, w,
                                    df_a, df_b,
                                    dg_a, dg_b,
                                    mu_a, mu_b,
                                    sig_a, sig_b,
                                    ww_a, ww_b,
                                    mp_a, mp_b,
                                    mpi_a, mpi_b);

    try {
      parallelFor(0, profile_len_a, matrix_profile);
    } catch (internal::InterruptedException &ex) {
      partial = true;
      Rcout << "Process AB terminated.\n";
    }

    // TODO: make a thread-safe interrupt check for stop the process and try to return the partial computation
    Rcpp::checkUserInterrupt();

    matrix_profile.set_ba();

    try {
      parallelFor(0, profile_len_b, matrix_profile);
    } catch (internal::InterruptedException &ex) {
      partial = true;
      Rcout << "Process BA terminated.\n";
    }

    // to do ed
    mp_a[mp_a > 1.0] = 1.0;
    mp_b[mp_b > 1.0] = 1.0;

    if (euclidean) { // correlation to ed
      mp_a = sqrt(2 * w * (1 - mp_a));
      mp_b = sqrt(2 * w * (1 - mp_b));
    }

    if (idxs) {
      return (List::create(
          Rcpp::Named("mp") = mp_a,
          Rcpp::Named("pi") = mpi_a,
          Rcpp::Named("mpb") = mp_b,
          Rcpp::Named("pib") = mpi_b,
          Rcpp::Named("partial") = partial
      ));
    } else {
      return (List::create(
          Rcpp::Named("mp") = mp_a,
          Rcpp::Named("mpb") = mp_b,
          Rcpp::Named("partial") = partial
      ));
    }
  } catch (...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
}
