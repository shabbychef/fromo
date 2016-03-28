// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sd3
NumericVector sd3(SEXP v, bool na_rm);
RcppExport SEXP fromo_sd3(SEXP vSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(sd3(v, na_rm));
    return __result;
END_RCPP
}
// skew4
NumericVector skew4(SEXP v, bool na_rm);
RcppExport SEXP fromo_skew4(SEXP vSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(skew4(v, na_rm));
    return __result;
END_RCPP
}
// kurt5
NumericVector kurt5(SEXP v, bool na_rm);
RcppExport SEXP fromo_kurt5(SEXP vSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(kurt5(v, na_rm));
    return __result;
END_RCPP
}
// cent_moments
NumericVector cent_moments(SEXP v, int max_order, int used_df, bool na_rm);
RcppExport SEXP fromo_cent_moments(SEXP vSEXP, SEXP max_orderSEXP, SEXP used_dfSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(cent_moments(v, max_order, used_df, na_rm));
    return __result;
END_RCPP
}
// run_sd3
NumericMatrix run_sd3(SEXP v, int winsize, int recoper, bool na_rm);
RcppExport SEXP fromo_run_sd3(SEXP vSEXP, SEXP winsizeSEXP, SEXP recoperSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type winsize(winsizeSEXP);
    Rcpp::traits::input_parameter< int >::type recoper(recoperSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(run_sd3(v, winsize, recoper, na_rm));
    return __result;
END_RCPP
}
// run_skew4
NumericMatrix run_skew4(SEXP v, int winsize, int recoper, bool na_rm);
RcppExport SEXP fromo_run_skew4(SEXP vSEXP, SEXP winsizeSEXP, SEXP recoperSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type winsize(winsizeSEXP);
    Rcpp::traits::input_parameter< int >::type recoper(recoperSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(run_skew4(v, winsize, recoper, na_rm));
    return __result;
END_RCPP
}
// run_kurt5
NumericMatrix run_kurt5(SEXP v, int winsize, int recoper, bool na_rm);
RcppExport SEXP fromo_run_kurt5(SEXP vSEXP, SEXP winsizeSEXP, SEXP recoperSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type winsize(winsizeSEXP);
    Rcpp::traits::input_parameter< int >::type recoper(recoperSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(run_kurt5(v, winsize, recoper, na_rm));
    return __result;
END_RCPP
}
// run_centered
NumericMatrix run_centered(SEXP v, int winsize, int recoper, bool na_rm);
RcppExport SEXP fromo_run_centered(SEXP vSEXP, SEXP winsizeSEXP, SEXP recoperSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type winsize(winsizeSEXP);
    Rcpp::traits::input_parameter< int >::type recoper(recoperSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(run_centered(v, winsize, recoper, na_rm));
    return __result;
END_RCPP
}
// run_scaled
NumericMatrix run_scaled(SEXP v, int winsize, int recoper, bool na_rm);
RcppExport SEXP fromo_run_scaled(SEXP vSEXP, SEXP winsizeSEXP, SEXP recoperSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type winsize(winsizeSEXP);
    Rcpp::traits::input_parameter< int >::type recoper(recoperSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(run_scaled(v, winsize, recoper, na_rm));
    return __result;
END_RCPP
}
// run_zscored
NumericMatrix run_zscored(SEXP v, int winsize, int recoper, bool na_rm);
RcppExport SEXP fromo_run_zscored(SEXP vSEXP, SEXP winsizeSEXP, SEXP recoperSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type winsize(winsizeSEXP);
    Rcpp::traits::input_parameter< int >::type recoper(recoperSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(run_zscored(v, winsize, recoper, na_rm));
    return __result;
END_RCPP
}
// run_tscored
NumericMatrix run_tscored(SEXP v, int winsize, int recoper, bool na_rm);
RcppExport SEXP fromo_run_tscored(SEXP vSEXP, SEXP winsizeSEXP, SEXP recoperSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type winsize(winsizeSEXP);
    Rcpp::traits::input_parameter< int >::type recoper(recoperSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    __result = Rcpp::wrap(run_tscored(v, winsize, recoper, na_rm));
    return __result;
END_RCPP
}
