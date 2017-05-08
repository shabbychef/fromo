// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sd3
NumericVector sd3(SEXP v, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_sd3(SEXP vSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(sd3(v, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// skew4
NumericVector skew4(SEXP v, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_skew4(SEXP vSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(skew4(v, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// kurt5
NumericVector kurt5(SEXP v, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_kurt5(SEXP vSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(kurt5(v, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// cent_moments
NumericVector cent_moments(SEXP v, int max_order, int used_df, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_cent_moments(SEXP vSEXP, SEXP max_orderSEXP, SEXP used_dfSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(cent_moments(v, max_order, used_df, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// std_moments
NumericVector std_moments(SEXP v, int max_order, int used_df, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_std_moments(SEXP vSEXP, SEXP max_orderSEXP, SEXP used_dfSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(std_moments(v, max_order, used_df, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// cent_cumulants
NumericVector cent_cumulants(SEXP v, int max_order, int used_df, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_cent_cumulants(SEXP vSEXP, SEXP max_orderSEXP, SEXP used_dfSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(cent_cumulants(v, max_order, used_df, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// std_cumulants
NumericVector std_cumulants(SEXP v, int max_order, int used_df, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_std_cumulants(SEXP vSEXP, SEXP max_orderSEXP, SEXP used_dfSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(std_cumulants(v, max_order, used_df, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// cent_sums
NumericVector cent_sums(SEXP v, int max_order, bool na_rm, SEXP wts, bool check_wts);
RcppExport SEXP fromo_cent_sums(SEXP vSEXP, SEXP max_orderSEXP, SEXP na_rmSEXP, SEXP wtsSEXP, SEXP check_wtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type wts(wtsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_wts(check_wtsSEXP);
    rcpp_result_gen = Rcpp::wrap(cent_sums(v, max_order, na_rm, wts, check_wts));
    return rcpp_result_gen;
END_RCPP
}
// join_cent_sums
NumericVector join_cent_sums(NumericVector ret1, NumericVector ret2);
RcppExport SEXP fromo_join_cent_sums(SEXP ret1SEXP, SEXP ret2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ret1(ret1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ret2(ret2SEXP);
    rcpp_result_gen = Rcpp::wrap(join_cent_sums(ret1, ret2));
    return rcpp_result_gen;
END_RCPP
}
// unjoin_cent_sums
NumericVector unjoin_cent_sums(NumericVector ret3, NumericVector ret2);
RcppExport SEXP fromo_unjoin_cent_sums(SEXP ret3SEXP, SEXP ret2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ret3(ret3SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ret2(ret2SEXP);
    rcpp_result_gen = Rcpp::wrap(unjoin_cent_sums(ret3, ret2));
    return rcpp_result_gen;
END_RCPP
}
// cent_cosums
NumericMatrix cent_cosums(SEXP v, int max_order, bool na_omit);
RcppExport SEXP fromo_cent_cosums(SEXP vSEXP, SEXP max_orderSEXP, SEXP na_omitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type na_omit(na_omitSEXP);
    rcpp_result_gen = Rcpp::wrap(cent_cosums(v, max_order, na_omit));
    return rcpp_result_gen;
END_RCPP
}
// cent_comoments
NumericMatrix cent_comoments(SEXP v, int max_order, int used_df, bool na_omit);
RcppExport SEXP fromo_cent_comoments(SEXP vSEXP, SEXP max_orderSEXP, SEXP used_dfSEXP, SEXP na_omitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< bool >::type na_omit(na_omitSEXP);
    rcpp_result_gen = Rcpp::wrap(cent_comoments(v, max_order, used_df, na_omit));
    return rcpp_result_gen;
END_RCPP
}
// join_cent_cosums
NumericMatrix join_cent_cosums(NumericMatrix ret1, NumericMatrix ret2);
RcppExport SEXP fromo_join_cent_cosums(SEXP ret1SEXP, SEXP ret2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ret1(ret1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ret2(ret2SEXP);
    rcpp_result_gen = Rcpp::wrap(join_cent_cosums(ret1, ret2));
    return rcpp_result_gen;
END_RCPP
}
// unjoin_cent_cosums
NumericMatrix unjoin_cent_cosums(NumericMatrix ret3, NumericMatrix ret2);
RcppExport SEXP fromo_unjoin_cent_cosums(SEXP ret3SEXP, SEXP ret2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ret3(ret3SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ret2(ret2SEXP);
    rcpp_result_gen = Rcpp::wrap(unjoin_cent_cosums(ret3, ret2));
    return rcpp_result_gen;
END_RCPP
}
// running_sd3
NumericMatrix running_sd3(SEXP v, SEXP window, bool na_rm, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_sd3(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_sd3(v, window, na_rm, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_skew4
NumericMatrix running_skew4(SEXP v, SEXP window, bool na_rm, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_skew4(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_skew4(v, window, na_rm, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_kurt5
NumericMatrix running_kurt5(SEXP v, SEXP window, bool na_rm, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_kurt5(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_kurt5(v, window, na_rm, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_cent_moments
NumericMatrix running_cent_moments(SEXP v, SEXP window, int max_order, bool na_rm, bool max_order_only, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_cent_moments(SEXP vSEXP, SEXP windowSEXP, SEXP max_orderSEXP, SEXP na_rmSEXP, SEXP max_order_onlySEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< bool >::type max_order_only(max_order_onlySEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_cent_moments(v, window, max_order, na_rm, max_order_only, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_std_moments
NumericMatrix running_std_moments(SEXP v, SEXP window, int max_order, bool na_rm, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_std_moments(SEXP vSEXP, SEXP windowSEXP, SEXP max_orderSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_std_moments(v, window, max_order, na_rm, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_cumulants
NumericMatrix running_cumulants(SEXP v, SEXP window, int max_order, bool na_rm, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_cumulants(SEXP vSEXP, SEXP windowSEXP, SEXP max_orderSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_cumulants(v, window, max_order, na_rm, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_apx_quantiles
NumericMatrix running_apx_quantiles(SEXP v, NumericVector p, SEXP window, int max_order, bool na_rm, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_apx_quantiles(SEXP vSEXP, SEXP pSEXP, SEXP windowSEXP, SEXP max_orderSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_apx_quantiles(v, p, window, max_order, na_rm, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_apx_median
NumericMatrix running_apx_median(SEXP v, SEXP window, int max_order, bool na_rm, int min_df, int used_df, int restart_period);
RcppExport SEXP fromo_running_apx_median(SEXP vSEXP, SEXP windowSEXP, SEXP max_orderSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP used_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< int >::type max_order(max_orderSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type used_df(used_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_apx_median(v, window, max_order, na_rm, min_df, used_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_centered
NumericMatrix running_centered(SEXP v, SEXP window, bool na_rm, int min_df, int lookahead, int restart_period);
RcppExport SEXP fromo_running_centered(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP lookaheadSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type lookahead(lookaheadSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_centered(v, window, na_rm, min_df, lookahead, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_scaled
NumericMatrix running_scaled(SEXP v, SEXP window, bool na_rm, int min_df, int lookahead, int restart_period);
RcppExport SEXP fromo_running_scaled(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP lookaheadSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type lookahead(lookaheadSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_scaled(v, window, na_rm, min_df, lookahead, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_zscored
NumericMatrix running_zscored(SEXP v, SEXP window, bool na_rm, int min_df, int lookahead, int restart_period);
RcppExport SEXP fromo_running_zscored(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP lookaheadSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type lookahead(lookaheadSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_zscored(v, window, na_rm, min_df, lookahead, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_sharpe
NumericMatrix running_sharpe(SEXP v, SEXP window, bool na_rm, bool compute_se, int min_df, int restart_period);
RcppExport SEXP fromo_running_sharpe(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP compute_seSEXP, SEXP min_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< bool >::type compute_se(compute_seSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_sharpe(v, window, na_rm, compute_se, min_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// running_tstat
NumericMatrix running_tstat(SEXP v, SEXP window, bool na_rm, int min_df, int restart_period);
RcppExport SEXP fromo_running_tstat(SEXP vSEXP, SEXP windowSEXP, SEXP na_rmSEXP, SEXP min_dfSEXP, SEXP restart_periodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type v(vSEXP);
    Rcpp::traits::input_parameter< SEXP >::type window(windowSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< int >::type min_df(min_dfSEXP);
    Rcpp::traits::input_parameter< int >::type restart_period(restart_periodSEXP);
    rcpp_result_gen = Rcpp::wrap(running_tstat(v, window, na_rm, min_df, restart_period));
    return rcpp_result_gen;
END_RCPP
}
// cent2raw
NumericVector cent2raw(NumericVector input);
RcppExport SEXP fromo_cent2raw(SEXP inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input(inputSEXP);
    rcpp_result_gen = Rcpp::wrap(cent2raw(input));
    return rcpp_result_gen;
END_RCPP
}
