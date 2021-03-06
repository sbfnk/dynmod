// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// tsir
Rcpp::List tsir(SEXP rcpp_model_params, SEXP rcpp_init, SEXP rcpp_t0, SEXP rcpp_tmax, SEXP rcpp_seed, SEXP rcpp_verbose, SEXP rcpp_threads, SEXP rcpp_summary);
RcppExport SEXP dynmod_tsir(SEXP rcpp_model_paramsSEXP, SEXP rcpp_initSEXP, SEXP rcpp_t0SEXP, SEXP rcpp_tmaxSEXP, SEXP rcpp_seedSEXP, SEXP rcpp_verboseSEXP, SEXP rcpp_threadsSEXP, SEXP rcpp_summarySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type rcpp_model_params(rcpp_model_paramsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rcpp_init(rcpp_initSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rcpp_t0(rcpp_t0SEXP);
    Rcpp::traits::input_parameter< SEXP >::type rcpp_tmax(rcpp_tmaxSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rcpp_seed(rcpp_seedSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rcpp_verbose(rcpp_verboseSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rcpp_threads(rcpp_threadsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rcpp_summary(rcpp_summarySEXP);
    __result = Rcpp::wrap(tsir(rcpp_model_params, rcpp_init, rcpp_t0, rcpp_tmax, rcpp_seed, rcpp_verbose, rcpp_threads, rcpp_summary));
    return __result;
END_RCPP
}
