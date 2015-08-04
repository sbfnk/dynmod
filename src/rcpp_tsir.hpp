// -*- compile-command: "cd .. ; R CMD INSTALL .; cd -"; -*-
/*! \file rcpp_tsir.cpp
  \brief Implementation of the wrapper for running the tsir model form R
*/
#ifndef _tsir_RCPP_TSIR_H
#define _tsir_RCPP_TSIR_H

#include <Rcpp.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

//! Wrapper for R function for running the TSIR model
/*!
  \param[in] _model_params List of model parameters (alpha, beta,
  tau1, tau2, theta, rho) 
  \param[in] _sim_params List of simulation parameters (verbose,
  startTime, endTime, allowRejection, use_seed, threads)
  \param[in] _model_options List of model option (nRecovered)
  \param[in] _summary Summary statistics to calculate
  \param[in] _init Initial values
  \param[in] _data List of model data (births, population sizes, city
  coordinates, vaccination rates, seasonality)
  \return A list of model outputs (summary statistics, whether the
  infection went extinct)
 */
RcppExport SEXP rcpp_tsir(SEXP _model_params, SEXP _model_options,
                          SEXP _sim_params, SEXP _summary, 
                          SEXP _init, SEXP _data, SEXP _seed);

#endif
