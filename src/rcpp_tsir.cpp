// -*- compile-command: "cd .. ; R CMD INSTALL .; cd -"; -*-
/*! \file rcpp_tsir.cpp
  \brief Implementation of the wrapper for running the tsir model form R
*/

#include <Rcpp.h>

#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp>
#include <iterator>
#include <vector>

#include "simulation.hpp"
#include "getseed.hpp"
#include "rcpp_tsir.hpp"


// [[Rcpp::export]]
Rcpp::List tsir(SEXP rcpp_model_params, SEXP rcpp_init,
                SEXP rcpp_t0, SEXP rcpp_tmax,
                SEXP rcpp_seed, SEXP rcpp_verbose,
                SEXP rcpp_threads, SEXP rcpp_summary)
{

  // Get R variables
  Rcpp::List model_params(rcpp_model_params); //!< model parameters

  // read simulation parameters
  unsigned int verbose = Rcpp::as<bool>(rcpp_verbose);
  unsigned int startTime = Rcpp::as<double>(rcpp_t0);
  unsigned int endTime = Rcpp::as<double>(rcpp_tmax);

  // check how many threads we want to run on
  unsigned int nThreads = Rcpp::as<unsigned int>(rcpp_threads);

  // if nThreads is 0, test the maximum of threads the hardware can
  // run in parallel
  if (nThreads == 0) {
    nThreads = boost::thread::hardware_concurrency();
  }

  // initialise gravity model
  GravityModel model;

  // read nRecovered classes (i.e., number of time steps one is
  // recovered -- if 0, this infinite)
  if (model_params.containsElementNamed("nRecovered")) {
    model.nRecovered = Rcpp::as<size_t>(model_params["nRecovered"]);
  } else {
    model.nRecovered = 0;
  }

  // assign model parameters
  std::vector<std::string> paramNames = model_params.names();
  for (std::vector<std::string>::iterator it = paramNames.begin();
       it != paramNames.end(); it++) {
      if (TYPEOF(model_params[*it]) == INTSXP ||
          TYPEOF(model_params[*it]) == REALSXP) {
          std::vector<double> paramVector =
              Rcpp::as<std::vector<double> >(model_params[*it]);
          // Here we're looking for scalar model parameters -- it's a
          // vector, we only give the first value to setParam take the
          // first value, the rest is ignored. SetParam will ignore
          // all parameters that are not found as scalar parameters in
          // the model
          bool found = model.setParam(*it, paramVector[0]);
          // If the model found it as a scalar parameter but we were
          // given a vector, give a warning
          if (found && paramVector.size() > 1) {
              Rcpp::stop("Parameter " + (*it) + "given as vector, expecting " +
                         "a scalar.");
          }
      }
  }

  // read initial state
  Rcpp::DataFrame init(rcpp_init);

  Rcpp::CharacterVector cityNames = init.attr("row.names");

  Rcpp::CharacterVector initVariables = init.names();
  std::vector<std::vector<int> > initValues;

  for (Rcpp::CharacterVector::iterator it = initVariables.begin();
       it != initVariables.end(); it++)
    initValues.push_back(init[std::string(*it)]);

  model.Init(cityNames, initValues);

  // birth data/parameter
  if (model_params.containsElementNamed("births")) {
    // row names of the passed data frame (or vector) are time points
    Rcpp::DataFrame births(model_params["births"]);
    std::vector<double> birthTimes =
      Rcpp::as<std::vector<double> >(births.attr("row.names"));

    // check if a number/vector or data frame was passed
    switch(TYPEOF(model_params["births"])) {

      // if births parameter is a number or vector, we're assuming
      // that it's one entry per time point, equal for all cities
    case INTSXP:
    case REALSXP: {
        std::vector<unsigned int> birthVector =
          Rcpp::as<std::vector<unsigned int> >(model_params["births"]);

        // loop over all cities and set to the same births vector
        for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
             it != model.cities.end(); it++) {
          it->births.setValues(birthTimes, birthVector);
        }

        break;
      }
      // if births parameter is a list (e.g., a data frame), we're assuming
      // that it's one column per city and one row per time point
    case VECSXP: {
        // loop over all cities and set to different births vectors
        for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
             it != model.cities.end(); it++) {
          std::vector<unsigned int> birthVector =
            Rcpp::as<std::vector<unsigned int> >(births[it->name]);
          it->births.setValues(birthTimes, birthVector);
        }
        break;
      }
    }
  } else {
    Rcpp::stop("Must supply births as a parameter");
  }

  // population size data / parameter
  if (model_params.containsElementNamed("N")) {
    // row names of the passed data frame (or vector) are time points
    Rcpp::DataFrame popSize(model_params["N"]);
    std::vector<double> popSizeTimes =
      Rcpp::as<std::vector<double> >(popSize.attr("row.names"));

    // check if a number/vector or data frame was passed
    switch(TYPEOF(model_params["N"])) {

      // if popSize parameter is a number or vector, we're assuming
      // that it's one entry per time point, equal for all cities
    case INTSXP:
    case REALSXP:
      {
        std::vector<unsigned int> popSizeVector =
          Rcpp::as<std::vector<unsigned int> >(model_params["N"]);

        // loop over all cities and set to the same popSize vector
        for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
             it != model.cities.end(); it++) {
          it->popSize.setValues(popSizeTimes, popSizeVector);
        }
      }

  break;
      // if popSize parameter is a list (e.g., a data frame), we're assuming
      // that it's one column per city and one row per time point
    case VECSXP:
      {
        // loop over all cities and set to different popSize vectors
        for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
             it != model.cities.end(); it++) {
          std::vector<unsigned int> popSizeVector =
            Rcpp::as<std::vector<unsigned int> >(popSize[it->name]);
          it->popSize.setValues(popSizeTimes, popSizeVector);
        }
        break;
      }
    }
  } else {
    Rcpp::stop("Must supply N as a model parameter.");
  }

  // vaccination data / parameter
  if (model_params.containsElementNamed("vaccination")) {
    // row names of the passed data frame (or vector) are time points
    Rcpp::DataFrame vaccination(model_params["vaccination"]);
    std::vector<double> vaccinationTimes =
      Rcpp::as<std::vector<double> >(vaccination.attr("row.names"));

    // check if a number/vector or data frame was passed
    switch(TYPEOF(model_params["vaccination"])) {

      // if vaccination parameter is a number or vector, we're assuming
      // that it's one entry per time point, equal for all cities
    case INTSXP:
    case REALSXP:
      {
        std::vector<double> vaccinationVector =
          Rcpp::as<std::vector<double> >(model_params["vaccination"]);

        // loop over all cities and set to the same vaccination vector
        for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
             it != model.cities.end(); it++) {
          it->vaccination.setValues(vaccinationTimes, vaccinationVector);
        }

        break;
      }
      // if vaccination parameter is a list (e.g., a data frame), we're assuming
      // that it's one column per city and one row per time point
    case VECSXP:
      {
        // loop over all cities and set to different vaccination vectors
        for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
             it != model.cities.end(); it++) {
          std::vector<double> vaccinationVector =
            Rcpp::as<std::vector<double> >(vaccination[it->name]);
          it->vaccination.setValues(vaccinationTimes, vaccinationVector);
        }
        break;
      }
    }
  } else {
    std::vector<double> vaccinationVector(1, 0.);
    std::vector<double> vaccinationTimes(1, 0.);
    for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
         it != model.cities.end(); it++) {
      it->vaccination.setValues(vaccinationTimes, vaccinationVector);
    }
  }

  // coordinate data -- this only works as a data frame
  if (model_params.containsElementNamed("coordinates") &&
      TYPEOF(model_params["coordinates"]) == VECSXP) {
    // row names of the passed data frame (or vector) are time points
    Rcpp::DataFrame coordinates(model_params["coordinates"]);

    // loop over all cities and set to different coordinates vectors
    for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
         it != model.cities.end(); it++) {
      std::vector<double> coordinateVector =
          Rcpp::as<std::vector<double> >(coordinates[it->name]);
      it->coordinates = std::make_pair(coordinateVector[0],
                                       coordinateVector[1]);
    }
  } else {
    // set all coordinates to (0,0) -- the gravity model will not be
    // applied, and all distances taken to be 1
    for (std::vector<GravityModel::city>::iterator it = model.cities.begin();
         it != model.cities.end(); it++) {
      it->coordinates = std::make_pair(0,0);
    }
  }

  // seasonality
  if (model_params.containsElementNamed("seasonality")) {
    Rcpp::NumericVector seasonality = model_params["seasonality"];
    model.seasonality = Rcpp::as<std::vector<double> >(seasonality);
  } else {
    model.seasonality = std::vector<double>(1, .1);
  }

  Rcpp::IntegerVector seed(rcpp_seed);
  if (R_IsNA(seed[0])) {
    seed = getSeed();
  } else {
    seed = seed[0];
  }

  Simulation sim(model, nThreads, verbose);

  sim.setRecordStep(.0);

  sim.setStartTime(startTime);
  sim.setEndTime(endTime);

  sim.setTimeStep("1/26");

  bool res = sim.run(seed[0]);

  Rcpp::CharacterVector summaryStrings(rcpp_summary);
  Rcpp::NumericVector summaryStats;
  sim.getSummaryStats(summaryStrings, summaryStats);

  Rcpp::List trajectory(sim.getData()[0].size());
  Rcpp::CharacterVector nameVector(sim.getData()[0].size());

  nameVector[0] = ("sim");
  nameVector[1] = ("time");

  for (size_t i = 0; i < sim.getData()[0].size(); ++i) {
    std::vector<double> column;
    for (size_t j = 0; j < sim.getData().size(); ++j) {
      column.push_back(sim.getData()[j][i]);
    }
    trajectory[i] = column;
    if (i > 1) {
      nameVector[i] = (i % 2 == 0 ? "S." : "I.") +
        sim.getCityName((i - 2) / 2);
    }
  }

  std::cout << sim.getData()[0].size() << " " << sim.getData().size()
            << std::endl;

  trajectory.attr("names") = nameVector;

  return Rcpp::List::create(Rcpp::Named("summaryStats") = summaryStats,
                            Rcpp::Named("extinction") = !res,
                            Rcpp::Named("trajectory") = Rcpp::DataFrame(trajectory));
}
