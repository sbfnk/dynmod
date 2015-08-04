// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file city_calculator.hpp
  \brief Header file for the city calculator

  Updates simulation data for a number of cities by a step.
*/

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>

#include "gravity_model.hpp"

#ifndef CITY_CALCULATOR_HPP
#define CITY_CALCULATOR_HPP

//! exit if no suitable number of new infected can be drawn within
//! this many attempts 
#define EXIT_CONSTANT 10000

//! Disease classes (susceptible, infected)
enum Classes {S,I};

//! Class for calculating updates of a range of cities
/*
  This takes a range of cities and updates their number of infected
  and births by one time step. It defines an operator () to use with
  boost::thread multithreading.
 */
class CityCalculator 
{

public:

  //! Constructor
  CityCalculator
  (size_t lower, 
   size_t upper,
   int seed,
   GravityModel& model,
   std::vector<unsigned int>& newInfections,
   std::vector<unsigned int>& newBirths,
   std::vector<std::vector<double> > const& distance,
   double time,
   double beta_prime,
   std::vector<double> const& tau2_powers,
   unsigned int verbose) :
    lower(lower), 
    upper(upper),
    seed(seed),
    model(model),
    newInfections(newInfections),
    newBirths(newBirths),
    distance(distance),
    time(time),
    beta_prime(beta_prime),
    tau2_powers(tau2_powers),
    verbose(verbose)
  {;}
    
  //! operator for boost threading
  void operator() ();

private:

  size_t
  //!< Lower bound of range of cities to update in this thread
    lower, 
  //!< Upper bound of range of cities to update in this thread
    upper;
  //!< Random seed
  int seed;
  //!< Model to run (defining parameters, variables etc.)
  GravityModel& model;
  //!< Container to store new infections in each city
  std::vector<unsigned int>& newInfections;
  //!< Container to store new births in each city
  std::vector<unsigned int>& newBirths;
  //!< The distances between cities (assumed to have been calculated
  //! before)
  std::vector<std::vector<double> > const& distance;
  //!< Current time
  double time;
  //!< Current beta_prime (modified by seasonality)
  double beta_prime;
  //! Pre-cached powers of tau2
  std::vector<double> const& tau2_powers;
  //! Verbosity level
  unsigned int verbose;

};

#endif
