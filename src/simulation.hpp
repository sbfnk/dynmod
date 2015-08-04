// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file simulation.hpp
  \brief Header file for the Simulation class
*/

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#define PI_VALUE 3.1416 //!< pi
#define EARTH_RADIUS 6371 //!< Earth radius in km

#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>

#include <ostream>
#include <cstring>

#include "gravity_model.hpp"
#include "city_calculator.hpp"

//! Simulation parameters
struct Simulation
{

public:

  //! Constructor
  /*!
    \param[in] model The model to run
    \param[in] nThreads Number of threads for parallelisation
    \param[in] verbose Verbosity level
  */
  Simulation(GravityModel& model,
             unsigned int nThreads = 1,
             unsigned int verbose = 0) :
    recordStep(1.),
    tStart(0.),
    tEnd(100.),
    timeStep(1),
    subStep(0),
    introductions(0),
    varyVaccination(false),
    verbose(verbose),
    nThreads(nThreads),
    model(model)
  { model.lastTime = tStart; }

  //! Set time step
  /*!
    \param[in] timeStep Timestep of simulating -- this can be given as
    a fraction (e.g., 1/26), this is then interpreted as the number of
    sub steps between time steps.
  */
  void setTimeStep(std::string timeStep);

  //! Run a simulation of the specified model
  /*!
    \param[out] out Output buffer
    \param[in] seed Random generator seed
    \param[in] simNb Simulation number
    \return true if the run was successful, false otherwise
  */
  bool run(int& seed, int simNb = 0);

  // Accessor functions
  void setRecordStep(double s) { recordStep = s; }
  void setStartTime(double t) { tStart = t; model.lastTime = t; }
  void setEndTime(double t) { tEnd = t; }
  void setIntroductionFrequency(double f) { introductions = f; }
  void setFutureVaccinationRate(double r) { futureVaccination = r; }
  void setVaryVaccination() { varyVaccination = true; }

  //! Calculate summary stats
  /*!
    \param[in] stats A vector of summary statistics to calculate
    \param[out] statValues A vector of calculated summary statistics
    (note that one statistic might calculate several values)
  */
  template <class stringVector, class numberVector>
  void getSummaryStats
  (stringVector& stats, numberVector& statValues) const;

  const std::vector<std::vector<double> >& getData() const;
  std::string getCityName(size_t i) const;
  //! Write data to output stream
  void print(std::ostream &os) const;

private:

  //! record (write) data at every nth timestep (0 for always)
  double recordStep;
  double tStart;     //!< starting time
  double tEnd;     //!< end time

  //! length of one time step (this should match data to fit)
  size_t timeStep;
  size_t subStep;     //!< number of sub-steps between timesteps

  //! frequency of introductions (as Poisson process)
  double introductions;
  //! Vaccination coverage after data finishes
  double futureVaccination;

  //! Variation in vaccination coverage
  //! true if vaccination coverage varies
  bool varyVaccination;

  unsigned int verbose; //!< Verbosity level

  unsigned int nThreads; //!< Number of concurrent threads

  GravityModel& model; //!< Model to simulate

  std::vector<std::vector<double> > data; //!< Simulated data

};

//! Stream operator for Simulation class (for printing the data)
std::ostream& operator<<(std::ostream& os, const Simulation& s);

//! Extract statistics: mean, standard deviation, max, min, stdev divided by mean
/*!
  \param[in] data Vector of doubles on which to calculate statistics
  \param[in] stats String of stats to extract (x=max, i=min, m=mean, v=standard deviation, d=stdev divided by mean)
  \param[out] statValues Vector of stats
*/
template <class T>
void extractStats(std::vector<double> const& data, std::string stats,
                  T& statValues);

//! Split an integer equal parts
/*!
  \param[in] parts The number of parts.
  \param[in] mem The integer to split.
  \return A vector of split points.
*/
std::vector<size_t> bounds(size_t parts, size_t mem);

template <class stringVector, class numberVector>
void Simulation::getSummaryStats(stringVector& stats, numberVector& statValues) const
{
  for (size_t i = 0; i < stats.size(); ++i)
    {
      std::string str(stats[i]);
      if (str.substr(0, 5) == "space")
        {

          // calculate correlation coefficient and likelihood
          std::vector<pearson> acc_correlations(model.cities.size());
          std::vector<double> correlations;

          if (model.london >= 0) // London exists in the data
            {
              // variable number of the infected compartment of London
              // in the data (first two variables are simulation run and
              // time, of the London variables the first one is
              // susceptibles, hence we add 2 and 1 to get the infected
              // compartment of London)
              size_t london_inf = (model.london * 2) + 2 + 1;

              for (std::vector<std::vector<double> >::const_iterator it =
                     data.begin(); it != data.end(); it++)
                {
		  for (size_t i = 0; i < model.cities.size(); ++i)
		    {
		      // infected compartment of city i in data
		      size_t city_inf = (i*2) + 2;
		      acc_correlations[i].add((*it)[city_inf],
					      (*it)[london_inf]);
		    }
                }

              // calculate acc_correlations
              for (size_t i = 0; i < model.cities.size(); ++i)
                {
                  double correlation = acc_correlations[i].calculate();
                  // test if a number
                  if (correlation == correlation)
                    correlations.push_back(correlation);
                }

              extractStats(correlations, str.substr(5), statValues);
            }
          else // didn't find London
            {
              for (size_t j = 5; j < stats.size(); ++j)
                statValues.push_back(-1);
            }
        }
      else if (str.substr(0,6) == "london")
        {
          if (model.london >= 0) // London exists in the data
            {
              // variable number of the infected compartment of London
              // in the data (first two variables are simulation run and
              // time, of the London variables the first one is
              // susceptibles, hence we add 2 and 1 to get the infected
              // compartment of London)
              size_t london_inf = (model.london * 2) + 2 + 1;

              std::vector<double> london_data;

              for (std::vector<std::vector<double> >::const_iterator it =
                     data.begin(); it != data.end(); it++)
                  london_data.push_back((*it)[london_inf]);

              extractStats(london_data, str.substr(6), statValues);
            }
          else // didn't find London
            {
              for (size_t j = 5; j < stats.size(); ++j)
                statValues.push_back(-1);
            }
        }
      else if (str.substr(0, 8) == "national")
        {
          std::vector<double> national_sums;

          for (std::vector<std::vector<double> >::const_iterator it =
                 data.begin(); it != data.end(); it++)
            {
              double sum = 0;
	      for (size_t i = 0; i < model.cities.size(); ++i)
	      {
		// infected compartment of city i in data
		size_t city_inf = (i*2) + 3;
		std::cout << (*it)[city_inf] << std::endl;
		sum += (*it)[city_inf];
	      }
              national_sums.push_back(sum);
            }

          extractStats(national_sums, str.substr(8), statValues);
        }
    }
}

template <class T>
void extractStats(std::vector<double> const& data,
                  std::string stats,
                  T& statValues)
{

  if (data.size() == 0) return;

  double sum =
    std::accumulate(data.begin(), data.end(), 0.0);

  double mean = sum / data.size();

  std::vector<double> diff(data.size());
  std::transform(data.begin(), data.end(), diff.begin(),
                 std::bind2nd(std::minus<double>(), mean));

  double sq_sum =
    std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / data.size());

//  auto minmax = std::minmax_element(data.begin(), data.end());
  double min = *std::min_element(data.begin(), data.end());
  std::vector<double>::const_iterator maxIt =
      std::max_element(data.begin(), data.end());
  double max = *maxIt;
  size_t maxPos = std::distance(data.begin(), maxIt);

  for (size_t j = 0; j < stats.size(); ++j) {
    if (stats[j] == 'm')
      statValues.push_back(mean);
    else if (stats[j] == 'v')
      statValues.push_back(stdev);
    else if (stats[j] == 'd')
      {
        double returnVal = 0.;
        if (stdev > 0)
          returnVal = mean / stdev;
        statValues.push_back(returnVal);
      }
    else if (stats[j] == 'i')
      statValues.push_back(min);
    else if (stats[j] == 'x')
      statValues.push_back(max);
    else if (stats[j] == 't')
	statValues.push_back(maxPos);
    else if (stats[j] == 'S')
      statValues.push_back(sum);
    else if (stats[j] == 's')
      {
        size_t index = 0;
        size_t seasonIndex = 0;
        std::vector<size_t> nEntries(26, 0);
        std::vector<double> seasonalAverages(26, 0);
        while (index < data.size()) {
          seasonalAverages[seasonIndex] += data[index];
          ++nEntries[seasonIndex];
          ++index;
          ++seasonIndex;
          if (seasonIndex == nEntries.size())
            seasonIndex = 0;
        }
        for (size_t i = 0; i < nEntries.size(); ++i) {
          seasonalAverages[i] /= nEntries[i];
          statValues.push_back(seasonalAverages[i]);
        }
      }
  }

}

#endif
