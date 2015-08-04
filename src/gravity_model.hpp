// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file gravity_model.hpp
  \brief Header file for the GravityModel class
*/

#ifndef GRAVITY_MODEL_HPP
#define GRAVITY_MODEL_HPP

#include <boost/program_options.hpp>
#include <iostream>

#include "param_container.hpp"
#include "pearson.hpp"

namespace po = boost::program_options;

/*! Class for measles data
  /parameters */
class GravityModel :
    public ParamContainer
{
public:
    //! Constructor
    /*!
      \param[in] dir Default directory from which to read parameter
      files
    */
    GravityModel(std::string dir = "", unsigned int verbose = 0);

    //! Read command line parameters.
    /*!
      This reads the model command line parameters and tables for the
      measles model.

      \param[in] vm The map of command line parameters
      \param[in] startTime Start time for the model
    */
    bool InitParams(po::variables_map& vm);

    //! Data of one city
    /*!
      This contains all the parameters and variables of a city, i.e. its
      name, birth and vaccination rates, population sizes, coordinates,
      initial numbers of suseptibles and infected, as well as numbers of
      susceptibles and infected.
    */
    struct city
    {
        //! Constructor (assigns a name)
        city(std::string n) : name(n) {;}

        std::string name; //!< name of the city

        //! Births over time.
        Data<unsigned int> births,
        //! Population size over time.
            popSize;
        //! Vaccination uptake over time.
        Data<double>  vaccination;

        //! lat/long coordinates of the city
        std::pair<double, double> coordinates;

        //!< initial numbers of susceptible/infected
        std::vector<int> init;

        //! The main variables: numbers of susceptible and infected
        std::vector<int> variables;

        //! cumulative number of infected (to record infections if
        //! reporting period <1)
        unsigned int cumulativeI;

    };

    std::vector<city> cities; //!< All the cities in the data set
    //! Seasonality of transmsission (modulates beta)
    std::vector<double> seasonality;
    //! Age-specific force of infection (currently not implemented)
    std::vector<double> foi;

    //! Base directory for parameter files
    std::string baseDir;

    //! Epiphenomenological exponent
    Parameter alpha,
    //! Infection rate
        beta,
    //! Scaling of gravitational attraction with size
        tau1,
    //! Scaling of emigration with size
        tau2,
    //! Spatial coupling strength
        theta,
    //! Decay of attraction with distance
        rho;
    //! number of recovered classes before going back to susceptible
    size_t nRecovered;

    //! index of London (often needed, so may as well store it)
    int london;

    //! do we want verbose output?
    unsigned int verbose;

    //! last timestep that was calculated
    double lastTime;

    // read initial conditions from provided data file
    bool ReadInitFile(std::string fileName);

    // initialise model with given initial values
    template <class StringVector>
    void Init(StringVector cityNames,
              std::vector< std::vector<int> > initValues);
};

template <class StringVector>
void GravityModel::Init(StringVector cityNames,
                        std::vector<std::vector<int> > initValues)
{
    // loop over all cities given in the initial conditions file and
    // initialise them with the given data
    for (size_t i = 0; i < cityNames.size(); ++i)
    {
        cities.push_back(city(std::string(cityNames[i])));
        for (std::vector<std::vector<int> >::iterator it = initValues.begin();
             it != initValues.end(); it++)
            cities.back().init.push_back((*it)[i]);

        if (cityNames[i] == "London") // found London
            london = i;
    }
}

#endif
