// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file gravity_model.cpp
  \brief Implementation file for the GravityModel class
*/

#include "gravity_model.hpp"

GravityModel::GravityModel(std::string dir, unsigned int verbose) :
    ParamContainer("gravity"),
    baseDir(dir),
    nRecovered(0),
    london(-1),
    verbose(verbose)
{
    //! This defines all the command line options.
    options.add_options()
        ("births,b",
         po::value<std::string>()->default_value
         (baseDir + "births_4494.csv"),
         "Table of birth rates")
        ("pop-sizes,p",
         po::value<std::string>()->default_value
         (baseDir + "popsizes_4494.csv"),
         "Table of poulation sizes")
        ("seasonality,s",
         po::value<std::string>()->default_value
         (baseDir + "seasonality.csv"),
         "Seasonality table")
        ("xytable,x",
         po::value<std::string>()->default_value
         (baseDir + "coordinates_4494.csv"),
         "Table of xy coordinates")
        ("init,i",
         po::value<std::string>()->default_value
         (baseDir + "init_4494.csv"),
         "Initial conditions")
        ("vaccination,n",
         po::value<std::string>()->default_value
         (baseDir + "vaccination_4494.csv"),
         "Vaccination coverage")
        ("params",
         po::value<std::string>()->default_value
         (baseDir + "params"),
         "Params file")
        ("foi",
         po::value<std::string>()->default_value
         (baseDir + "foi"),
         "Force of infection")
        ("nrecovered",
         po::value<unsigned int>()->default_value
         (0),
         "Number of recovered classes (0 for lifelong recovery)")
        ;

    //! Then it fills the params vector.
    params.push_back(ParamInfo("alpha", "Epiphenomenological exponent", &alpha));
    params.push_back(ParamInfo("beta", "Infection rate", &beta));
    params.push_back(ParamInfo("tau1", "Scaling of attraction with size", &tau1));
    params.push_back(ParamInfo("tau2", "Scaling of emigration with size", &tau2));
    params.push_back(ParamInfo("theta", "Spatial coupling strength", &theta, true));
    params.push_back(ParamInfo("rho", "Decay of attraction", &rho));

    //! Then it adds the parameters in the params vector to the command
    //! line options.
    for (std::vector<ParamInfo>::iterator it = params.begin();
         it != params.end(); it++)
    {
        options.add_options()
            ((it->option).c_str(), po::value<double>(),
             it->description.c_str())
            ((it->option+"_low").c_str(), po::value<double>(),
             it->description.c_str())
            ((it->option+"_high").c_str(), po::value<double>(),
             it->description.c_str());
    }
}

bool GravityModel::InitParams(po::variables_map &vm)
{
    //! This reads model parameters and data files from the given
    //! command line parameters.
    if (vm.count("init"))
    {
        if (cities.size() == 0)
        {
            if (!ReadInitFile(vm["init"].as<std::string>()))
            {
                std::cerr << "ERROR: error reading initial conditions "
                          << std::endl;
                return false;
            }
        }
    }
    else
    {
        std::cerr << "ERROR: No initial conditions" << std::endl;
        return false;
    }
    if (verbose > 0) {
        std::cout << "Read initial conditions for " << cities.size()
                  << " cities from " << vm["init"].as<std::string>()
                  << std::endl;
    }

    // read birth rate data
    if (vm.count("births"))
    {
        std::vector<std::vector<unsigned int> > tempBirths;
        std::vector<std::string> birthCities;
        std::vector<double> birthTimes;
        // read birth rates from provided data file
        if (ReadTable(vm["births"].as<std::string>(), birthCities,
                      birthTimes, tempBirths, true))
        {
            // assign birth rates to city variables
            size_t cities_read(0);
            for (size_t i = 0; i < birthCities.size(); ++i)
            {
                std::vector<city>::iterator it = cities.begin();
                while (it != cities.end() && it->name != birthCities[i])
                    ++it;

                if (it != cities.end())
                {
                    it->births.setValues(birthTimes, tempBirths[i]);
                    // we add a last value for interpolation
                    it->births.addLastValue(tempBirths[i].back());
                    ++cities_read;
                }
            }
            if (verbose > 0)
            {
                std::cout << "Read birth statistics for "
                          << cities_read << " cities from "
                          << vm["births"].as<std::string>() << std::endl;
            }

        }
        else // ReadTable failed
        {
            std::cerr << "ERROR: error reading birth statistics from "
                      << vm["births"].as<std::string>() << std::endl;
            return false;
        }
    }
    else
    {
        std::cerr << "ERROR: No births file" << std::endl;
        return false;
    }

    // read population sizes data
    if (vm.count("pop-sizes"))
    {
        std::vector<std::vector<unsigned int> > tempPopSizes;
        std::vector<std::string> popSizeCities;
        std::vector<double> popSizeTimes;
        // read population sizes from provided data file
        size_t pop_read(0);
        if (ReadTable(vm["pop-sizes"].as<std::string>(), popSizeCities,
                      popSizeTimes, tempPopSizes, true))
        {
            // assign population sizes to city variables
            for (size_t i = 0; i < popSizeCities.size(); ++i)
            {
                std::vector<city>::iterator it = cities.begin();
                while (it != cities.end() && it->name != popSizeCities[i])
                    ++it;
                if (it != cities.end())
                {
                    it->popSize.setValues(popSizeTimes, tempPopSizes[i]);
                    ++pop_read;
                }
            }
            if (verbose > 0)
            {
                std::cout << "Read population sizes for "
                          << popSizeCities.size() << " cities from "
                          << vm["pop-sizes"].as<std::string>() << std::endl;
            }
        }
        else // ReadTable failed
        {
            std::cerr << "ERROR: error reading population sizes from "
                      << vm["pop-sizes"].as<std::string>() << std::endl;
            return false;
        }
    }
    else
    {
        std::cerr << "ERROR: No population sizes file" << std::endl;
        return false;
    }

    // read seasonality data
    if (vm.count("seasonality"))
    {
        std::vector<std::vector<double> > tempSeasonality;
        std::vector<double> times;
        std::vector<int> doubleWeeks;
        // read seasonality from provided data file
        if (ReadTable(vm["seasonality"].as<std::string>(), doubleWeeks,
                      times, tempSeasonality))
        {
            if (verbose > 0)
            {
                std::cout << "Read " << tempSeasonality.front().size()
                          << " seasonality values from "
                          << vm["seasonality"].as<std::string>() << std::endl;
            }
            seasonality = tempSeasonality.front();
        }
        else // ReadTable failed
        {
            std::cerr << "ERROR: error reading seasonality from "
                      << vm["seasonality"].as<std::string>() << std::endl;
            return false;
        }

    }
    else
    {
        std::cerr << "ERROR: No seasonality file" << std::endl;
        return false;
    }

    // read lat/long coordinates data
    if (vm.count("xytable"))
    {
        std::vector<std::vector<double> > tempCoordinates;
        std::vector<std::string> columnHeaders;
        std::vector<std::string> rowHeaders;
        // read dat/long coordinates from provided data file
        if (ReadTable(vm["xytable"].as<std::string>(), columnHeaders,
                      rowHeaders, tempCoordinates, true))
        {
            // assign lat/long coordinates to city variables
            size_t coordinates_read(0);
            for (size_t i = 0; i < columnHeaders.size(); ++i)
            {
                std::vector<city>::iterator it = cities.begin();
                while (it != cities.end() && it->name != columnHeaders[i])
                    ++it;

                if (it != cities.end())
                {
                    it->coordinates =
                        std::make_pair(tempCoordinates[i][0],
                                       tempCoordinates[i][1]);
                    ++coordinates_read;
                }
            }

            if (verbose > 0)
            {
                std::cout << "Read " << coordinates_read
                          << " coordinates from "
                          << vm["xytable"].as<std::string>() << std::endl;
            }
        }
        else // ReadTable failed;
        {
            std::cerr << "ERROR: error reading coordinates from "
                      << vm["xytable"].as<std::string>() << std::endl;

            return false;
        }
    }
    else
    {
        std::cerr << "ERROR: No coordinate file" << std::endl;
        return false;
    }

    // read vaccination data
    if (vm.count("vaccination"))
    {
        std::vector<std::vector<double> > tempVaccination;
        std::vector<std::string> vaccinationCities;
        std::vector<double> vaccinationTimes;
        // read vaccination rates from provided data file
        size_t vaccination_read(0);
        if (ReadTable(vm["vaccination"].as<std::string>(),
                      vaccinationCities, vaccinationTimes,
                      tempVaccination, true))
        {
            // assign vaccination rates to city variables
            for (size_t i = 0; i < vaccinationCities.size(); ++i)
            {
                std::vector<city>::iterator it = cities.begin();
                while (it != cities.end() && it->name != vaccinationCities[i])
                    ++it;

                if (it != cities.end())
                {
                    it->vaccination.setValues(vaccinationTimes,
                                              tempVaccination[i]);
                    ++vaccination_read;
                }

            }
            if (verbose > 0)
            {
                std::cout << "Read vaccination data for "
                          << vaccinationCities.size() << " cities from "
                          << vm["xytable"].as<std::string>() << std::endl;
            }
        }
        else // ReadTable failed
        {
            std::cerr << "ERROR: error reading vaccination data from "
                      << vm["vaccination"].as<std::string>() << std::endl;
            return false;
        }
    }
    else
    {
        std::cerr << "WARNING: No vaccination file" << std::endl;
        return false;
    }

    // read parameter file
    if (vm.count("params"))
    {
        std::ifstream in(vm["params"].as<std::string>().c_str());
        if (!in.is_open())
        {
            std::cerr << "Could not open params file "
                      << vm["params"].as<std::string>() << std::endl;
            return false;
        }
        else
        {
            // store parameters in vm variable -- this can later be read
            // out by ParamContainer::ReadParams
            po::store(po::parse_config_file(in, options), vm);
            po::notify(vm);
            in.close();
            if (verbose > 0)
            {
                std::cout << "Read parameters from "
                          << vm["params"].as<std::string>() << std::endl;
            }
        }
    }

    // read age-specific force of infection data
    if (vm.count("foi"))
    {
        std::vector<std::vector<double> > tempFOI;
        std::vector<unsigned int> foiAges;
        std::vector<double> foiTimes;
        // read age-specific force of infection from provided data file
        if (ReadTable(vm["foi"].as<std::string>(), foiAges,
                      foiTimes, tempFOI))
        {
            foi.clear();
            unsigned int currentAge = 0;
            double currentFOI = .0;

            for (size_t j = 0; j < foiAges.size(); ++j)
            {
                while (currentAge < foiAges[j])
                {
                    foi.push_back(currentFOI);
                    ++currentAge;
                }
                currentFOI = tempFOI.front()[j];
            }

            if (verbose > 0)
            {
                std::cout << "Read " << tempFOI.front().size()
                          << " age-specific forces of infection from "
                          << vm["foi"].as<std::string>() << std::endl;
            }

        }
        else
        {
            std::cerr << "WARNING: Error opening foi file" << std::endl;
        }
    }

    // assign number of recovered classes
    if (vm.count("nrecovered"))
    {
        nRecovered = vm["nrecovered"].as<unsigned int>();
    }
    return ParamContainer::ReadParams(vm);
}

bool GravityModel::ReadInitFile(std::string fileName)
{
    std::vector<std::string> variables;
    std::vector<std::string> cityNames;
    std::vector<std::vector<int> > tempInit;
    // read initial conditions from provided data file
    if (ReadTable(fileName, variables, cityNames, tempInit, true))
    {
        Init(cityNames, tempInit);
    }
    else // ReadTable failed
    {
        std::cerr << "ERROR: error reading initial contiions from "
                  << fileName << std::endl;
        return false;
    }

    return true;
}
