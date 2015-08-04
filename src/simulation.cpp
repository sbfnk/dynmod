/*! \file simulation.cpp
  \brief Implementation file for the Simulation class
*/

#include <algorithm>
#include "simulation.hpp"

void Simulation::setTimeStep(std::string timeStep)
{
    std::vector<std::string> strings;
    boost::split(strings, timeStep,
                 boost::is_any_of("/"));
    if (strings.size() > 1)
    {
        std::istringstream s(strings[1]);
        s >> subStep;
        std::istringstream t(strings[0]);
        t >> timeStep;
    }
    else
    {
        subStep = 1;
        std::istringstream i(strings[0]);
        i >> timeStep;
    }
}

bool Simulation::run(int& seed, int simNb)
{
    // clear data
    data.clear();

    // Random generator initialised with given seed
    boost::mt19937 gen(seed);

    // uniform random generator
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
        var_uniform(gen, boost::uniform_real<> (0,1));

    // set parameter values
    model.setRandomParams(seed);

    if (verbose)
        model.Print();

    //!< for caching powers of tau2, to speed up the code
    std::vector<double> tau2_powers(1, .0);

    //! Local variable to access city data
    std::vector<GravityModel::city>& cities = model.cities;

    //! Distances between cities
    std::vector<std::vector<double> > distance
        (cities.size(),
         std::vector<double>(cities.size()));

    // convert lat/long coordinates to radians
    for (std::vector<GravityModel::city>::iterator
             it = cities.begin(); it != cities.end(); it++)
    {
        it->coordinates.first *= PI_VALUE / 180.;
        it->coordinates.second *= PI_VALUE / 180.;
    }

    // Calculate distances between cities in km

    for (size_t i = 0; i < cities.size(); ++i)
    {
        for (size_t j = i; j < cities.size(); ++j)
        {
            distance[i][j] = distance[j][i] =
                pow(
                    sqrt(pow((cities[i].coordinates.first -
                              cities[j].coordinates.first) *
                             cos((cities[i].coordinates.first +
                                  cities[j].coordinates.first) / 2), 2) +
                         pow(cities[i].coordinates.second -
                             cities[j].coordinates.second, 2)) *
                    EARTH_RADIUS,
                    model.rho.getValue()
                    );
        }
    }

    // This is for the threads later, bounds for splitting the city
    // dataset
    std::vector<size_t> bnd = bounds(nThreads, cities.size());

    size_t currentSeason = 0;

    double time = tStart; // current time

    //! This initialised the model variables with the given initial
    //! conditions.

    for (std::vector<GravityModel::city>::iterator it = cities.begin();
         it != cities.end(); it++)
    {
        it->variables.resize(model.nRecovered + 2, 0);
        it->variables[I] = it->init[0];
        if (it->init.size() == 1)
            it->variables[S] =
                static_cast<int>(floor(it->popSize.getValue(time)/30.));
        else
            it->variables[S] = it->init[1];

        it->cumulativeI = it->variables[I];
    }

    double tn = recordStep; //!< Next time the state of the system
                            //!will be recorded
    size_t currentTimeStep = 0; //!< current time step
    size_t currentSubStep = 0; //!< current sub time step

    double saveTime = time; //! Previous time step
    double recordTime = time; //! Previous recorded time step

    double nextIntroduction = 0; //! Time of next introduction

    //! Then it calculates the time of first introduction (with
    //! exponential waiting time, according to the given rate)

    if (introductions > 0)
    {
        double rnd = var_uniform();
        nextIntroduction = time - log(rnd) / introductions;
        if (verbose > 1)
        {
            std::cout << "Next introduction at " << nextIntroduction << std::endl;
        }
    }

    //! Do we need to worry about future vaccination rates?
    if (futureVaccination >= 0)
    {
        if (varyVaccination)
        {
            // if we vary vaccination, we first put the cities
            // in random order and then assign vaccination rates
            // one-by-one, adjusting as we go along to match the
            // given population vaccination rate

            //! Cities in random order
            std::vector<size_t> random_order;

            //! Size of the population to whom have already
            //! assigned a vaccination rate
            size_t doneSize = 0;
            size_t totalSize = 0; //! Total population size
            //! Current overall vaccination rate of the part of
            //! the population to which we have already assigned
            //! vaccination rates
            double currentVaccination = .0;

            //! Currently required vaccination rate for the
            //! population that has not yet been assigned a
            //! vaccination rate to match the population
            //! vaccination rate
            double currentTarget = futureVaccination;

            // Put cities in random order (and record total
            // population size)
            for (size_t i = 0; i < cities.size(); ++i)
            {
                random_order.push_back(i);
                totalSize += cities[i].popSize.getLastValue();
            }

            std::random_shuffle(random_order.begin(),random_order.end());

            // Go through cities one-by-one and assign
            // vaccination rates
            for (std::vector<size_t>::iterator it = random_order.begin();
                 it != random_order.end(); it++)
            {
                double rnd = var_uniform();
                double uptake;

                // if currentTarget, the overall vaccination
                // rate which the part of the population which
                // has not yet been assigned a vaccination rate
                // has to have for the whole population to match
                // the overall vaccination rate is >100%, we
                // assign 100% to all remaining populations
                if (currentTarget > 100)
                {
                    uptake = 100;
                }
                // if currentTarget is >50%, we assign a random
                // vaccination rate between 100% and
                // (200-2*currentTarget))
                else if (currentTarget > 50)
                {
                    uptake = 100 - (200-2*currentTarget)*rnd;
                }
                // if currentTarget is <50%, we assign a random
                // vaccination rate between 0% and
                // (2*currentTarget))
                else
                {
                    uptake = (2*currentTarget)*rnd;
                }

                cities[*it].vaccination.addLastValue(uptake);

                doneSize += cities[*it].popSize.getLastValue();
                currentVaccination +=
                    uptake * cities[*it].popSize.getLastValue();

                // update currentTarget given what has been
                // assigned to the city currently under consideration
                currentTarget =
                    (totalSize * futureVaccination - currentVaccination) /
                    (totalSize - doneSize + .0);
            }
        }
        else // vaccination rates are not being varied, every
            // city gets the same (futureVaccination)
            // vaccination rate
        {
            for (std::vector<GravityModel::city>::iterator it =
                     cities.begin(); it != cities.end(); it++)
            {
                it->vaccination.addLastValue(futureVaccination);
            }
        }
    }

    // save initial state
    // add an empty row
    data.push_back(std::vector<double>());
    // add data to row
    data.back().push_back(simNb);
    data.back().push_back(recordTime);
    for (std::vector<GravityModel::city>::iterator it = cities.begin();
         it != cities.end(); it++)
    {
        data.back().push_back(it->variables[S]);
        data.back().push_back(it->cumulativeI);
        it->cumulativeI = 0;
    }

    // Then it runs the simulation until the given end time
    while (time <= tEnd)
    {
        // update timestep
        ++currentSubStep;
        if (currentSubStep == subStep)
        {
            currentTimeStep += timeStep;
            currentSubStep = 0;
        }
        saveTime = time;
        time = tStart + currentTimeStep +
            currentSubStep*timeStep/(subStep + .0);

        // do we need to advance the season?
        double intpart;
        // if the real part of time has incrased (new timestep), we
        // set currentSeason to zero
        if (modf(model.lastTime, &intpart) <
            (currentSeason / (model.seasonality.size() + .0)))
            currentSeason = 0;

        // advance currentSeason until it matches the real part of time
        while (fabs(((currentSeason+1) / (model.seasonality.size() + .0)) -
                    modf(model.lastTime, &intpart)) <
               fabs((((currentSeason) / (model.seasonality.size() + .0)) -
                     modf(model.lastTime, &intpart))))
            ++currentSeason;

        double beta_prime = model.beta.getValue() *
            model.seasonality[currentSeason];

        // Do we have a new introduction?
        if (introductions > 0)
        {
            unsigned int pop_sum = 0;
            while (nextIntroduction < time)
            {
                // introduce infection proportional to susceptible city size
                if (pop_sum == 0)
                {
                    for (std::vector<GravityModel::city>::iterator it =
                             cities.begin(); it != cities.end(); it++)
                    {
                        pop_sum += it->popSize.getValue(time);
                    }
                }

                // choose random city to get infected, proportional to city size
                unsigned int end_sum =
                    static_cast<int>(pop_sum * var_uniform());
                std::vector<GravityModel::city>::iterator it = cities.begin();
                unsigned int running_pop_sum = it->popSize.getValue(time);

                while (running_pop_sum < end_sum)
                {
                    it++;
                    running_pop_sum += it->popSize.getValue(time);
                }

                // do we have any susecptibles in the chosen city? Or,
                // if the overall population is zero, we skip this
                // introduction
                if (it->variables[S] > 0 || pop_sum == 0)
                {
                    // if we have any susceptibles in the chosen city,
                    // we introduce infection there
                    if (it->variables[S] > 0)
                    {
                        if (verbose > 1)
                        {
                            std::cout << "Infection introduced into "
                                      << it->name << std::endl;
                        }
                        --(it->variables[S]);
                        ++(it->variables[I]);
                        ++(it->cumulativeI);
                    }

                    // determine time of next introduction
                    nextIntroduction -= log(var_uniform()) / introductions;
                    if (verbose > 1)
                    {
                        std::cout << "Next introduction at " << nextIntroduction
                                  << std::endl;
                    }
                }
            }
        }

        // Pre-cache powers to \tau_2. These take a long time to compute
        // so we want to calculate every one of them only once.

        for (std::vector<GravityModel::city>::iterator it = cities.begin();
             it != cities.end(); it++)
        {
            if (it->variables[I] > 0)
            {
                if (tau2_powers.size() <
                    static_cast<unsigned int>(((it->variables[I]) + 1)))
                {
                    tau2_powers.resize((it->variables[I]) + 1, .0);
                }
                if (tau2_powers[it->variables[I]] == 0)
                {
                    tau2_powers[it->variables[I]] =
                        pow(it->variables[I], model.tau2.getValue());
                }
            }
        }

        //! Temporary variable to store the number of new births in
        //! every city
        std::vector<unsigned int> newBirths(cities.size(), 0);

        //! Temporary variable to store the number of new infections in
        //! every city
        std::vector<unsigned int> newInfections(cities.size(), 0);

        //! Group of threads for parallelising city updates
        boost::thread_group threads;

        // Launch threads
        for (size_t i = 0; i < bnd.size() - 1; ++i)
        {
            CityCalculator c(bnd[i],
                             bnd[i+1],
                             gen(),
                             model,
                             newInfections,
                             newBirths,
                             distance,
                             time,
                             beta_prime,
                             tau2_powers,
                             verbose);

            if (nThreads > 1)
                threads.add_thread(new boost::thread(c));
            else
                c();
        }

        // Join all threads
        if (nThreads > 1)
            threads.join_all();

        // update variables from what's been reported back from the
        // threads
        for (size_t i = 0; i < cities.size(); ++i)
        {
            if (verbose > 1)
            {
                std::cout << "Updating susceptibles: "
                          << cities[i].variables[S] << " + "
                          << newBirths[i] << " - "
                          << newInfections[i];
            }
            cities[i].variables[S] =
                cities[i].variables[S] + newBirths[i] - newInfections[i];
            if (verbose > 1)
            {
                std::cout << " = " << cities[i].variables[S] << std::endl;
            }
            if (model.nRecovered > 0)
            {
                cities[i].variables[S] +=
                    cities[i].variables[model.nRecovered + 1];
            }
            for (size_t j = model.nRecovered + 1; j > I; --j)
            {
                cities[i].variables[j] = cities[i].variables[j-1];
            }
            cities[i].variables[I] = newInfections[i];
            cities[i].cumulativeI += newInfections[i];

            if (verbose > 1)
            {
                std::cout << "t = " << time;
                std::cout << ", I[" << i << "] = " << cities[i].variables[I];
                std::cout << ", S[" << i << "] = " << cities[i].variables[S];
                std::cout << ", cumulativeI[" << i << "] = "
                          << cities[i].cumulativeI;
                std::cout << ", Inew[" << i << "] = " << newInfections[i];
                std::cout << ", B[" << i << "] = " << newBirths[i];
                std::cout << std::endl;
            }
        }

        // save state of the system if the current time is greater than
        // or equal to the calculated time this should happen next
        if (time >= tn && recordStep >= 0)
        {
            // add an empty row
            data.push_back(std::vector<double>());
            // add data to row
            data.back().push_back(simNb);
            data.back().push_back(recordTime);
            for (std::vector<GravityModel::city>::iterator it = cities.begin();
                 it != cities.end(); it++)
            {
                data.back().push_back(it->variables[S]);
                data.back().push_back(it->cumulativeI);
                it->cumulativeI = 0;
            }

            recordTime = time;
            tn += recordStep;
        }
        model.lastTime = time;
    }

    // // if last time step is greater than the last recorded time
    // if (model.lastTime > data.back()[1])
    //   {
    //     // add final state of system to data
    //     data.push_back(std::vector<double>());
    //     data.back().push_back(simNb);
    //     data.back().push_back(time);
    //     for (std::vector<GravityModel::city>::iterator it = cities.begin();
    //          it != cities.end(); it++)
    //       {
    //         data.back().push_back(it->variables[S]);
    //         data.back().push_back(it->cumulativeI);
    //         it->cumulativeI = 0;
    //       }
    //   }

    double sum = 0;
    for (size_t i = 0; i < cities.size(); ++i)
    {
        sum += data.back()[3+2*i];
    }
    if (sum > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Simulation::print(std::ostream &os) const
{
    for (std::vector<std::vector<double> >::const_iterator it = data.begin();
         it != data.end(); it++)
    {
        bool firstColumn = true;
        for (std::vector<double>::const_iterator it2 = it->begin();
             it2 != it->end(); it2++)
        {
            if (firstColumn)
            {
                firstColumn = false;
            }
            else
            {
                os << ",";
            }
            os << *it2;
        }
        os << std::endl;
    }
}

std::string Simulation::getCityName(size_t i) const
{
  return model.cities[i].name;
}

const std::vector<std::vector<double> >& Simulation::getData() const
{
  return data;
}

std::ostream& operator<<(std::ostream& os, const Simulation& s)
{
    s.print(os);
    return os;
}

std::vector<size_t> bounds(size_t parts, size_t mem)
{
    std::vector<size_t>bnd;
    size_t delta = mem / parts;
    size_t remainder = mem % parts;
    size_t N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (size_t i = 0; i < parts; ++i)
    {
        N2 = N1 + delta;
        if (i == parts - 1)
            N2 += remainder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}
