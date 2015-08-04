// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file city_calculator.cpp
  \brief Implementation file for the city calculator

  Updates simulation data for a number of cities by a step.
*/

#include "city_calculator.hpp"

void CityCalculator::operator() ()
{

    boost::mt19937 gen(seed);

    // Define local variables
    std::vector< GravityModel::city >& cities = model.cities;

    const double& alpha = model.alpha.getValue();
    const double& tau1 = model.tau1.getValue();
    const double& theta = model.theta.getValue();

    std::vector<double> sum(cities.size(), .0);

    //! This calculates the number of new infections in one city.

    for (size_t j = lower; j < upper; ++j)
    {
        double iota;
        double lambda;

        for (size_t k = 0; k < cities.size(); ++k)
        {
            if (j != k)
            {
                sum[j] += tau2_powers[cities[k].variables[I]] /
                    (1 + distance[j][k]);
            }
        }

        double m = theta *
            pow(cities[j].popSize.getValue(time), tau1) *
            sum[j];

        if (verbose > 1)
        {
            std::cout << "m[" << j << "] = " << theta << " * "
                      << (cities[j].popSize.getValue(time))
                      << "^" << tau1 << " * " << sum[j]
                      << " = " << m << ", ";
        }

        if (m > 0)
        {
            boost::random::gamma_distribution<> gd(m);
            boost::variate_generator
                <boost::mt19937&, boost::random::gamma_distribution<> >
                var_gamma(gen, gd);
            iota = var_gamma(m);
        }
        else
        {
            iota = 0;
        }

        if (cities[j].popSize.getValue(time) == 0)
            cities[j].variables[S] = 0;

        //!< Current force of infection
        double currentFOI =
            beta_prime *
            pow(cities[j].variables[I] + iota, alpha) /
            cities[j].popSize.getValue(time);

        if (model.foi.size() > 0)
        {
            // adjust force of infection according to the average age of infection
            size_t avgAge = 0;
            if (cities[j].births.getValue(time) > 0)
            {
                avgAge = static_cast<size_t>
                    (cities[j].popSize.getValue(time) /
                     cities[j].births.getValue(time) /
                     (15 * (1-cities[j].vaccination.getValue(time))));
            }

            double modFOI;
            if (avgAge > model.foi.size())
            {
                modFOI = model.foi.back();
            }
            else
            {
                modFOI = model.foi[avgAge];
            }

            currentFOI *= modFOI;

            if (verbose > 1)
            {
                std::cout << time << " " << j <<  "," << cities[j].name << ": avg age "
                          << avgAge << ", FOI modulator: " << modFOI << std::endl;
            }
        }

        lambda = cities[j].variables[S] * (1 - exp(-currentFOI));

        if (verbose > 1)
        {
            std::cout << "lambda[" << j << "] = " << cities[j].variables[S]
                      << " * (1 - exp("
                      << beta_prime << " * ("
                      << cities[j].variables[I] << " + "
                      << iota << ")^" << alpha << ") / "
                      << cities[j].popSize.getValue(time) << ") = "
                      << lambda << std::endl;
        }

        double totalI = cities[j].variables[I] + iota;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
            var_uniform(gen, boost::uniform_real<> (0,1));
        if (totalI > 0)
        {

            unsigned int c = 0;
            double tempNewInfections = cities[j].variables[S] + 1;
            while (tempNewInfections > cities[j].variables[S] &&
                   c++ < EXIT_CONSTANT)
            {
                tempNewInfections =
                    floor(quantile
                          (boost::math::negative_binomial_distribution<>
                           (totalI, totalI/(totalI+lambda)),
                           var_uniform()
                           )
                          );
                if (tempNewInfections > cities[j].variables[S])
                {
                    tempNewInfections = cities[j].variables[S];
                }
            }

            if (tempNewInfections > cities[j].variables[S])
            {
                std::cerr << std::endl;
                std::cerr << "ERROR: Couldn't get newInfected < S in "
                          << EXIT_CONSTANT << " iterations. Bailing out."
                          << std::endl;
                std::cerr << std::endl;
                std::cerr << "Variable values:" << std::endl;
                std::cerr << "================" << std::endl;
                std::cerr << "i = " << (j) << std::endl;
                std::cerr << "S[i] = " << cities[j].variables[S] << std::endl;
                std::cerr << "I[i] = " << cities[j].variables[I] << std::endl;
                std::cerr << std::endl;
                exit(1);
            }

            newInfections[j] = static_cast<unsigned int>(tempNewInfections);

        }

        newBirths[j] = cities[j].births.getValue(time, model.lastTime);

    }
}
