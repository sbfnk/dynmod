// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file param_container.cpp
  \brief Implementation file for the ParamContainer class
*/

#include <iostream>

#include "param_container.hpp"

bool ParamContainer::ReadParams(po::variables_map  &vm) 
{
  for (std::vector<ParamInfo>::iterator it = params.begin();
       it != params.end(); it++) 
    {
      if (vm.count(it->option)) 
        {
          // command line parameter has been specified, assign to
          // model variable 
          if (it->logTransform)
            it->param->setFixedValue(exp(vm[it->option].as<double>()));
          else
            it->param->setFixedValue(vm[it->option].as<double>());
        }
      if (vm.count(it->option + "_low"))
        {
          if (it->logTransform)
            it->param->setLowerLimit(vm[it->option + "_low"].as<double>());
          else
            it->param->setLowerLimit(exp(vm[it->option + "_low"].as<double>()));
        }
      if (vm.count(it->option + "_high"))
        {
          if (it->logTransform)
            it->param->setUpperLimit(vm[it->option + "_high"].as<double>());
          else
            it->param->setUpperLimit(exp(vm[it->option + "_high"].as<double>()));
        }
    }
  return true;
}

int ParamContainer::ReadInputValues(std::vector<std::string> readParams, 
                                    std::vector<double> values)
{

  int seed = 0; //!< random seed
  size_t index = 0; //!< index in values vector

  if (values.size() > readParams.size()) // first entry of values vector is seed
    {
      seed = static_cast<int>(values[0]);
      ++index;
    }

  for (size_t i = 0; i < readParams.size(); ++i)
    {
      bool found = false;
      size_t paramIndex = 0;

      while (!found && paramIndex < params.size()) 
        {
          if (params[paramIndex].option == readParams[i])
            // we've found the option
            {
              found = true;
              if (params[paramIndex].logTransform)
                values[index] = exp(values[index]);
              params[paramIndex].param->setFixedValue(values[index]);
            }

          ++paramIndex;
        }

      if (!found)
        {
          std::cerr << "WARNING: Parameter " << readParams[index]
                    << " not specified in model." << std::endl;
        }
      ++index;
    }

  return seed;
}

void ParamContainer::setRandomParams(int& seed)
{
  // Random generator initialised with given seed
  boost::mt19937 gen(seed);

  // uniform random generator
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    var_uniform(gen, boost::uniform_real<> (0,1));

  for (std::vector<ParamInfo>::iterator it = params.begin();
       it != params.end(); it++)
    it->param->setRandomValue(var_uniform());
}

bool ParamContainer::setParam(std::string name, double value)
{
  for (std::vector<ParamInfo>::iterator it = params.begin();
       it != params.end(); it++)
    if (it->option == name)
      // we've found the option
      {
        if (it->logTransform)
          it->param->setFixedValue(exp(value));
        else
          it->param->setFixedValue(value);
        
        return true;
      }
  
  return false;
}

void ParamContainer::Print() const 
{
  for (std::vector<ParamInfo>::const_iterator it = params.begin();
       it != params.end(); it++) 
    {
      std::cout << it->option << ": " << it->param->getValue() << std::endl;
    }
}
