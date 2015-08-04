// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file parameter.hpp
  \brief Header file for the Parameter class
*/

#ifndef PARAMETER_HPP
#define PARAMETER_HPP

//! Class holding a parameter with limits (so we can shoose randomly,
//! if we want)
struct Parameter
{
  //! Default constructor
  Parameter() :
    value(0), limits(std::make_pair(0,0))
  {;}
  //! Constructor
  /*!
    \param[in] v Value of the parameter (limits assumed to be fixed)
   */  
  Parameter(double v) :
    value(v), limits(std::make_pair(v,v))
  {;}
  //! Constructor
  /*!
    \param[in] v Value of the parameter
    \param[in] l Upper and lower limits of the parameter
  */  
  Parameter(double v, std::pair<double, double> l) :
    value(v), limits(l)
  {;}

  //! Accessor for the parameter value
  double getValue() const { return value; }

  //! Set the parameter to a random value
  /*!
    \param[in] r Random number between 0 and 1
   */
  void setRandomValue(double r) 
  { value = limits.first + r * (limits.second - limits.first); }

  //! Set the parameter to a fixed value
  /*!
    \param[in] m Value to set the parameter to
  */
  void setFixedValue(double m) 
  { value = limits.first = limits.second = m; }

  //! Set the lower limit of the parameter.
  void setLowerLimit(double m){ limits.first = m;}
  //! Set the upper limit of the parameter.
  void setUpperLimit(double m) { limits.second = m;}

private:

  double value; //!< current value of the parameter
  //!< Upper and lower limits if parameter is varied
  std::pair<double, double> limits; 

};

#endif
