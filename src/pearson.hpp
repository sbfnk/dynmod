// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file pearson.hpp
  \brief Header file for the pearson class
*/

#ifndef PEARSON_HPP
#define PEARSON_HPP

#include <numeric>
#include <math.h>

//! Structure to calculate the Pearson's correlation coefficient
//! between two vectors x and y
class pearson 
{

public:
  pearson(): xy(.0), x(.0), y(.0), x2(.0), y2(.0), n(0) {;}

  //! Add a data point
  /*!
    \param[in] addX x value to add to correlatoin calculatoin
    \param[in] addY y value to add to correlatoin calculatoin
   */
  void add(double addX, double addY);

  //! Calculate Pearson's correlation coefficient
  double calculate() const;
  
private:
  double xy; //!< \\sum x*y
  double x;  //!< \\sum x
  double y;  //!< \\sum y
  double x2; //!< \\sum x^2
  double y2; //!< \\sum y^2
  int n;  //!< counter
};

#endif
