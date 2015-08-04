// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file pearson.cpp
  \brief Implementation file for the pearson class
*/

#include "pearson.hpp"

void pearson::add(double addX, double addY)
{
  x += addX;
  y += addY;
  xy += addX * addY;
  x2 += pow(addX, 2);
  y2 += pow(addY, 2);
  n++;
}

double pearson::calculate() const
{
  return (n * xy - x * y) / 
    (sqrt(n * x2 - pow(x, 2)) * sqrt(n * y2 - pow(y, 2)));
}

