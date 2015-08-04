// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file getseed.hpp
  \brief Header file for the getting a random seed
*/

#ifndef GETSEED_HPP
#define GETSEED_HPP

//! Obtain random number generator seed from /dev/urandom.
/*!
  \return A random integer
*/
int getSeed();

#endif
