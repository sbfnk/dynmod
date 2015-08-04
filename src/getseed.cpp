// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file getseed.cpp
  \brief Implementation of getting a random seed
*/

#include <fstream>

int getSeed()
{
  std::ifstream rand("/dev/urandom");
  char tmp[sizeof(int)];
  rand.read(tmp,sizeof(int));
  rand.close();
  int* number = reinterpret_cast<int*>(tmp);
  return (*number);
}

