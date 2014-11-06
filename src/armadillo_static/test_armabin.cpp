//============================================================================
// Name        : legr_core.cpp
// Version     : 0.1
// Copyright   : Your copyright notice
// Description : HDF5 test.
//
//   $ Hyunwoo J. Kim $  $ 2014/11/05 19:18:43 (CST) $
//============================================================================
#include <iostream>
#include <armadillo>
#include <boost/filesystem.hpp> 
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "spd_funcs.h"
#include "gr_spd.h"

using namespace std;
using namespace arma;
namespace fs=boost::filesystem;

int main(int argc, char** argv)
{
  printf("Hello World.\n");
  mat A = randu<mat>(10,5);
  A.print();
  A.save("A.bin", arma_binary);

  return 0;
}
