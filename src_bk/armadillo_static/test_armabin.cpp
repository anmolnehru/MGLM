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
#include <limits>
#include <boost/filesystem.hpp> 
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


#include "spd_funcs.h"
#include "gr_spd.h"

using namespace std;
using namespace arma;
namespace fs=boost::filesystem;

int main(int argc, char** argv)
{
  int seed = time(NULL);
  srand(seed);
  printf("Hello World.\n");
  mat A = randu<mat>(10,5);
  A.print();
  A.save("A.bin", arma_binary);
  mat B;
  B.load("A.bin",arma_binary);
  mat C;
  C = A-B;
  double error = sum(sum(abs(C)));
  if(error == 0.0){
    printf("No error.\n");
  }
  else{
    printf("Error is %lf.\n",error);
  }

  return 0;
}
