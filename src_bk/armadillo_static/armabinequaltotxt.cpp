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
#define SZBUF 300
using namespace std;
using namespace arma;
namespace fs=boost::filesystem;

int main(int argc, char** argv)
{
    int seed = time(NULL);
    srand(seed);
	// Let's consider binary file read and write later
    if(argc < 3){
        printf("Usage: >> armabinequaltotxt filename.bin filename.txt");
        return 0;
    }
    char fbinname[SZBUF+1];
    char ftxtname[SZBUF+1];
    if(argc >= 3){
        strncpy(fbinname, argv[1], SZBUF);
    }
    if(argc >= 3){
        strncpy(ftxtname, argv[2], SZBUF);
    }

    mat A, B;
    A.load(fbinname,arma_binary);
    B.load(ftxtname,raw_ascii);
    mat C;
    C = A-B;
    double error = sum(sum(abs(C)));
    if(error == 0.0){
      printf("No error.\n");
    }
    else{
      printf("Error is %lf.\n",error);
      A.print("A:");
      B.print("B:");
    }

    return 0;
}
