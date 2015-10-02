//============================================================================
// Name        : legr_core.cpp
// Author      : Hyunwoo J. Kim
// Version     : 0.1
// Date        : 07.16.2014
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <iostream>
#include <armadillo>
#include <boost/filesystem.hpp> 
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "spd_funcs.h"
#include "gr_spd.h"

#define DIM_DTI 6
// Index of elements in save files
#define iDxx 0
#define iDxy 1
#define iDxz 2
#define iDyy 3
#define iDyz 4
#define iDzz 5
#define SZBUF 1000
using namespace std;
using namespace arma;
namespace fs=boost::filesystem;

int main(int argc, char** argv)
{
    fs::path input_dir;
    fs::path output_dir;
    fs::path shared_dir;

	// Let's consider binary file read and write later
    if(argc == 1){
        printf("Usage: >> show_arma_mx filename.bin");
        return 0;
    }
    char finname[SZBUF+1];
    if(argc >= 2){
        strncpy(finname, argv[1], SZBUF);
    }
    mat M;
//    printf("Inputfile: %s\n",finname);
    M.load(finname, arma_binary);
    M.print();

    return 1;
}
