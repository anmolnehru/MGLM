//============================================================================
// Name        : legr_core.cpp
// Author      : Hyunwoo J. Kim
// Version     : 0.1
// Date        : 07.16.2014
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <iostream>
#include <fstream>
#include <armadillo>
#include <boost/filesystem.hpp> 
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include "spd_funcs.h"
#include "gr_spd.h"
#define SZBUF 300
using namespace std;
using namespace arma;
namespace fs=boost::filesystem;
typedef std::numeric_limits< double > dbl;
int main(int argc, char** argv)
{
    fs::path input_dir;
    fs::path output_dir;
    fs::path shared_dir;

	// Let's consider binary file read and write later
    if(argc < 2){
        printf("Usage: >> armabin2ascii filename.bin");
        printf("Usage: >> armabin2ascii filename.bin filename.txt");
        return 0;
    }
    char finname[SZBUF+1];
    char foutname[SZBUF+1];
    if(argc >= 2){
        strncpy(finname, argv[1], SZBUF);
    }
    if(argc >= 3){
        strncpy(foutname, argv[2], SZBUF);
    }

    mat M;
//    printf("Inputfile: %s\n",finname);
    M.load(finname, arma_binary);
//    M.print();
    int i=0,j=0;
    
    if(argc == 2){ 
        cout.precision(dbl::digits10+2);
        for(i=0 ; i < M.n_rows;i++){
            for(j=0;j <M.n_cols-1;j++){
                cout << scientific << M(i,j) << " ";
            }
            cout << scientific << M(i,j) <<endl;
        }
    }
    if(argc == 3){
        ofstream myfile;
        myfile.open(foutname);
        myfile.precision(dbl::digits10+2);
        for(i=0 ; i < M.n_rows;i++){
            for(j=0;j <M.n_cols-1;j++){
                myfile << scientific << M(i,j) << " ";
            }
            myfile << scientific << M(i,j) <<endl;
        }
        myfile.close();
    }
                
    return 1;
}
