//============================================================================
// Name        : legr_core.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <iostream>
#include <armadillo>
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

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{

	// Let's consider binary file read and write later
    mat X;
    X.load("Xs_arma.mat",raw_ascii);

    mat Yv;
    string input_dir="";
    string output_dir="";
    if(argc > 1){
        input_dir = argv[1];
    }
    if(argc == 3){
        output_dir = argv[2];
    }else{
        output_dir = input_dir;
    }
    cout << "Input : "+ input_dir+"Ys_arma.mat" <<endl;
    Yv.load(input_dir+"Ys_arma.mat",raw_ascii);
    // Convert into cube
    unsigned int nsubjects = Yv.n_cols;
    cube Y(3,3,nsubjects);

    imat idx_dti;
    idx_dti.load("idx_dti_arma.mat",raw_ascii);
    unsigned int nperms = idx_dti.n_rows;

    imat mask_job;
    mask_job.load(input_dir+"mask_job_arma.mat",raw_ascii);
    unsigned int nvoxels = mask_job.n_rows;


    mat ErrMx(nvoxels, nperms);
    ErrMx = ErrMx.zeros();
    // Extract one voxel and reshape
    unsigned int ivoxel;
    for(ivoxel = 0; ivoxel < nvoxels; ivoxel++){
    	getY(Y,Yv,ivoxel);
    	//cout<<Y<<endl;
    	GR_legr_spd_perm(ErrMx, X, Y, idx_dti,ivoxel);
    }

    ErrMx.save(output_dir+"result.mat",raw_ascii);
    cout << "Output : "+ output_dir+"Ys_arma.mat" <<endl;
  return 0;
}
