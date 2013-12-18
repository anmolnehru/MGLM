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
	X.load("X_arma.mat",raw_ascii);

    mat Yv;
    Yv.load("Ys_arma.mat",raw_ascii);
    // Convert into cube
    unsigned int nsubjects = Yv.n_cols;
    cube Y(3,3,nsubjects);

    imat idx_dti;
    idx_dti.load("idx_dti_test_int_arma.mat",raw_ascii);
    unsigned int nperms = idx_dti.n_rows;

    imat mask_job;
    mask_job.load("mask_job_arma.mat",raw_ascii);
    unsigned int nvoxels = mask_job.n_rows;
    nvoxels=1;

    mat ErrMx(nvoxels, nperms);
    ErrMx = ErrMx.zeros();
    // Extract one voxel and reshape
    unsigned int ivoxel;
    for(ivoxel = 0; ivoxel < nvoxels; ivoxel++){
    	getY(Y,Yv,ivoxel);
    	cout<<Y<<endl;
    	GR_legr_spd_perm(ErrMx, X, Y, idx_dti,ivoxel);
    }

    cout << ErrMx <<endl;






  return 0;
}
