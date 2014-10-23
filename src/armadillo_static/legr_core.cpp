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

using namespace std;
using namespace arma;
namespace fs=boost::filesystem;

int main(int argc, char** argv)
{

    fs::path input_dir;
    fs::path output_dir;
    fs::path shared_dir;

	// Let's consider binary file read and write later
    if(argc >= 2){
        input_dir = argv[1];
    }
    if(argc >= 3){
        output_dir = argv[2];
    }else{
        output_dir = input_dir;
    }
    if(argc >= 4){
        shared_dir = argv[3];
    }else{
        shared_dir = "./";
    }

    mat X;
    fs::path Xname = "Xs_arma.mat";

    X.load((shared_dir/Xname).string(), raw_ascii);

    fs::path cur_dir(fs::current_path());
    cout << "Current directory : " << cur_dir <<endl;

    mat Yv;
    fs::path Yname = "Ys_arma.mat";
    cout << "Input : "+ (input_dir/Yname).string() <<endl;
    Yv.load((input_dir/Yname).string(),raw_ascii);
    // Convert into cube
    unsigned int nsubjects = Yv.n_cols;
    cube Y(3,3,nsubjects);

    imat idx_dti;
    fs::path idx_name = "idx_dti_arma.mat";
    idx_dti.load((shared_dir/idx_name).string(), raw_ascii);
    unsigned int nperms = idx_dti.n_rows;

    imat mask_job;
    fs::path maskname = "mask_job_arma.mat";
    mask_job.load((input_dir/maskname).string(),raw_ascii);
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
    if(output_dir.string().length() != 0){
        if (!fs::exists(output_dir)){
            cout<< "Not exists " << output_dir << endl;
            fs::create_directories(output_dir);
        }
    }
    fs::path resname = "result.mat";
    ErrMx.save((output_dir/resname).string(),raw_ascii);
    cout << "Output : "+ (output_dir/resname).string()<<endl;
  return 0;
}
