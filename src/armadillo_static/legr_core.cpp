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

    mat X; //declaring a matrix armadillo code 

    mat Yv;

    fs::path Xname = "Xs_arma.mat";

    X.load((shared_dir/Xname).string(), raw_ascii);  //load into the variable X

    fs::path cur_dir(fs::current_path());
    cout << "Current directory : " << cur_dir <<endl;


    fs::path Yname = "Ys_arma.mat";
    cout << "Input : "+ (input_dir/Yname).string() <<endl;
    Yv.load((input_dir/Yname).string(),raw_ascii); // load into vY

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


    mat ErrMx1(nvoxels, nperms);
    
    mat ErrMx2(nvoxels, nperms);


    mat ErrMx1 = ErrMx1.zeros(); //initialize
    mat ErrMx2 = ErrMx2.zeros();
   

	mat X_full=X;
	mat X_part=X.shed_cols[n_cols-1]; 


    // Extract one voxel and reshape
    unsigned int ivoxel;
    for(ivoxel = 0; ivoxel < nvoxels; ivoxel++){
    	getY(Y,Yv,ivoxel);
    	//cout<<Y<<endl;
    	GR_legr_spd_perm(ErrMx1, X_part, Y, idx_dti,ivoxel);
    	GR_legr_spd_perm(ErrMx2, X_full, Y, idx_dti,ivoxel);
    }


    mat ErrMxfinal = abs(ErrMx1-ErrMx2); //difference/improvement in the errors


//Anmol's changes

	//get the ascii format of ErrMxfinal somehow, need to know ErrMxfinal format and convert to this
	mat ErrMxfinal_ascii = ErrMxfinal; //this is the asciied version, the 0th value should be the one being compared to //TO DO

    size_t length=nperms;
    float *p_value=(float*)malloc(nvoxels*sizeof(float)); //creates a p_value vector of type float for all voxels
    size_t count=0;

    for(ivoxel = 0; ivoxel < nvoxels; ivoxel++){
        count=0;
        while(length--)
            {
                if(ErrMxfinal_ascii[ivoxel][length]>ErrMxfinal_ascii[ivoxel][0])
                count++; //find out values greater than ref value
            }

     p_value[ivoxel]= float(count/nperms)*100; //typecasting

    }   


    //Writing the p_value vector to a .txt file. Alternately armabinascii maybe used

    ofstream fout("p_value.txt"); //opening an output stream for file p_value.txt
    /*checking whether file could be opened or not. If file does not exist or don't have write permissions, file
  stream could not be opened.*/
  if(fout.is_open())
    {
    //file opened successfully so we are here
    cout << "File Opened successfully!!!. Writing data from p_value to file" << endl;

        for(int i = 0; i<nperms; i++)
        {
            fout << p_value[i]; //writing ith character of p_value in the file
        }
    cout << "p_value data successfully saved into the file p_value.txt" << endl;
    }
    else //file could not be opened
    {
        cout << "File could not be opened." << endl;
    }



/*
    if(output_dir.string().length() != 0){
        if (!fs::exists(output_dir)){
            cout<< "Not exists " << output_dir << endl;
            fs::create_directories(output_dir);
        }
    }
    fs::path resname = "result.bin";
    ErrMx.save((output_dir/resname).string(),arma_binary);
    cout << "Output : "+ (output_dir/resname).string()<<endl;

    */
 // return p_value;

    return 0;
}
