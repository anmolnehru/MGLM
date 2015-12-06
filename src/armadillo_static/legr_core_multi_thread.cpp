//============================================================================
// Name        : legr_core.cpp
// Author      : Hyunwoo J. Kim
//	       : Anmol Mohanty
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
#include <omp.h>

//These would be the main header files including all the 'meat'. Scramble these
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

using namespace arma; //for the armadillo library, this is separate from boost

namespace fs=boost::filesystem; //namespace is being declared

int main(int argc, char** argv)
{

    //creating 3 paths named below
    fs::path input_dir;
    fs::path output_dir;
    fs::path shared_dir;

	// Let's consider binary file read and write later ?? ------------ What is this ask Hyunwoo


//interesting-figure out why this is being done?
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

    mat X; //declaring a matrix armadillo code || Armadillo API 

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

    //Could this be causing the issue?
    cube Y(3,3,nsubjects);

    imat idx_dti;
    fs::path idx_name = "idx_perm_arma.mat";
    idx_dti.load((shared_dir/idx_name).string(), raw_ascii);
    float nperms = idx_dti.n_rows;

    imat mask_job;
    fs::path maskname = "mask_job_arma.mat";
    mask_job.load((input_dir/maskname).string(),raw_ascii);
    unsigned int nvoxels = mask_job.n_rows;

    mat ErrMx1(nvoxels, nperms);
    mat ErrMx2(nvoxels, nperms);

    ErrMx1 = ErrMx1.zeros(); //initialize
    ErrMx2 = ErrMx2.zeros();
   
    unsigned int num_cols = X.n_cols-1;
    unsigned int num_rows = X.n_rows-1;

	//mat X_full=X;
	//mat X_part=X.shed_cols(num_cols);   //From API of arma
//cout<<X.row(0)<<endl;
    // Extract one voxel and reshape
    unsigned int ivoxel;

#pragma omp parallel for
    for(ivoxel = 0; ivoxel < nvoxels; ivoxel++){
    	getY(Y,Yv,ivoxel);


//----------------- Anmol Initial Modifications------------------//
//cout<<X.cols(0,num_cols-1)<<endl; //some issue there
//cout<<X<<endl;
//	//system("PAUSE"); -- Windows only

	//getchar();
//----------------- Anmol Initial Modifications------------------//
    	//cout<<Y<<endl;
    //	GR_legr_spd_perm(ErrMx1, X_part, Y, idx_dti,ivoxel);
    //	GR_legr_spd_perm(ErrMx2, X_full, Y, idx_dti,ivoxel);
//----------------- Anmol Initial Modifications------------------//

//simpler method //just pass truncated X when doing the func call

//testing beginning -- Anmol  11/10




	
        GR_legr_spd_perm(ErrMx1, X.row(0), Y, idx_dti,ivoxel);   //generalize this one

        GR_legr_spd_perm(ErrMx2, X, Y, idx_dti,ivoxel);
	
    }

    mat ErrMxfinal = ErrMx1-ErrMx2; //difference/improvement in the errors


//Anmol's changes

	//get the ascii format of ErrMxfinal somehow, need to know ErrMxfinal format and convert to this
	//mat ErrMxfinal = ErrMxfinal; //this is the asciied version, the 0th value should be the one being compared to //TO DO

    size_t length=0;//initializing length as a counter
    float *p_value=(float*)malloc(nvoxels*sizeof(float)); //creates a p_value vector of type float for all voxels
    float count=0;

//	cout<<endl<<"ErrMxfina"<<endl;

//    for(ivoxel = 0; ivoxel < nvoxels; ivoxel++){
//    	length=nperms;
//    	while(length--)
//    		cout<<ErrMxfinal(ivoxel,length)<<" ";
//	    	cout<<endl;
//	}


cout<<endl<<endl<<"PVALUES"<<endl;

//should automatically spawn the right amount of threads
////disabling the below in hope that the error will go away
#pragma omp parallel for
    for(ivoxel = 0; ivoxel < nvoxels; ivoxel++){
        count=0;
    	length=nperms;
        while(length--)
            {
                if(ErrMxfinal(ivoxel,length)>ErrMxfinal(ivoxel,0))
	             count++; //find out values greater than ref value
            }
//		cout<<"count-"<<count<<" nperms-"<<nperms;
    	p_value[ivoxel]= count/nperms; //typecasting
    

    	//introduce logic for early termination et al
//    	cout<<":: p_value="<<p_value[ivoxel]<<endl;
    }   


    //Writing the p_value vector to a .txt file. Alternately armabinascii maybe used

	//Checking if file exists and deletes it
	
	if( remove( "p_value.txt" ) == 0 ) //file exists
	    perror( "File existed and has been cleaned up" );  //recheck if this actually removes the file?
	else puts("Writing p_values to file");

    ofstream fout("p_value.txt"); //opening an output stream for file p_value.txt
    /*checking whether file could be opened or not. If file does not exist or don't have write permissions, file
  stream could not be opened.*/
  if(fout.is_open())
    {
    //file opened successfully so we are here
    cout << "File Opened successfully!!!. Writing data from p_value to file p_value.txt" << endl;
		//this line below is incorrect, p_values are computed for a voxel and not over permutations
    //    for(int i = 0; i<nperms; i++)
        //this update should make it correct, however note that this cannot be parallelized, as values need be printed in seq
	for(int i=0;i<nvoxels;i++)
        {
            fout << p_value[i]; //writing ith character of p_value in the file
            fout <<",\n";
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
