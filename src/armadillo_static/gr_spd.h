//Anmol, 12/11

#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;

#include "spd_funcs.h"

/*
 *  Explicit allocation need here.
 *  X is centered.
 */


void GR_legr_spd_perm(mat& ErrMx, const mat& X, const cube& Y, const imat & idx_dti, unsigned ithvox){

	unsigned int ndata = X.n_cols;
	///cout<<ndata<<endl;
	unsigned int nperms = idx_dti.n_rows;
	unsigned int niter = 100;



	mat p(3,3);
	karcher_mean_spd(p, Y, niter);
	mat sqrtp = zeros<mat>(3,3);
	mat invg = zeros<mat>(3,3);
	mat g = zeros<mat>(3,3);
	get_g_invg(g,invg,p);
	sqrtm(sqrtp, p);

	cube logY(3,3,ndata);
	logmap_pt2array_spd(logY, p, Y);

	mat Yv(6,ndata);
	cube S(3,3,ndata);

	embeddingR6_vecs(Yv, S, p, logY);
	cube logYvhat_perm(6,ndata,nperms);

	unsigned int iperm;
	mat Xc(X.n_rows, X.n_cols);
	mat L;
	mat logYv_hat(6,ndata);
	cube V_hat(3,3,ndata);


	for(iperm=0; iperm < nperms; iperm++){
		mxpermute(Xc, X, idx_dti, iperm);
		//    A \ B	  	solve(A,B)
		L = solve(Xc.t(),Yv.t());
		if(iperm == 0){
		}
		logYv_hat = L.t()*Xc;
		logYvhat_perm.slice(iperm) = logYv_hat;
	}

	// Check distance
    unsigned int idata;
    for(idata =0; idata <ndata;idata++){
    	dist_M_pt2array(ErrMx,p, sqrtp, g, invg, Y, logYvhat_perm, idata, ithvox);
    }
}
