#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;

#include "spd_funcs.h"

/*
 *  Explicit allocation need here.
 *  X is centered.
 */
void GR_legr_spd_perm(vec& ErrMx, const mat& X, const cube& Y, const mat & idx_dti){

	unsigned int ndata = X.n_cols;
	unsigned int nperms = idx_dti.n_rows;
	unsigned int niter = 100;
	unsigned int ndimx = X.n_rows;

	mat p(3,3);
	karcher_mean_spd(p, Y, niter);

	cube logY(3,3,ndata);
	logmap_pt2array_spd(logY, p, Y);

	mat Yv(6,ndata);
	cube S(3,3,ndata);

	embeddingR6_vecs(Yv, S, p, logY);
	cube Yhat_perm(3,3,ndata);
	unsigned int iperm;

	mat Xc(ndimx, ndata);
	mat L;
	mat logYv_hat(6,ndata);
	cube V_hat(3,3,ndata);
	for(iperm=0; iperm < nperms; iperm++){
		mxpermute(Xc,X, idx_dti.row(iperm));
		//    A \ B	  	solve(A,B)
		L = solve(Yv.t(),Xc.t()).t();
		logYv_hat = L*Xc;


	}





}
