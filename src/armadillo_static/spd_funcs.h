/*
 * spd_funcs.h
 *
 *  Created on: Jan 08, 2014
 *      Author: Hyunwoo
 *
 */
#include <iostream>
#include <armadillo>
using namespace arma;
using namespace std;

#ifndef SPD_FUNCS_H_
#define SPD_FUNCS_H_

/*
 * get_g_invg is the common part in logarithm map and exponential map related functions.
 */

void get_g_invg(mat& g, mat& invg, const mat &p){
    vec d;
    mat U;
    eig_sym(d, U, p);
    mat D = diagmat(sqrt(d));
    g = U*D;
    int i;
    for(i=0;i<3;i++){
    	D(i,i) = 1./D(i,i);
    }
    invg = D*U.t();
}

mat expmap_spd(const mat & P, const mat& X){
	mat exp_p_x;
    if( norm(X,2) < 1e-18){
    	exp_p_x = P;
    	return exp_p_x;
    }

    mat g,invg;
    get_g_invg(g, invg, P);
    mat Y = invg*X*invg.t();
    vec s;
    mat V;
    eig_sym(s, V, Y);
    mat gv = g*V;
    exp_p_x = gv*diagmat(exp(s))*gv.t();

	return exp_p_x;
}


/*
 * V is a set of tangent vectors log_{P}X_{i}
 * X is a set of points on M.
 */
void logmap_pt2array_spd(cube& V, const mat& p,const cube& X){

	vec d;
	mat U;
	mat g,invg;
	get_g_invg(g, invg, p);
	unsigned int i=0;
	mat tmp(3,3);
	mat Y;
	mat H;
	for(i=0;i<X.n_slices;i++){
		if(norm(p-X.slice(i),2) <1e-18){
			V.slice(i) = tmp.zeros();
			continue;
		}
		Y = invg*X.slice(i)*invg.t();
		vec s;
		mat U;
		eig_sym(s,U,Y);
		H = g*U;
		V.slice(i) = H*diagmat(log(s))*H.t();
	}
}
/*
 * Arithmetic mean of matrices in cube.
 */
mat mean_cube(const cube& X){
	unsigned int i = 0;
	mat mV(X.n_rows,X.n_cols);
	mV = mV.zeros();
	for(i=0; i<X.n_slices; i++){
		mV += X.slice(i);
	}
	return mV/X.n_slices;
}

/*
 *
 */
void karcher_mean_spd(mat& xbar, const cube& X, const int niter){
	int iter =0;
	cube V(X.n_rows, X.n_cols, X.n_slices);

	xbar = X.slice(0);
	mat phi(3,3);

    for(iter=0; iter < niter; iter++){
    	logmap_pt2array_spd(V, xbar, X);
    	phi = mean_cube(V);
    	if(norm(phi,2) < 1e-18)break;
    	xbar = expmap_spd(xbar, phi);
    }
}

/*
 * S is a memory pool to avoid frequent dynamic memory allocation for speed up.
 * Need to be checked.
 */
void embeddingR6_vecs(mat& Vnew, cube& S, const mat& p, const cube& V){
	int nmx = V.n_slices;
	vec w(6);
	w[0] = 1;
	w[1] = sqrt(2); w[2] = sqrt(2);
	w[3]= 1; w[4] = sqrt(2); w[5] = 1;

	int i;
	mat U;
	vec d;
	eig_sym(d,U,p);

	mat sqrtinvp;
	d = sqrt(d);
	for(i=0;i<3;i++)d(i) = 1/d(i);
	sqrtinvp = U*diagmat(d)*U.t();

	for(i=0;i<nmx;i++){
		S.slice(i) = sqrtinvp*V.slice(i)*sqrtinvp;
		Vnew(0,i) = w(0)*S(0,0,i);
		Vnew(1,i) = w(1)*S(0,1,i);
		Vnew(2,i) = w(2)*S(0,2,i);
		Vnew(3,i) = w(3)*S(1,1,i);
		Vnew(4,i) = w(4)*S(1,2,i);
		Vnew(5,i) = w(5)*S(2,2,i);
	}
}

/*
 * Vectors 2 symmetric matrices.
 * M is a symmetric matrix.
 */
void invembeddingR6_vecs(cube& Vnew, const mat& p,const mat& V){
	int nmx = V.n_cols;
	int i;
	vec w(6);
	w[0]= 1; w[1] = sqrt(2); w[2] = sqrt(3);
	w[3]= 1; w[4] = sqrt(2); w[5] = 1;
    mat U;
    vec d;
    mat sqrtp = U*diagmat(sqrt(d))*U.t();

	for(i=0;i < nmx;i++){
		Vnew(0,0,i) = 1/w(0)*V(0,i);
		Vnew(0,1,i) = 1/w(1)*V(1,i);
		Vnew(1,0,i) = 1/w(1)*V(1,i);
		Vnew(0,2,i) = 1/w(2)*V(2,i);
		Vnew(2,0,i) = 1/w(2)*V(2,i);
		Vnew(1,1,i) = 1/w(3)*V(3,i);
		Vnew(1,2,i) = 1/w(4)*V(4,i);
		Vnew(2,1,i) = 1/w(4)*V(4,i);
		Vnew(2,2,i) = 1/w(5)*V(5,i);

		Vnew.slice(i) = sqrtp*Vnew.slice(i)*sqrtp;
	}
}
/*
 * For speedup, we assume that most of variables are allocated outside of the function.
 * It allows us to minimize the number of dynamic allocation.
 */
void dist_M_pt2array(mat& ErrMx, const mat&p, const mat& sqrtp, const mat& g, const mat& invg, const cube& Y, const cube& logYvhat_perm, unsigned int idata, unsigned int ithvox){

	vec w(6);
	w[0]= 1; w[1] = sqrt(2); w[2] = sqrt(2);
	w[3]= 1; w[4] = sqrt(2); w[5] = 1;


	mat Yvihat = zeros<mat>(3,3);
	mat Yihat = zeros<mat>(3,3);
	mat Z = zeros<mat>(3,3);
	mat V = zeros<mat>(3,3);
	vec s = zeros<vec>(3);
	mat gV = zeros<mat>(3,3);

	// idata
	mat U = zeros<mat>(3,3);
	vec d = zeros<vec>(3);
	eig_sym(d,U,Y.slice(idata));
	unsigned int i;
	d = sqrt(d);
	for(i=0;i<3;i++)d(i) = 1/d(i);
	mat invg_Yi = U*diagmat(d)*U.t();

    unsigned int nperms = logYvhat_perm.n_slices;
    unsigned iperm;
    for(iperm=0; iperm<nperms; iperm++){
    	// invembeddingR6_vecs
    	Yvihat(0,0) = 1/w(0)*logYvhat_perm(0,idata,iperm);
		Yvihat(0,1) = 1/w(1)*logYvhat_perm(1,idata,iperm);
		Yvihat(1,0) = 1/w(1)*logYvhat_perm(1,idata,iperm);
		Yvihat(0,2) = 1/w(2)*logYvhat_perm(2,idata,iperm);
		Yvihat(2,0) = 1/w(2)*logYvhat_perm(2,idata,iperm);
		Yvihat(1,1) = 1/w(3)*logYvhat_perm(3,idata,iperm);
		Yvihat(1,2) = 1/w(4)*logYvhat_perm(4,idata,iperm);
		Yvihat(2,1) = 1/w(4)*logYvhat_perm(4,idata,iperm);
		Yvihat(2,2) = 1/w(5)*logYvhat_perm(5,idata,iperm);

		// Tangent vector
		Z = sqrtp*Yvihat*sqrtp;
		if(norm(Z,2) <1e-18){
			Yihat = p;
		}else{
			// exponential map
			Z = invg*Z*invg.t();
			eig_sym(s,V,Z);
			gV = g*V;
			Yihat = gV*diagmat(exp(s))*gV.t();
		}
		eig_sym(s,V, invg_Yi*Yihat*invg_Yi.t());
    	ErrMx(ithvox,iperm)+= sum(square(log(s))); // Doubt about the data type of return from pow function.
    }
}
/*
 * Extract Y matrices from concatenated vector forms of data.
 */
void getY(cube&Y, const mat&Yv,unsigned int ivoxel){
	unsigned int start_idx = 6*ivoxel;
	unsigned int nsubjects = Yv.n_cols;
	unsigned int i;
	for(i=0; i < nsubjects; i++){
		Y(0,0,i) = Yv(start_idx,i);
		Y(0,1,i) = Yv(start_idx+1,i);
		Y(1,0,i) = Yv(start_idx+1,i);
		Y(0,2,i) = Yv(start_idx+2,i);
		Y(2,0,i) = Yv(start_idx+2,i);
		Y(1,1,i) = Yv(start_idx+3,i);
		Y(1,2,i) = Yv(start_idx+4,i);
		Y(2,1,i) = Yv(start_idx+4,i);
		Y(2,2,i) = Yv(start_idx+5,i);
	}
}

/*
 * Index is a row vector. The type is a submatrix of a matrix
 */
void mxpermute(mat& Xperm,const mat& X, const imat& idx, unsigned int iperm){


	unsigned int i;
	for(i=0; i < idx.n_cols; i++){
		Xperm.col(i) = X.col(idx(iperm,i)-1);
	}
}

void sqrtm(mat& outmx,const mat& inmx){
	mat U;
	vec d;
	eig_sym(d,U,inmx);
	outmx = U*diagmat(sqrt(d))*U.t();
}

#endif /* SPD_FUNCS_H_ */

