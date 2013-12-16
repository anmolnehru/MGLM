/*
 * spd_funcs.h
 *
 *  Created on: Dec 16, 2013
 *      Author: hyunwoo
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
	int i=0;
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

mat karcher_mean_spd(&mat X,mat W,niter){

}

#endif /* SPD_FUNCS_H_ */
