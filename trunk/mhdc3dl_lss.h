// -*- C++ -*-

#ifndef MHDCL_LSS_H
#define MHDCL_LSS_H

#include "input.h"

#include "basic_fields.h"
#include "linops.h"
#include "block_vector.h"

#include <lass.h>

#include <iostream>
#include <string>

#include <cat.h>
#include <goops.h>

using namespace std;
//using namespace cg;
using namespace cat;

//Pointers to global objects

input * input_obj;

spectral * spectral_obj;

basic_fields * basic;

a_nought * a_nought_obj;
a_nought_adjoint * a_nought_adjoint_obj;
a_one * a_one_obj;

precond * precond_obj;
precond * precond_adjoint_obj;

//averages
cat::tvector<double,3> av_b_k_gamma_ij_vel[4][2][2];
cat::tvector<double,3> av_b_k_gamma_ij_mag[4][2][2];

void save_block_vector(CBVF & field,const string & fname);

void solve_zero(CBVF * s,
		CSF * s_p,
		const CBVF & rhs_constant_zero,
		const int & i);
void solve_one(CBVF & g,
	       CBVF * s,
	       CSF * s_p,
	       const int & i,const int & j);

cat::array<double,2> eval_ep(const cat::tvector<double,2> & q);
cat::array<double,2> eval_e(const cat::tvector<double,2> & q);

void diag(double & lambda1,double & lambda2,const cat::array<double,2> & matrix);


#endif
