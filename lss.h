// -*- C++ -*-

#ifndef MHDCL_H
#define MHDCL_H

#include "input.h"

#include "block_vector.h"

#include "spectral.h"

#include "basic.h"
#include "linops.h"

#include <cat.h>

#include <lass.h>

#include <iostream>
#include <string>

#include "spectral.h"

//LARGE SCALE STABILITY
class lss
{
private:
  //Reference to input object
  input & input_obj;
	
  Spectral spectral_obj;
  Basic basic;
  a_nought a_nought_obj;
  a_nought_adjoint a_nought_adjoint_obj;
  a_one a_one_obj;
  precond precond_obj;
  precond precond_adjoint_obj;

  //averages
  cat::tvector<double,3> av_b_k_gamma_ij_vel[4][2][2];
  cat::tvector<double,3> av_b_k_gamma_ij_mag[4][2][2];

public:
  //Ctor
  lss(input & input_obj_);
  //Dtor
  ~lss();
private:
  lss();
  lss(const lss &);


public:
  //Run
  void run(double & lambda_minimal,
	   double & lambda_maximal);
  
private:

  void save_BlockVector(CBVF & field,
			 const string & fname);
	
	void rawsave_aux_field_hat	(const string & fname,const CBVF & field);
	void vtksave_aux_field_real	(const string & fname,const CBVF & field);
	void rawload_aux_field_hat(const string & filename,CBVF & field,const int & lr_n1,const int & lr_n2,const int & lr_n3);
	
  void solve_zero(CBVF * s,
		  CSF * s_p,
		  const CBVF & rhs_constant_zero,
		  const int & i);
  void solve_one(CBVF & gamma,
		 CBVF * s,
		 CSF * s_p,
		 const int & i,const int & j);
  cat::array<double,2> eval_ep(const cat::tvector<double,2> & q);
  cat::array<double,2> eval_e(const cat::tvector<double,2> & q);
  void diag(double & lambda1,
	    double & lambda2,
	    const cat::array<double,2> & matrix);
  
};



#endif
