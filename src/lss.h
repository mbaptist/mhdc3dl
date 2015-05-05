// -*- C++ -*-
/*
 
Copyright 2004,2005,2006 Manuel Baptista
 
This file is part of MHDC3DL
 
MHDC3DL is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
 
MHDC3DL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
*/



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
    //Reference to Input object
    Input & input_obj;

    int & n1;
    int & n2;
    int & n3;

    double & l1;
    double & l2;
    double & l3;

    Spectral spectral_obj;
    Basic basic;
    ANought ANought_obj;
    ANoughtAdjoint ANoughtAdjoint_obj;
    AOne AOne_obj;
    precond precond_obj;
    precond precond_adjoint_obj;

    //averages
    cat::Tvector<double,3> av_b_k_gamma_ij_vel[4][2][2];
    cat::Tvector<double,3> av_b_k_gamma_ij_mag[4][2][2];

  public:
    //Ctor
    lss(Input & input_obj_);
    //Dtor
    ~lss();
  private:
    lss();
    lss(const lss &);

  public:
    //Run
    void run(double & theta_min,std::complex<double> & lambda_min,double & theta_max,std::complex<double> & lambda_max);
    //void run(double & theta_min,double & lambda_min,double & theta_max,double & lambda_max);

  private:

    void solve_zero(CBVF * s,CSF * s_p,const CBVF & rhs_constant_zero,const int & i);
    void solve_one(CBVF & gamma,CBVF * s,CSF * s_p,const int & i,const int & j);
    cat::Array<double,2> eval_ep(const cat::Tvector<double,2> & q);
    cat::Array<double,2> eval_e(const cat::Tvector<double,2> & q);
    //void diag(std::complex<double> & lambda1,std::complex<double> & lambda2,const double & theta);
		void diag(std::complex<double> & lambda1,
									 double & alpha1_1,double & alpha1_2,
					std::complex<double> & lambda2,
		 double & alpha2_1,double & alpha2_2,
	 const double & theta);
    void rawsave_aux_field_hat	(const string & fname,const CBVF & field);
    void vtksave_aux_field_real	(const string & fname,const CBVF & field);
    void rawload_aux_field_hat(const string & filename,CBVF & field,const int & lr_n1,const int & lr_n2,const int & lr_n3);

  };

#endif
