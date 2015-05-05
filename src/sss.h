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


#ifndef MHDCL_SSS_H
#define MHDCL_SSS_H

#include "input.h"

#include <cat.h>

#include "spectral.h"

#include "linops.h"
#include "basic.h"

#include <iostream>

#include "vzdeigen/vzdeigen.h"


class sss
  {
  private:
    //Members

    //Reference to Input object
    Input & input_obj;

    Spectral spectral_obj;
    Basic basic;
    ANought ANought_obj;

  public:
    //Ctor
    sss(Input & input_obj_);
    //Dtor
    ~sss();
  private:
    sss();
    sss(const sss &);

    //Public Methods
  public:
    //Run
    void run(double & xp,double & eim);
    void external_prodx(int * nit,double * vi, double * vo, int * m);

    //Private Methods
  private:
    void gen_random_v(double * v,const int & sym,const int & seed);
    void bv2v(double * v,const CBVF & bv);
    void v2bv(CBVF & bv,const double * v);
  };

//This class is a singleton containing a reference to the linear operator
//and a function to apply the operator to the vector. It is "instatiated"
//each time the function external_prodx_ is called
class ExternalProdx
  {
  private:
    //Members
    sss * sss_ptr;
  private:
    //Ctors
    ExternalProdx()
    {}
    ;
  public:
    //Dtor
    ~ExternalProdx()
    {}
    ;
  private:
    //Forbidden Constructors
    ExternalProdx(const ExternalProdx &);
    ExternalProdx operator=(const ExternalProdx &);
  public:
    //Instantiate
    static ExternalProdx & instance()
    {
      static ExternalProdx the_ExternalProdx;
      return the_ExternalProdx;
    }
    void switch_ptr(sss * sss_ptr__)
    {
      sss_ptr=sss_ptr__;
    };
    void operator()(int * nit,double * vi, double * vo, int * m)
    {
      sss_ptr->external_prodx(nit,vi,vo,m);
    }
  };

//External function for performing the product of the linear operator by the vector
//This function is called by the function PRODX in vzdeigen
//NOTE: In the original version of vzdeigen.f, Vlad implemented PRODX directly in fortran.
//For this version, if your linear operator is implemented in fortran, you must write the function
//external_prodx(integer nit,real(8) vi, real(8) vo, integer m) that invokes your operator.
extern"C"
  {
    void external_prodx_(int * nit,double * vi, double * vo, int * m);
  }


#endif
