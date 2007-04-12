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



#ifndef BASIC_H
#define BASIC_H

#include "globals.h"

#include <cat.h>

#include <complex>

using namespace std;
using namespace cat;

//Forward declarations
class Spectral;
class Basic;
class Input;


////////////////////////////////
// Class Basic declaration
////////////////////////////////

//This is the base class for basic fields
//Derived classes add the specific initialisation method

class Basic
  {
    //Members
  protected:
    const Input  & input_obj;
    Spectral & spectral_obj;
    const int & n1;
    const int & n2;
    const int & n3;
    //fields
    RVF vel_;
    RVF mag_;
    RSF temp_;
    RVF curl_vel_;
    RVF curl_mag_;
    RVF curl_curl_mag_;
    RVF grad_temp_;
    //Accessors
  public:
    RVF & vel();
    RVF & mag();
    RSF & temp();
    RVF & curl_vel();
    RVF & curl_mag();
    RVF & curl_curl_mag();
    RVF & grad_temp();
    //Constructors/destructor
  public:
    //Constructors
    Basic(const Input & input_obj__,Spectral & spectral_obj__);
    //Destructor
    virtual ~Basic();
    //Public methods
  public:
    void load(const string & fielname);
    void rawload_hat(const string & filename,const int & lr_n1,const int & lr_n2,const int & lr_n3);
    void rawsave_hat(const string & filename);
    void vtksave_real(const string & filename);
    void save_energ_spec(const string & filename);
    //Private methods
  protected:
    void eval_derivatives();//evaluates the derivatives of basic fields
  };

#endif


