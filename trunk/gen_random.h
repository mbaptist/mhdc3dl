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



#ifndef GEN_RANDOM_H
#define GEN_RANDOM_H

#include "globals.h"

#include <cat.h>

#include <complex>

using namespace std;
using namespace cat;

//Forward declarations
class Spectral;
class basic_fields;
class Input;

class gen_random
  {
  private:
    //random number generator
    random_generator random;
    const Input & input_obj;
    Spectral & spectral_obj;
  public:
    //Constructors
    explicit gen_random(const Input & input_obj__,
                        Spectral & spectral__);
    explicit gen_random(const Input & input_obj__,
                        Spectral & spectral__,
                        const int & seed_);
    //Destructor
    ~gen_random();
    //Forbidden constructors
  private:
    gen_random();
    gen_random(const gen_random &);
    //Private methods
  public:
    //generate random fields
    void gen_random_field_hat
    (CSF & field_hat,
     const int & ki, const int & kf,
     const double & alpha,const double & p,
     const bool & kind,const bool & sym,
     const std::string & br_spectrum);
    //fourier coefficients for vector field
    void gen_random_field_hat
    (CVF & field_hat,
     const int & ki, const int & kf,
     const double & alpha,const double & p,
     const bool & kind,const bool & sym,
     const std::string & br_spectrum);
    //scalar field in real space
    void gen_random_field
    (RSF & field,
     const int & ki, const int & kf,
     const double & alpha,const double & p,
     const bool & kind,const bool & sym,
     const std::string & br_spectrum);
    //vector field in real space
    void gen_random_field
    (RVF & field,
     const int & ki, const int & kf,
     const double & alpha,const double & p,
     const bool & kind,const bool & sym,
     const std::string & br_spectrum);
  };

#endif


