// -*- C++ -*-

#ifndef BASIC_H
#define BASIC_H

#include "globals.h"

#include <cat.h>

#include <complex>

using namespace std;
using namespace cat;

//Forward declarations
class spectral;
class Basic;
class input;


////////////////////////////////
// Class Basic declaration
////////////////////////////////

//This is the base class for basic fields
//Derived classes add the specific initialisation method

class Basic
{    
  //Members
protected:
  const input  & input_obj;
  spectral & spectral_obj;
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
  Basic(const input & input_obj__,spectral & spectral_obj__);
  //Destructor
  virtual ~Basic();
  //Public methods
public:
  void load(const string & vel_fname,const string & mag_fname,const string & temp_fname);//loads basic fields
  void save(const string & vel_fname,const string & mag_fname,const string & temp_fname);//saves basic fields  
  //Private methods
protected:
  void eval_derivatives();//evaluates the derivatives of basic fields
};

#endif


