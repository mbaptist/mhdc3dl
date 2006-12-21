// -*- C++ -*-

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
class input;

class gen_random
{
private:
  //random number generator
	random_generator random;
	const input & input_obj;
	Spectral & spectral_obj;
public:
  //Constructors
	explicit gen_random(const input & input_obj__,
	                    Spectral & spectral__);
	explicit gen_random(const input & input_obj__,
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


