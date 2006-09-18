// -*- C++ -*-

#ifndef LINOPS_H
#define LINOPS_H

#include "input.h"

#include "block_vector.h"

#include "globals.h"
#include "block_vector.h"

#include "basic.h"

#include <cat.h>

#include <complex>


class linops_base;


//Class linops_tmp
//Contains temporaries to evaluate non-linearities
//Implemented as a singleton
class linops_tmp
{
  friend class linops_base;
  //Members
private:
  int n1;
  int n2;
  int n3;
protected:
  //Temporary Variables for nonlinear calculations
  //Fields (perturbation fields) in real space
  //velocity field in real space
  RVF vel;
  //magnetic field in real space
  RVF mag;
  //temperature field in real space
  RSF temp;
  //To evaluate derivatives of fields for nonlinear terms
  RVF curl_vel;
  RVF curl_mag;
  RVF grad_temp;
  CVF curl_vel_hat;
  CVF curl_mag_hat;
  CVF grad_temp_hat;
  //To evaluate nonlinear terms
  RSF s_nonlinear;
  RVF nonlinear;
  CSF s_nonlinear_hat;
  CVF nonlinear_hat;
private:
  //Constructor
  linops_tmp(input & input_obj__);
public:
  //Destructor
  ~linops_tmp();
private:
  //Forbidden Constructors
  linops_tmp();
  linops_tmp(const linops_tmp &);
  linops_tmp operator=(const linops_tmp &);
public:
  //Instantiate
  static linops_tmp & instance(input & input_obj__)
  {
    static linops_tmp the_linops_tmp(input_obj__);
    return the_linops_tmp;
  }
};


class linops_base
{
  //Members
protected:

  //input parameters
  int n1,n2,n3;//Spatial grid
  Real l1,l2,l3;//Spatial extent

  //physical Parameters
  Real visc,omegaz,compress,g,diff,deltat,econd,tcond;

  //reference to Spectral object
  Spectral & spectral_obj;

  //reference to basic fields
  Basic & basic;

  //Symmetry subspace used
  int sym_sub_;

  //References to Temporary Variables for nonlinear calculations
  //Fields (perturbation fields) in real space
  //velocity field in real space
  RVF & vel;
  //magnetic field in real space
  RVF & mag;
  //temperature field in real space
  RSF & temp;
  //To evaluate derivatives of fields for nonlinear terms
  RVF & curl_vel;
  RVF & curl_mag;
  RVF & grad_temp;
  CVF & curl_vel_hat;
  CVF & curl_mag_hat;
  CVF & grad_temp_hat;
  //To evaluate nonlinear terms
  RSF & s_nonlinear;
  RVF & nonlinear;
  CSF & s_nonlinear_hat;
  CVF & nonlinear_hat;

  //Constructor/destructor
public:
  //Constructor
  linops_base(input & input_obj__,
	      Spectral & spectral_obj__,
	      Basic & basic__);
  //Destructor
  ~linops_base();

public:

  //Symmetry subspace accessor
  int & sym_sub();
  const int & sym_sub() const;

  //Scalar Product
  Real scalar_prod(const CBVF & xx,
		     const CBVF & yy) const;

};



class a_nought: public linops_base
{
private:
  using linops_base::n1;
  using linops_base::n2;
  using linops_base::n3;
  using linops_base::l1;
  using linops_base::l2;
  using linops_base::l3;
  using linops_base::visc;
  using linops_base::omegaz;
  using linops_base::compress;
  using linops_base::g;
  using linops_base::diff;
  using linops_base::deltat;
  using linops_base::econd;
  using linops_base::tcond;
  using linops_base::spectral_obj;
  using linops_base::basic;
  using linops_base::sym_sub_;
  
  using linops_base::vel;
  using linops_base::mag;
  using linops_base::temp;
  using linops_base::curl_vel;
  using linops_base::curl_mag;
  using linops_base::grad_temp;
  using linops_base::curl_vel_hat;
  using linops_base::curl_mag_hat;
  using linops_base::grad_temp_hat;
  using linops_base::s_nonlinear;
  using linops_base::nonlinear;
  using linops_base::s_nonlinear_hat;
  using linops_base::nonlinear_hat;

public:

  using linops_base::sym_sub;
  using linops_base::scalar_prod;

public:
  //Constructor
  a_nought(input & input_obj__,
	   Spectral & spectral_obj__,
	   Basic & basic__);
  //Destructor
  ~a_nought();

private:
  //Forbidden Constructors
  a_nought();
  a_nought(const a_nought &);
  
  //Public Methods
public:
  //Apply operator
  CBVF operator()
    (const CBVF& w); 
  //const;

};


class a_nought_adjoint: public linops_base
{
private:
  using linops_base::n1;
  using linops_base::n2;
  using linops_base::n3;
  using linops_base::l1;
  using linops_base::l2;
  using linops_base::l3;
  using linops_base::visc;
  using linops_base::omegaz;
  using linops_base::compress;
  using linops_base::g;
  using linops_base::diff;
  using linops_base::deltat;
  using linops_base::econd;
  using linops_base::tcond;
  using linops_base::spectral_obj;
  using linops_base::basic;
  using linops_base::sym_sub_;
  
  using linops_base::vel;
  using linops_base::mag;
  using linops_base::temp;
  using linops_base::curl_vel;
  using linops_base::curl_mag;
  using linops_base::grad_temp;
  using linops_base::curl_vel_hat;
  using linops_base::curl_mag_hat;
  using linops_base::grad_temp_hat;
  using linops_base::s_nonlinear;
  using linops_base::nonlinear;
  using linops_base::s_nonlinear_hat;
  using linops_base::nonlinear_hat;

public:

  using linops_base::sym_sub;
  using linops_base::scalar_prod;

public:
  //Constructor
  a_nought_adjoint(input & input_obj__,
	   Spectral & spectral_obj__,
	   Basic & basic__);
  //Destructor
  ~a_nought_adjoint();

private:
  //Forbidden Constructors
  a_nought_adjoint();
  a_nought_adjoint(const a_nought &);
  
  //Public Methods
public:
  //Apply operator
  CBVF operator()
    (const CBVF& w);
  //const;

};


class a_one: public linops_base
{
private:
  using linops_base::n1;
  using linops_base::n2;
  using linops_base::n3;
  using linops_base::l1;
  using linops_base::l2;
  using linops_base::l3;
  using linops_base::visc;
  using linops_base::omegaz;
  using linops_base::compress;
  using linops_base::g;
  using linops_base::diff;
  using linops_base::deltat;
  using linops_base::econd;
  using linops_base::tcond;
  using linops_base::spectral_obj;
  using linops_base::basic;
  using linops_base::sym_sub_;
  
  using linops_base::vel;
  using linops_base::mag;
  using linops_base::temp;
  using linops_base::curl_vel;
  using linops_base::curl_mag;
  using linops_base::grad_temp;
  using linops_base::curl_vel_hat;
  using linops_base::curl_mag_hat;
  using linops_base::grad_temp_hat;
  using linops_base::s_nonlinear;
  using linops_base::nonlinear;
  using linops_base::s_nonlinear_hat;
  using linops_base::nonlinear_hat;

public:

  using linops_base::sym_sub;
  using linops_base::scalar_prod;

public:
  //Constructor
  a_one(input & input_obj__,
	   Spectral & spectral_obj__,
	   Basic & basic__);
  //Destructor
  ~a_one();

private:
  //Forbidden Constructors
  a_one();
  a_one(const a_nought &);
  
  //Public Methods
public:
  //Apply operator
  CBVF operator()
    (const CBVF& w,
     const int & index); //const;

};




class precond
{
  //Members
private:
  Spectral & spectral_obj;
  Real qq_;
public:
  //Constructor
	precond(Spectral & spectral_obj__,float qq__);
  //Destructor
  ~precond();
private:
	//Forbidden ctors	
		precond();
  //Pulic Methods
public:
  CBVF operator()
    (const CBVF& x) const;
};










//////////////////////////////////////////////////////////////////////

#if 0



class linops
{
  //Members
private:

  //input parameters
  int n1,n2,n3;//Spatial grid
  Real l1,l2,l3;//Spatial extent
  //physical Parameters
  Real visc,omegaz,compress,g,diff,deltat,econd,tcond;

  //reference to Spectral object
  Spectral & spectral_obj;

  //reference to basic fields
  Basic & basic;

  //Temporary Variables for nonlinear calculations
  //Fields (perturbation fields) in real space
  //velocity field in real space
  RVF vel;
  //magnetic field in real space
  RVF mag;
  //temperature field in real space
  RSF temp;
  //To evaluate derivatives of fields for nonlinear terms
  RVF curl_vel;
  RVF curl_mag;
  RVF grad_temp;
  CVF curl_vel_hat;
  CVF curl_mag_hat;
  CVF grad_temp_hat;
  //To evaluate nonlinear terms
  RSF s_nonlinear;
  RVF nonlinear;
  CSF s_nonlinear_hat;
  CVF nonlinear_hat;
    
  //BlockVector to hold the gradient part of a_nought
  CBVF gradient_;

  //Symmetry subspace used
  int sym_sub_;


  //Constructor/destructor
public:
  //Constructor
  linops(input & input_obj__,Spectral & spectral_obj__,Basic & basic__);
  //Destructor
  ~linops();


  //Public methods
public:
  //Linear Operators
  //a_nought
  CBVF eval_a_nought(const CBVF& w);
  //a_nought_star
  CBVF eval_a_nought_star(const CBVF & w);
  //b
  CBVF eval_b(const CBVF & field,const int & index);

  //Accessor to the gradient part of a_nought
  CBVF & gradient();//accessor to BlockVector that holds the gradient part of a_nought

  //Test the adjoint by its definition
  void test_adjoint(CBVF & a,
		    CBVF & b);

  //Symmetry subspace accessor
  int & sym_sub();
  const int & sym_sub() const;

  //preconditioner for cgsolver
  CBVF 
  precond(const CBVF& x,
	  const Real & qq) const;


};


class a_nought
{
public:
  a_nought(linops & linops_obj_,Spectral & spectral_obj_):linops_obj(linops_obj_),spectral_obj(spectral_obj_){};
  ~a_nought(){};
  CBVF apply_direct
  (const CBVF & xx) const;
  CBVF apply_adjoint
  (const CBVF & xx) const; 
  Real scalar_prod(const CBVF & xx,
		     const CBVF & yy) const;
  CBVF precond(const CBVF & xx,
					 Real & qq) const;
private:
  linops & linops_obj;
  Spectral & spectral_obj;
};


class a_nought_star
{
public:
  a_nought_star(linops & linops_obj_,Spectral & spectral_obj_):linops_obj(linops_obj_),spectral_obj(spectral_obj_){};
  ~a_nought_star(){};
  CBVF apply_direct
  (const CBVF & xx) const;
  CBVF apply_adjoint
  (const CBVF & xx) const; 
  Real scalar_prod(const CBVF & xx,
		     const CBVF & yy) const;
  CBVF precond(const CBVF & xx,
					 Real & qq) const;
private:
  linops & linops_obj;
  Spectral & spectral_obj;
};

#endif











#endif
