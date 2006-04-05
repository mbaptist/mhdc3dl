// -*- C++ -*-



#ifndef SPECTRAL_H
#define SPECTRAL_H

#include "globals.h"

#include <goops.h>

#include <cat.h>


class spectral : public spectral_fourier_layer
{
 public:
  using spectral_fourier_layer::dealias;
  using spectral_fourier_layer::poisson_hat; 
  using spectral_fourier_layer::d_dhorizontal_hat; 
  using spectral_fourier_layer::grad_hat; 
  using spectral_fourier_layer::div_hat;
  using spectral_fourier_layer::curl_hat;
  using spectral_fourier_layer::lap_hat;
  using spectral_fourier_layer::remove_gradient;
  using spectral_fourier_layer::scalar_prod;
  using spectral_fourier_layer::pnvh;

 public:

  spectral(const int & n1__,const int & n2__,const int & n3__,
	   const Real & l1__,const Real & l2__,const Real & l3__);

  ~spectral();
    
  void dealias(CBVF & field)  const;//block_vector

  ///Laplacian of block_vector
  CBVF lap_hat(const CBVF & field);
  
  //remove gradient part after of a block_vector
  CBVF remove_gradient(CBVF & bfield,const bool kind);	
  
  Real scalar_prod(const CBVF & xx,const CBVF & yy) const ;
  
  void pnvh(const CBVF & field);

};

#endif
