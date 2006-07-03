// -*- C++ -*-



#ifndef SPECTRAL_H
#define SPECTRAL_H

#include "globals.h"

#include <goops.h>

#include <cat.h>


class spectral : public SpectralFourierLayer
{

public:
	using SpectralFourierLayer::wv;
	using SpectralFourierLayer::wv2;	
	using SpectralFourierLayer::wnmax;
	using SpectralFourierLayer::nwn;
	using SpectralFourierLayer::wnstep;
	
 public:
  using SpectralFourierLayer::dealias;
  using SpectralFourierLayer::poisson_hat;
  using SpectralFourierLayer::d_dx_index_hat;
  using SpectralFourierLayer::grad_hat;
  using SpectralFourierLayer::div_hat;
  using SpectralFourierLayer::curl_hat;
  using SpectralFourierLayer::lap_hat;
  using SpectralFourierLayer::remove_gradient;
  using SpectralFourierLayer::scalar_prod;
	using SpectralFourierLayer::eval_energ_spec;
	using SpectralFourierLayer::pnvh_hat;
  using SpectralFourierLayer::pnvh;

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
  
  void pnvh_hat(const CBVF & field);

};

#endif
