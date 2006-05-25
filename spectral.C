

#include "globals.h"

#include "spectral.h"

#include <goops.h>

#include <cat.h>

using namespace cat;
//using namespace goops;

spectral::spectral(const int & n1__,const int & n2__,const int & n3__,
		   const Real & l1__,const Real & l2__,const Real & l3__):
  spectral_fourier_layer(n1__,n2__,n3__,
			 l1__,l2__,l3__)
{
}

spectral::~spectral()
{
}

//block_vectors
void spectral::dealias(CBVF & field) const
{
  spectral_fourier_layer::dealias(field.vel());
  spectral_fourier_layer::dealias(field.mag());
  spectral_fourier_layer::dealias(field.temp());
}


//Laplacian of block_vector
CBVF spectral::lap_hat(const CBVF & field)
{
  CBVF out(field);
  out.vel()=spectral_fourier_layer::lap_hat(field.vel());
  out.mag()=spectral_fourier_layer::lap_hat(field.mag());
  out.temp()=spectral_fourier_layer::lap_hat(field.temp());
  return out;
}

//block_vector
CBVF spectral::remove_gradient(CBVF & bfield,const bool kind)
{
  CBVF out(bfield);
  //Extract gradient part from velocity
  out.vel()=spectral_fourier_layer::remove_gradient(bfield.vel(),kind);
  //Extract gradient part from magnetic field
  out.mag()=spectral_fourier_layer::remove_gradient(bfield.mag(),kind);
  //Temperature has no gradient part!
  out.temp()=0;
  return out;
}

//block vectors kind ccs/ccs/s
Real spectral::scalar_prod(const CBVF & x,const CBVF & y) const
{
  return Real(spectral_fourier_layer::scalar_prod(x.vel(),y.vel())+
	      spectral_fourier_layer::scalar_prod(x.mag(),y.mag())+
	      spectral_fourier_layer::scalar_prod(x.temp(),y.temp()));
}

//block vector
void spectral::pnvh(const CBVF & field)
{
  cout << "vel" << endl;
  spectral_fourier_layer::pnvh(field.vel());
  cout << "mag" << endl;
  spectral_fourier_layer::pnvh(field.mag());
  cout << "temp" << endl;
  spectral_fourier_layer::pnvh(field.temp());
}

