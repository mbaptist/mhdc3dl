

#include "globals.h"

#include "spectral.h"

#include <goops.h>

#include <cat.h>

using namespace cat;
//using namespace goops;

Spectral::Spectral(const int & n1__,const int & n2__,const int & n3__,
		   const Real & l1__,const Real & l2__,const Real & l3__):
  SpectralFourierLayer(n1__,n2__,n3__,
			 l1__,l2__,l3__)
{
}

Spectral::~Spectral()
{
}

//BlockVectors
void Spectral::dealias(CBVF & field) const
{
  SpectralFourierLayer::dealias(field.vel());
  SpectralFourierLayer::dealias(field.mag());
  SpectralFourierLayer::dealias(field.temp());
}


//Laplacian of BlockVector
CBVF Spectral::lap_hat(const CBVF & field)
{
  CBVF out(field);
  out.vel()=SpectralFourierLayer::lap_hat(field.vel());
  out.mag()=SpectralFourierLayer::lap_hat(field.mag());
  out.temp()=SpectralFourierLayer::lap_hat(field.temp());
  return out;
}

//BlockVector
CBVF Spectral::remove_gradient(CBVF & bfield,const bool kind)
{
  CBVF out(bfield);
  //Extract gradient part from velocity
  out.vel()=SpectralFourierLayer::remove_gradient(bfield.vel(),kind);
  //Extract gradient part from magnetic field
  out.mag()=SpectralFourierLayer::remove_gradient(bfield.mag(),kind);
  //Temperature has no gradient part!
  out.temp()=0;
  return out;
}

//block vectors kind ccs/ccs/s
Real Spectral::scalar_prod(const CBVF & x,const CBVF & y,const bool & kind) const
{
  return Real(SpectralFourierLayer::scalar_prod(x.vel(),y.vel(),kind)+
	      SpectralFourierLayer::scalar_prod(x.mag(),y.mag(),kind)+
	      SpectralFourierLayer::scalar_prod(x.temp(),y.temp(),kind));
}

cat::array<cat::tvector<double,4>,1> Spectral::eval_energ_spec(const CBVF & field,const bool & kind)
{
	cat::array<cat::tvector<double,4>,1> out(nwn);
	out[1]=eval_energ_spec(field.vel(),kind);
	out[2]=eval_energ_spec(field.mag(),kind);
	out[3]=eval_energ_spec(field.temp(),kind);
	for(int i=0;i<nwn;++i)
		out(i)[0]=i*wnstep;
	return out;
}

//block vector
void Spectral::pnvh_hat(const CBVF & field)
{
  cout << "vel" << endl;
	SpectralFourierLayer::pnvh_hat(field.vel());
  cout << "mag" << endl;
  SpectralFourierLayer::pnvh_hat(field.mag());
  cout << "temp" << endl;
  SpectralFourierLayer::pnvh_hat(field.temp());
}

