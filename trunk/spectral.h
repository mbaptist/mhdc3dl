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





#ifndef SPECTRAL_H
#define SPECTRAL_H

#include "globals.h"

#include <goops.h>

#include <cat.h>


class Spectral : public SpectralFourierLayer
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

    Spectral(const int & n1__,const int & n2__,const int & n3__,
             const Real & l1__,const Real & l2__,const Real & l3__);

    ~Spectral();

    void dealias(CBVF & field)  const;//BlockVector

    //Laplacian of BlockVector
    CBVF lap_hat(const CBVF & field);

    //remove gradient part after of a BlockVector
    CBVF remove_gradient(CBVF & bfield,const bool kind);

    Real scalar_prod(const CBVF & xx,const CBVF & yy,const bool & kind) const ;

    cat::Array<cat::Tvector<double,4>,1> eval_energ_spec(const CBVF & field,const bool & kind);

    void pnvh_hat(const CBVF & field);

  };

#endif
