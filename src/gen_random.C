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


#include "input.h"

#include "globals.h"

#include "gen_random.h"

#include "block_vector.h"
#include "linops.h"

#include <iostream>
#include <string>
#include <fstream>

#include <iomanip>

#include <cat.h>

#include "spectral.h"

using namespace std;
using namespace cat;

//Constructor (using timer value as seed)
gen_random::gen_random(const Input & input_obj__,
                       Spectral & spectral__):
    random(),
    input_obj(input_obj__),
    spectral_obj(spectral__)
{}

//Constructor (using a specified seed)
gen_random::gen_random(const Input & input_obj__,
                       Spectral & spectral__,
                       const int & seed__):
    random(seed__),
    input_obj(input_obj__),
    spectral_obj(spectral__)
{}

gen_random::~gen_random()
{}

//generate random fields
//fourier coefficients for scalar field for sine/cossine representation
void gen_random::gen_random_field_hat(CSF & field_hat,
                                      const int & ki, const int & kf,
                                      const double & alpha,const double & p,
                                      const bool & kind,
                                      const bool & sym,
                                      const std::string & br_spectrum)
{

  //	cout << ki << kf << kind <<sym << endl;
  //cout << spectral_obj.wnmax << endl;
  //cout << spectral_obj.wnstep << endl;

  const int s1(field_hat.shape()[0]);
  const int s2(field_hat.shape()[1]);
  const int s3(field_hat.shape()[2]);

  //Randomly generate field components in Fourier space
  //symmetry about the z axis is already partially imposed (see comment below)
  cat::Array<Complex,3>::iterator field_hat_iterator(field_hat);
  cat::Array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
  cat::Array<cat::Tvector<Real,3>,3>::iterator wv_iterator(spectral_obj.wv);
  for(field_hat_iterator=field_hat.begin(),
      wv2_iterator=spectral_obj.wv2.begin(),
      wv_iterator=spectral_obj.wv.begin();
      field_hat_iterator!=field_hat.end(),
      wv2_iterator!=spectral_obj.wv2.end(),
      wv_iterator!=spectral_obj.wv.end();
      ++field_hat_iterator,
      ++wv2_iterator,
      ++wv_iterator)
    {
      int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
      //cout << index << " " << *wv2_iterator  << " " << *wv_iterator << endl;
      if(index>=ki && index<kf)
        {
          if(sym==1)
            *field_hat_iterator=complex<double>(random(-1.,1.),0.);
          else
            *field_hat_iterator=complex<double>(0.,random(-1.,1.));
        }
      else
        *field_hat_iterator=0.;
    }

  //cout	<< "pnvh" << endl;
  //spectral_obj.pnvh_hat(field_hat);

  //symmetry about z axis
  //combining symmetry about the z axis and hermitian symmetry, we obtain that
  //symetric fields are real and anti-symetric fields are imaginary;
  //this condition is imposed above, as we generated the coefficients;
  //any two of the above three conditions specify the symmetry completely;
  //since only half of the harmonics are stored, 2 conditions are being used, except
  //in the plane ky=0; therefore, we impose below the condition of symmetry
  //about the z axis for that plane.
  for(int i=1;i<s1/2+1;++i)
    for(int k=0;k<s3;++k)
      {
        field_hat(s1-i,0,k)=field_hat(i,0,k);
        field_hat(s1-i,0,k)*=(sym?1:-1);
        //field_hat(s1-i,0,k)=(sym?1:-1)*field_hat(i,0,k);
      }
  //for anti-simmetric fields the values at the line kx=0, ky=0 must all be zero.
  if (sym==0)
    for(int k=0;k<s3;++k)
      field_hat(0,0,k)=0;
  //Ensure that kz=0 terms are zero in sine representation
  if (kind==0)
    for(int i=0;i<s1;++i)
      for(int j=0;j<s2;++j)
        {
          field_hat(i,j,0)=0.;
          field_hat(i,j,s3-1)=0.;
        }

  //Ensure that fields have zero average
  field_hat(0,0,0)=0.;

  //dealiasing
  spectral_obj.dealias(field_hat);

  //Eval spectrum between ki and kf
  cat::Array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));

  //cout << setprecision(20) << "Before" << "\n" << energ_spec << endl;

  // 	field_hat*=2;
  //
  // 	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind); //re-eval spectrum
  //
  // 	cout << "After" << "\n" << energ_spec << endl;
  //
  // 	exit(0);

  //Normalise to spectrum between ki and kf
  for(field_hat_iterator=field_hat.begin(),
      wv2_iterator=spectral_obj.wv2.begin();
      field_hat_iterator!=field_hat.end(),
      wv2_iterator!=spectral_obj.wv2.end();
      ++field_hat_iterator,
      ++wv2_iterator)
    {
      int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
      if(index<spectral_obj.nwn)
        {
          if(energ_spec(index)!=0)
            {
              //double factor=sqrt(pow(sqrt(*wv2_iterator),-alpha)/energ_spec(index));
              //double factor=sqrt(pow((*wv2_iterator),-alpha/2.)/energ_spec(index));


              double factor;
              if (br_spectrum=="power")
                factor=sqrt(pow(index,-alpha)/energ_spec(index));
              else if (br_spectrum=="exp")
                factor=sqrt(exp(-alpha*index)/energ_spec(index));

              //double factor=sqrt(exp(-index*alpha)/energ_spec(index));

              //cout << "||wv||^2=" << (*wv2_iterator) << " ||wv||=" << sqrt(*wv2_iterator) << "SphSh= " << index << " factor=" << factor << endl;
              (*field_hat_iterator)*=factor;
            }
        }
    }

  //dealiasing
  spectral_obj.dealias(field_hat);

  // //renormalise to RMS(field)=p - in fact for total energy equal to p/2
  energ_spec=spectral_obj.eval_energ_spec(field_hat,kind); //re-eval spectrum

  //cout << "After" << "\n" << energ_spec << endl;
  //exit(0);

  double av_energ=sum(energ_spec)*spectral_obj.wnstep;
  field_hat*=p/sqrt(2.*av_energ);

}


//fourier coefficients for vector field
void gen_random::gen_random_field_hat
(CVF & field_hat,
 const int & ki, const int & kf,
 const double & alpha,const double & p,
 const bool & kind,const bool & sym,
 const std::string & br_spectrum)
{
  CSF aux(field_hat.shape());
  aux=0;
  gen_random_field_hat(aux,ki,kf,alpha,p,(kind?0:1),(sym?0:1),br_spectrum);
  field_hat[0]=aux;
  aux=0;
  gen_random_field_hat(aux,ki,kf,alpha,p,(kind?0:1),(sym?0:1),br_spectrum);
  field_hat[1]=aux;
  aux=0;
  gen_random_field_hat(aux,ki,kf,alpha,p,(kind?1:0),(sym?1:0),br_spectrum);
  field_hat[2]=aux;

  //Make the field solenoidal
  spectral_obj.remove_gradient(field_hat,kind);

  //dealiasing
  spectral_obj.dealias(field_hat);

  //Eval spectrum between ki and kf
  cat::Array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));

  //Normalise to spectrum between ki and kf
  cat::Array<cat::Tvector<complex<Real>,3>,3>::iterator field_hat_iterator(field_hat);
  cat::Array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
  for(field_hat_iterator=field_hat.begin(),
      wv2_iterator=spectral_obj.wv2.begin();
      field_hat_iterator!=field_hat.end(),
      wv2_iterator!=spectral_obj.wv2.end();
      ++field_hat_iterator,
      ++wv2_iterator)
    {
      int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
      if(index<spectral_obj.nwn)
        {
          if(energ_spec(index)!=0)
            {
              //double factor=sqrt(pow(sqrt(*wv2_iterator),-alpha)/energ_spec(index));
              //double factor=sqrt(pow((*wv2_iterator),-alpha/2.)/energ_spec(index));

              double factor;
              if (br_spectrum=="power")
                factor=sqrt(pow(index,-alpha)/energ_spec(index));
              else if (br_spectrum=="exp")
                factor=sqrt(exp(-alpha*index)/energ_spec(index));

              //double factor=sqrt(exp(-index*alpha)/energ_spec(index));

              //cout << "||wv||^2=" << (*wv2_iterator) << " ||wv||=" << sqrt(*wv2_iterator) << "SphSh= " << index << " factor=" << factor << endl;
              (*field_hat_iterator)*=factor;
            }
        }


    }

  //dealiasing
  spectral_obj.dealias(field_hat);

  //renormalise to RMS(field)=p
  energ_spec=spectral_obj.eval_energ_spec(field_hat,kind); //re-eval spectrum
  double av_energ=sum(energ_spec)*spectral_obj.wnstep;
  field_hat*=p/sqrt(2.*av_energ);
}


//scalar field in real space
void gen_random::gen_random_field(RSF & field,
                                  const int & ki, const int & kf,
                                  const double & alpha,const double & p,
                                  const bool & kind,const bool & sym,
                                  const std::string & br_spectrum)
{
  const int s1(field.shape()[0]);
  const int s2(field.shape()[1]);
  const int s3(field.shape()[2]);

  CSF field_hat(s1,s2/2+1,s3);
  gen_random_field_hat(field_hat,ki,kf,alpha,p,kind,sym,br_spectrum);

  //transform fields to real space
  if (kind)
    spectral_obj.sfft_c.inverse_transform(field,field_hat);
  else
    spectral_obj.sfft_s.inverse_transform(field,field_hat);

}

//vector field in real space
void gen_random::gen_random_field(RVF & field,
                                  const int & ki, const int & kf,
                                  const double & alpha,const double & p,
                                  const bool & kind,const bool & sym,
                                  const std::string & br_spectrum)
{
  const int s1(field.shape()[0]);
  const int s2(field.shape()[1]);
  const int s3(field.shape()[2]);

  CVF field_hat(s1,s2/2+1,s3);
  gen_random_field_hat(field_hat,ki,kf,alpha,p,kind,sym,br_spectrum);

  //transform fields to real space
  if(kind)
    spectral_obj.fft_ssc.inverse_transform(field,field_hat);
  else
    spectral_obj.fft_ccs.inverse_transform(field,field_hat);

}
