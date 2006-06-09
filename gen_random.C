
#include "input.h"

#include "globals.h"

#include "gen_random.h"

#include "block_vector.h"
#include "linops.h"

#include <iostream>
#include <string>
#include <fstream>

#include <cat.h>

#include "spectral.h"

using namespace std;
using namespace cat;

//Constructor (using timer value as seed)
gen_random::gen_random(const input & input_obj__,
                       spectral & spectral__):
random(),
spectral_obj(spectral__)
{
}

//Constructor (using a specified seed)
gen_random::gen_random(const input & input_obj__,
                       spectral & spectral__,
                       const int & seed__):
random(seed__),
spectral_obj(spectral__)
{
}

gen_random::~gen_random()
{
}

//generate random scalar field in fourier space
//generates random fourier coefficients for the required symmetry
void gen_random::gen_random_scalar_field_hat(CSF & field_hat,const int & ki, const int & kf, const bool & kind, const bool & sym)
{
	
	const int s1(field_hat.shape()[0]);
	const int s2(field_hat.shape()[1]);
	const int s3(field_hat.shape()[2]);
	
//Eval energy spectrum
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));
	double wvstep(sqrt(max(spectral_obj.wv2))/(energ_spec.size()-1));

	double dl=4./9.*max(spectral_obj.wv2);
	
  //Randomly generate field components in Fourier space
	//symmetry about the z axis is already partially imposed (see comment below)
	cat::array<Complex,3>::iterator field_hat_iterator(field_hat);
	cat::array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		int index=static_cast<int>(sqrt(*wv2_iterator)/wvstep);
		//if(index>=ki && index<kf)
		if((*wv2_iterator)<dl)
		{
			if(sym==1)
				*field_hat_iterator=complex<double>(random(-1.,1.),0.);
			else
				*field_hat_iterator=complex<double>(0.,random(-1.,1.));
		}
		else
			*field_hat_iterator=0.;
	}
	
	  
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
			field_hat(s1-i,0,k)=(sym?1:-1)*field_hat(i,0,k);
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
}




//generate random fields
//fourier coefficients for scalar field for sine/cossine representation
void gen_random::gen_random_field_hat(CSF & field_hat,
                                      const double & a,const double & b,
                                      const int & ki, const int & kf,
                                      const double & alpha,const double & p,
                                      const bool & kind,
                                      const bool & sym)
{
	const int s1(field_hat.shape()[0]);
	const int s2(field_hat.shape()[1]);
	const int s3(field_hat.shape()[2]);
	
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));
	double wvstep(sqrt(max(spectral_obj.wv2))/(energ_spec.size()-1));
	
	gen_random_scalar_field_hat(field_hat,ki,kf,kind,sym);
	
	 //Ensure that fields have zero average
	field_hat(0,0,0)=0.;
	
  //dealiasing
	spectral_obj.dealias(field_hat);
	
  //Eval spectrum between ki and kf
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
  //Normalise to spectrum between ki and kf
	cat::array<Complex,3>::iterator field_hat_iterator(field_hat);
	cat::array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		int index=static_cast<int>(sqrt(*wv2_iterator)/wvstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=	energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);
	}
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
  //Normalise for RMS vel,mag,temp = p
	//field_hat*=(p/sqrt(2*sum(energ_spec)*wvstep));
	double normfac=(p/sqrt((4*M_PI*M_PI*M_PI)*spectral_obj.scalar_prod(field_hat,field_hat,kind)));
	field_hat*=normfac;
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
}


//fourier coefficients for vector field
void gen_random::gen_random_field_hat
(CVF & field_hat,
 const double & a,const double & b,
 const int & ki, const int & kf,
 const double & alpha,const double & p,
 const bool & kind,const bool & sym)
{
	const int s1(field_hat.shape()[0]);
	const int s2(field_hat.shape()[1]);
	const int s3(field_hat.shape()[2]);
	
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));
	double wvstep(sqrt(max(spectral_obj.wv2))/(energ_spec.size()-1));
	
	CSF aux(field_hat.shape());
	aux=0;
	gen_random_scalar_field_hat(aux,ki,kf,!kind,!sym);
	field_hat[0]=0;
	aux=0;
	gen_random_scalar_field_hat(aux,ki,kf,!kind,!sym);
	field_hat[1]=aux;
	aux=0;
	gen_random_scalar_field_hat(aux,ki,kf,kind,sym);
	field_hat[2]=aux;
	
	 //Ensure that fields have zero average
	field_hat(0,0,0)=0.;
	
	//Make the field solenoidal
	spectral_obj.remove_gradient(field_hat,kind);
	
  //dealiasing
	spectral_obj.dealias(field_hat);
	
  //Eval spectrum between ki and kf
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	wvstep=sqrt(max(spectral_obj.wv2))/(energ_spec.size()-1);
	
  //Normalise to spectrum between ki and kf
	cat::array<cat::tvector<complex<Real>,3>,3>::iterator field_hat_iterator(field_hat);
	cat::array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		int index=static_cast<int>(sqrt(*wv2_iterator)/wvstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=	energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);
	}
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
  //Normalise for RMS vel,mag,temp = p
	//field_hat*=(p/sqrt(2*sum(energ_spec)));
	double normfac=(p/sqrt((4*M_PI*M_PI*M_PI)*spectral_obj.scalar_prod(field_hat,field_hat)));
	field_hat*=normfac;
	
//Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
}


//scalar field in real space
void gen_random::gen_random_field(RSF & field,
                                  const double & a,const double & b,
                                  const int & ki, const int & kf,
                                  const double & alpha,const double & p,
                                  const bool & kind,const bool & sym)
{
	const int s1(field.shape()[0]);
	const int s2(field.shape()[1]);
	const int s3(field.shape()[2]);
	
	CSF field_hat(s1,s2/2+1,s3);
	gen_random_field_hat(field_hat,a,b,ki,kf,alpha,p,kind,sym);
	
  //transform fields to real space
	if (kind)
		spectral_obj.sfft_c.inverse_transform(field,field_hat);
	else
		spectral_obj.sfft_s.inverse_transform(field,field_hat);
	
	//Test Symmetry
	
	
// 	for(int i=1;i<s1;++i)
// 		for(int j=1;j<s2/2+1;++j)
// 			for(int k=0;k<s3;++k)
// 				field(i,j,k)=(sym?1:-1)*field(s1-i,s2-j,k);
	
	
}

//vector field in real space
void gen_random::gen_random_field(RVF & field,
                                  const double & a,const double & b,
                                  const int & ki, const int & kf,
                                  const double & alpha,const double & p,
                                  const bool & kind,const bool & sym)
{
	const int s1(field.shape()[0]);
	const int s2(field.shape()[1]);
	const int s3(field.shape()[2]);
	
	CVF field_hat(s1,s2/2+1,s3);
	gen_random_field_hat(field_hat,a,b,ki,kf,alpha,p,kind,sym);
	
  //transform fields to real space
	if(kind)
		spectral_obj.fft_ssc.inverse_transform(field,field_hat);
	else
		spectral_obj.fft_ccs.inverse_transform(field,field_hat);
}
