
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

//generate random fields
//fourier coefficients for scalar field for sine/cossine representation
void gen_random::gen_random_field_hat(CSF & field_hat,
                                      const int & ki, const int & kf,
                                      const double & alpha,const double & p,
                                      const bool & kind,
                                      const bool & sym)
{
	const int s1(field_hat.shape()[0]);
	const int s2(field_hat.shape()[1]);
	const int s3(field_hat.shape()[2]);
	
  //Eval energy spectrum
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat));
	double wv2step(max(spectral_obj.wv2)/(energ_spec.size()-1));
	
  //Randomly generate field components in Fourier space
	cat::array_iterator<Complex,3>
		field_hat_iterator(field_hat);
	cat::array_iterator<Real,3> wv2_iterator(spectral_obj.wv2);
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		int index=static_cast<int>((*wv2_iterator)/wv2step);
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
	
  //Ensure that fields have zero average
	field_hat(0,0,0)=0.;
	
  //Ensure that kz=0 terms are zero in sine representation
	if (kind==0)
		for(int i=0;i<s1;++i)
			for(int j=0;j<s2;++j)
			{
				field_hat(i,j,0)=0.;
				field_hat(i,j,s3-1)=0.;
			}
	
  //symmetry about z axis
	for(int i=1;i<s1/2+1;++i)
		//for(int j=0;j<s2/2+1;++j)
		for(int k=0;k<s3;++k)
			field_hat(s1-i,0,k)=(sym?1:-1)*field_hat(i,0,k);
	
  //dealiasing
	spectral_obj.dealias(field_hat);
	
  //Eval spectrum between ki and kf
	energ_spec=spectral_obj.eval_energ_spec(field_hat);
	
  //Normalise to spectrum between ki and kf
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		int index=static_cast<int>((*wv2_iterator)/wv2step);
		if(index>=ki && index<kf)
		{
			double power=pow(double(*wv2_iterator),-alpha);
			(*field_hat_iterator)*=sqrt(power/energ_spec(index));
		}
	}
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat);
	
  //Normalise for RMS vel,mag,temp = p
	field_hat*=(p/sqrt(2*sum(energ_spec)));
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat);	
	
}


//fourier coefficients for vector field
void gen_random::gen_random_field_hat
(CVF & field_hat,
 const int & ki, const int & kf,
 const double & alpha,const double & p,
 const bool & kind,const bool & sym)
{
	CSF aux(field_hat.shape());
	aux=0;
gen_random_field_hat(aux,ki,kf,alpha,p,(kind?0:1),(sym?0:1));
	field_hat[0]=aux;
	aux=0;
gen_random_field_hat(aux,ki,kf,alpha,p,(kind?0:1),(sym?0:1));
	field_hat[1]=aux;
	aux=0;
gen_random_field_hat(aux,ki,kf,alpha,p,(kind?1:0),(sym?1:0));
	field_hat[2]=aux;
	
	spectral_obj.remove_gradient(field_hat,kind);
	
  //dealiasing
	spectral_obj.dealias(field_hat);
	
  //Eval spectrum between ki and kf
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat));
	double wv2step(max(spectral_obj.wv2)/(energ_spec.size()-1));
	
  //Normalise to spectrum between ki and kf
	cat::array_iterator<cat::tvector<complex<Real>,3>,3>
		field_hat_iterator(field_hat);
	cat::array_iterator<Real,3> wv2_iterator(spectral_obj.wv2);
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		int index=static_cast<int>((*wv2_iterator)/wv2step);
		if(index>=ki && index<kf)
		{
			double power=pow(double(*wv2_iterator),-alpha);
			(*field_hat_iterator)*=sqrt(power/energ_spec(index));
		}
	}
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat);
	
  //Normalise for RMS vel,mag,temp = p
	field_hat*=(p/sqrt(2*sum(energ_spec)));
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat);
	
}


//scalar field in real space
void gen_random::gen_random_field(RSF & field,
                                  const int & ki, const int & kf,
                                  const double & alpha,const double & p,
                                  const bool & kind,const bool & sym)
{
	const int s1(field.shape()[0]);
	const int s2(field.shape()[1]);
	const int s3(field.shape()[2]);
	
	CSF field_hat(s1,s2/2+1,s3);
	gen_random_field_hat(field_hat,ki,kf,alpha,p,kind,sym);
	
  //transform fields to real space
	if (kind)
		spectral_obj.sfft_c.inverse_transform(field,field_hat);
	else
		spectral_obj.sfft_s.inverse_transform(field,field_hat);

	//Test Symmetry

// 	for(int i=1;i<s1/2+1;++i)
// 		for(int j=1;j<s2/2+1;++j)
// 			for(int k=0;k<s3;++k)
// 			{
// 				cout << i << j <<k <<endl;
// 					assert(field(s1-i,s2-j,k)==(sym?1:-1)*field(i,j,k));
// 			}

	for(int i=1;i<s1/2+1;++i)
		for(int j=0;j<s2;++j)
			for(int k=0;k<s3;++k)
				field(s1-i,j,k)=(sym?1:-1)*field(i,j,k);
	for(int i=0;i<s1;++i)
		for(int j=1;j<s2/2+1;++j)
			for(int k=0;k<s3;++k)
				field(i,s2-j,k)=(sym?1:-1)*field(i,j,k);
				
}

//vector field in real space
void gen_random::gen_random_field(RVF & field,
                                  const int & ki, const int & kf,
                                  const double & alpha,const double & p,
                                  const bool & kind,const bool & sym)
{
	const int s1(field.shape()[0]);
	const int s2(field.shape()[1]);
	const int s3(field.shape()[2]);
	
	CVF field_hat(s1,s2/2+1,s3);
	gen_random_field_hat(field_hat,ki,kf,alpha,p,kind,sym);
	
  //transform fields to real space
	if(kind)
		spectral_obj.fft_ssc.inverse_transform(field,field_hat);
	else
		spectral_obj.fft_ccs.inverse_transform(field,field_hat);
}
