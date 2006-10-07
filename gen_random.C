
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
                       Spectral & spectral__):
random(),
spectral_obj(spectral__)
{
}

//Constructor (using a specified seed)
gen_random::gen_random(const input & input_obj__,
                       Spectral & spectral__,
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

	//	cout << ki << kf << kind <<sym << endl;
	//cout << spectral_obj.wnmax << endl;
	//cout << spectral_obj.wnstep << endl;
		
	const int s1(field_hat.shape()[0]);
	const int s2(field_hat.shape()[1]);
	const int s3(field_hat.shape()[2]);
	
  //Randomly generate field components in Fourier space
	//symmetry about the z axis is already partially imposed (see comment below)
	cat::array<Complex,3>::iterator field_hat_iterator(field_hat);
	cat::array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
	cat::array<cat::tvector<Real,3>,3>::iterator wv_iterator(spectral_obj.wv);
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
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));
	
  //Normalise to spectrum between ki and kf
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);
	}
		
// //renormalise to RMS(field)=p - in fact for total energy equal to p/2
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind); //re-eval spectrum
	double tenerg=sum(energ_spec)*spectral_obj.wnstep;
	field_hat*=sqrt(p/(2.*tenerg));
	
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

	//Make the field solenoidal
	spectral_obj.remove_gradient(field_hat,kind);
	
  //dealiasing
	spectral_obj.dealias(field_hat);
	
  //Eval spectrum between ki and kf
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));
	
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
		int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);	
	}
	
		//renormalise to RMS(field)=p - in fact for total energy equal to p/2
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind); //re-eval spectrum
	double tenerg=sum(energ_spec)*spectral_obj.wnstep;
	field_hat*=sqrt(p/(2.*tenerg));

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
