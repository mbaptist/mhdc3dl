
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

<<<<<<< .mine
//generate random fields
//fourier coefficients for scalar field for sine/cossine representation
void gen_random::gen_random_field_hat(CSF & field_hat,
                                      const int & ki, const int & kf,
                                      const double & alpha,const double & p,
                                      const bool & kind,
                                      const bool & sym)
=======
void gen_random::gen_random_fcoeffs_scalar_field_hat(CSF & field_hat)
>>>>>>> .r49
{
<<<<<<< .mine

	const int s1(field_hat.shape()[0]);
	const int s2(field_hat.shape()[1]);
	const int s3(field_hat.shape()[2]);
	
  //Randomly generate field components in Fourier space
	//symmetry about the z axis is already partially imposed (see comment below)
=======
	gen_random_fourier_coeffs(field_hat,0,spectral_obj.max_shell);
}

void gen_random::gen_random_fourier_coeffs(CSF & field_hat,const int & initial_shell, const int & final_shell)
{
>>>>>>> .r49
	cat::array<Complex,3>::iterator field_hat_iterator(field_hat);
	cat::array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
<<<<<<< .mine
		int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
=======
>>>>>>> .r49
		if(index>=ki && index<kf)
				*field_hat_iterator=complex<double>(random(-1.,1.),random(-1.,1.);
		else
			*field_hat_iterator=0.;
	}
<<<<<<< .mine
  //symmetry about z axis
=======
}

//Applies symmetry about the z-axis to a scalar field
void gen_random::sym_scalar_field_hat(CSF & field_hat, const bool & kind, const bool & sym)
{
	// Imaginary(or real) parts must vanish
	cat::array<Complex,3>::iterator field_hat_iterator(field_hat);
	cat::array<Real,3>::iterator wv2_iterator(spectral_obj.wv2);
	for(field_hat_iterator=field_hat.begin(),
	    wv2_iterator=spectral_obj.wv2.begin();
	    field_hat_iterator!=field_hat.end(),
	    wv2_iterator!=spectral_obj.wv2.end();
	    ++field_hat_iterator,
	    ++wv2_iterator)
	{
		if(sym==1)
			field_hat_iterator->imag()=0;
		else
			field_hat_iterator->real()=0;
} 
   //symmetry about z axis
>>>>>>> .r49
	//combining symmetry about the z axis and hermitian symmetry, we obtain that
	//symetric fields are real and anti-symetric fields are imaginary;
	//this condition is imposed above;
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
	gen_random_fcoeffs_scalar_field_hat(field);
	sym_scalar_field_hat(field,kind,sym);
	 //Ensure that fields have zero average
	field_hat(0,0,0)=0.;
  //dealiasing
	spectral_obj.dealias(field_hat);
	//Normalise to spectrum 
	spectral_obj.normalise_to_spectrum(field,spectrum_function);
	
<<<<<<< .mine
  //Eval spectrum between ki and kf
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));
	
=======
>>>>>>> .r49
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
<<<<<<< .mine
		int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);
=======
		int index=static_cast<int>(sqrt(*wv2_iterator)/wvstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=	energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);
>>>>>>> .r49
	}
		
// //renormalise to RMS(field)=p - in fact for total energy equal to p/2
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind); //re-eval spectrum
	double tenerg=sum(energ_spec)*spectral_obj.wnstep;
	field_hat*=sqrt(p/(2.*tenerg));
	
<<<<<<< .mine
=======
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
  //Normalise for RMS vel,mag,temp = p
	//field_hat*=(p/sqrt(2*sum(energ_spec)*wvstep));
	double normfac=(p/sqrt((4*M_PI*M_PI*M_PI)*spectral_obj.scalar_prod(field_hat,field_hat,kind)));
	field_hat*=normfac;
	
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
>>>>>>> .r49
}


//fourier coefficients for vector field
void gen_random::gen_random_field_hat
(CVF & field_hat,
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
<<<<<<< .mine
gen_random_field_hat(aux,ki,kf,alpha,p,(kind?0:1),(sym?0:1));
=======
	gen_random_scalar_field_hat(aux,ki,kf,!kind,!sym);
>>>>>>> .r49
	field_hat[0]=aux;
	aux=0;
<<<<<<< .mine
gen_random_field_hat(aux,ki,kf,alpha,p,(kind?0:1),(sym?0:1));
=======
	gen_random_scalar_field_hat(aux,ki,kf,!kind,!sym);
>>>>>>> .r49
	field_hat[1]=aux;
	aux=0;
<<<<<<< .mine
gen_random_field_hat(aux,ki,kf,alpha,p,(kind?1:0),(sym?1:0));
=======
	gen_random_scalar_field_hat(aux,ki,kf,kind,sym);
>>>>>>> .r49
	field_hat[2]=aux;
	
	 //Ensure that fields have zero average
	field_hat(0,0,0)=0.;
	
	//Make the field solenoidal
	spectral_obj.remove_gradient(field_hat,kind);
	
  //dealiasing
	spectral_obj.dealias(field_hat);
	
  //Eval spectrum between ki and kf
<<<<<<< .mine
	cat::array<Real,1> energ_spec(spectral_obj.eval_energ_spec(field_hat,kind));
=======
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	wvstep=sqrt(max(spectral_obj.wv2))/(energ_spec.size()-1);
>>>>>>> .r49
	
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
<<<<<<< .mine
		int index=static_cast<int>(sqrt(*wv2_iterator)/spectral_obj.wnstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);
=======
		int index=static_cast<int>(sqrt(*wv2_iterator)/wvstep);
		double power=pow(sqrt(*wv2_iterator),-alpha);
		double espec=	energ_spec(index);
		if(espec!=0)
			(*field_hat_iterator)*=sqrt(power/espec);
>>>>>>> .r49
	}
	
<<<<<<< .mine
		//renormalise to RMS(field)=p - in fact for total energy equal to p/2
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind); //re-eval spectrum
	double tenerg=sum(energ_spec)*spectral_obj.wnstep;
	field_hat*=sqrt(p/(2.*tenerg));
=======
  //Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
>>>>>>> .r49
	
<<<<<<< .mine
=======
  //Normalise for RMS vel,mag,temp = p
	//field_hat*=(p/sqrt(2*sum(energ_spec)));
 	double normfac=(p/sqrt((4*M_PI*M_PI*M_PI)*spectral_obj.scalar_prod(field_hat,field_hat)));
	field_hat*=normfac;
	
//Re-eval energ_spec
	energ_spec=spectral_obj.eval_energ_spec(field_hat,kind);
	
>>>>>>> .r49
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
