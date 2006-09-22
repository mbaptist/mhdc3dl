
#include "linops.h" 

#include "globals.h"
#include "block_vector.h"
#include "input.h"
#include "basic.h"

#include <cat.h>
#include "spectral.h"

#include <complex>

#include <iostream>

using namespace std;
using namespace cat;


///////////////////////////////////////
//LINOPS_TMP
///////////////////////////////////////

linops_tmp::linops_tmp(input & input_obj__):  
n1(input_obj__.n1),
n2(input_obj__.n2),
n3(input_obj__.n3),
  //temporary variables
vel(n1,n2,n3),
mag(n1,n2,n3),
temp(n1,n2,n3),
curl_vel(n1,n2,n3),
curl_mag(n1,n2,n3),
grad_temp(n1,n2,n3),
curl_vel_hat(n1,n2/2+1,n3),
curl_mag_hat(n1,n2/2+1,n3),
grad_temp_hat(n1,n2/2+1,n3),
s_nonlinear(n1,n2,n3),
nonlinear(n1,n2,n3),
s_nonlinear_hat(n1,n2/2+1,n3),
nonlinear_hat(n1,n2/2+1,n3)
{
}

linops_tmp::~linops_tmp()
{
}

/////////////////////////////////////////////////////////////////////



///////////////////////////////////////////
//LINOPS_BASE
///////////////////////////////////////////

//Constructor
linops_base::linops_base(input & input_obj__,
                         Spectral & spectral_obj__,
                         Basic & basic__):
  //references for input parameters
n1(input_obj__.n1),
n2(input_obj__.n2),
n3(input_obj__.n3),
l1(input_obj__.l1),
l2(input_obj__.l2),
l3(input_obj__.l3),
visc(input_obj__.visc),
omegaz(input_obj__.omegaz),
compress(input_obj__.compress),
g(input_obj__.g),
diff(input_obj__.diff),
deltat(input_obj__.deltat),
econd(input_obj__.econd),
tcond(input_obj__.tcond),
  //reference to Spectral object
spectral_obj(spectral_obj__),
  //refrence to basic fields
basic(basic__),
sym_sub_(2),
  //references to temporary variables
  //temporaries are encapsulated in the singleton linops_tmp
  //and the appropriate references are created
vel(linops_tmp::instance(input_obj__).vel),
mag(linops_tmp::instance(input_obj__).mag),
temp(linops_tmp::instance(input_obj__).temp),
curl_vel(linops_tmp::instance(input_obj__).curl_vel),
curl_mag(linops_tmp::instance(input_obj__).curl_mag),
grad_temp(linops_tmp::instance(input_obj__).grad_temp),
curl_vel_hat(linops_tmp::instance(input_obj__).curl_vel_hat),
curl_mag_hat(linops_tmp::instance(input_obj__).curl_mag_hat),
grad_temp_hat(linops_tmp::instance(input_obj__).grad_temp_hat),
s_nonlinear(linops_tmp::instance(input_obj__).s_nonlinear),
nonlinear(linops_tmp::instance(input_obj__).nonlinear),
s_nonlinear_hat(linops_tmp::instance(input_obj__).s_nonlinear_hat),
nonlinear_hat(linops_tmp::instance(input_obj__).nonlinear_hat)
{
}


//Destructor
linops_base::~linops_base()
{
}

//Symmetry subspace accessor
int & linops_base::sym_sub()
{
	return sym_sub_;
}
const int & linops_base::sym_sub() const
{
	return sym_sub_;
}

//Scalar Product
Real linops_base::scalar_prod(const CBVF & xx,const CBVF & yy) const
{
	return (spectral_obj.scalar_prod(xx.vel(),yy.vel())
		+spectral_obj.scalar_prod(xx.mag(),yy.mag())
	        +spectral_obj.scalar_prod(xx.temp(),yy.temp()))/xx.size();
}

////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////
//A_NOUGHT
/////////////////////////////////////////////

//Ctor
a_nought::a_nought(input & input_obj__,
                   Spectral & spectral_obj__,
                   Basic & basic__):
linops_base(input_obj__,spectral_obj__,basic__)
{
}

//Dtor
a_nought::~a_nought()
{
}

//a_nought(w_hat)
CBVF a_nought::operator()(const CBVF & w_hat)
// const
{
	
  //create temporary block vector for output
	CBVF out(w_hat.shape());
	
  //transform fields into real space
	spectral_obj.fft_ccs.inverse_transform(vel,w_hat.vel());
	spectral_obj.fft_ccs.inverse_transform(mag,w_hat.mag());
	spectral_obj.sfft_s.inverse_transform(temp,w_hat.temp());
	
  //evaluate derivatives in Fourier space
	curl_vel_hat=spectral_obj.curl_hat(w_hat.vel(),0);
	curl_mag_hat=spectral_obj.curl_hat(w_hat.mag(),0);
	grad_temp_hat=spectral_obj.grad_hat(w_hat.temp(),0);
	
  //transform derivatives into real space
	spectral_obj.fft_ssc.inverse_transform(curl_vel,curl_vel_hat);
	spectral_obj.fft_ssc.inverse_transform(curl_mag,curl_mag_hat);
	spectral_obj.fft_ssc.inverse_transform(grad_temp,grad_temp_hat);
	
  //evaluate nonlinear term for velocity
	nonlinear=cross_product(basic.vel(),curl_vel)
		-cross_product(basic.curl_vel(),vel)
		-cross_product(basic.mag(),curl_mag)
		+cross_product(basic.curl_mag(),mag);
	
  //transform nonlinear into Fourier space
	spectral_obj.fft_ccs.direct_transform(nonlinear_hat,nonlinear);
	
  //evaluate velocity component
	CVF aux(out.vel().shape());
	aux=w_hat.temp();
	out.vel()=visc*spectral_obj.lap_hat(w_hat.vel())+compress*cat::tvector<Real,3>(0,0,g)*aux+nonlinear_hat;
	//Coriolis term
  //out.vel()+=(-cross_product(cat::tvector<Real,3>(0,0,omegaz),w_hat.vel()));
	
  //evaluate nonlinear term for magnetic field
	nonlinear=cross_product(basic.vel(),mag)
		-cross_product(basic.mag(),vel);
	
  //transform nonlinear into Fourier space
	spectral_obj.fft_ssc.direct_transform(nonlinear_hat,nonlinear);
	
  //magnetic field component
	out.mag()=diff*spectral_obj.lap_hat(w_hat.mag())+spectral_obj.curl_hat(nonlinear_hat,1);
	
  //evaluate nonlinear term (s_nonlinear) for temperature
	s_nonlinear=(-dot_product(vel,basic.grad_temp()));
  //Joule term
  //RSF rsaux(s_nonlinear.shape());
  //rsaux=dot_product(basic.curl_mag(),curl_mag);
  // rsaux*=econd;
  //s_nonlinear+=rsaux;
	s_nonlinear+=(-dot_product(basic.vel(),grad_temp));
	
  //transform s_nonlinear into Fourier space
	spectral_obj.sfft_s.direct_transform(s_nonlinear_hat,s_nonlinear);
	
  //temperature component
	out.temp()=-deltat*w_hat.vel()[2]+tcond*spectral_obj.lap_hat(w_hat.temp())+s_nonlinear_hat;
	
  //dealiasing
	spectral_obj.dealias(out);
	
  //remove gradient
	spectral_obj.remove_gradient(out,0);
	
	return out;
	
}


//A_NOUGHT_ADJOINT

//Ctor
a_nought_adjoint::a_nought_adjoint(input & input_obj__,
                                   Spectral & spectral_obj__,
                                   Basic & basic__):
linops_base(input_obj__,spectral_obj__,basic__)
{
}

//Dtor
a_nought_adjoint::~a_nought_adjoint()
{
}

//a_nought_adjoint(w_hat)
CBVF a_nought_adjoint::operator()(const CBVF & w_hat)
//const
{
	
  //create temporary block vector for output
	CBVF out(w_hat.shape());
	
  //transform fields into real space
	spectral_obj.fft_ccs.inverse_transform(vel,w_hat.vel());
	spectral_obj.fft_ccs.inverse_transform(mag,w_hat.mag());
	spectral_obj.sfft_s.inverse_transform(temp,w_hat.temp());
	
  //evaluate derivatives in Fourier space
	curl_vel_hat=spectral_obj.curl_hat(w_hat.vel(),0);
	curl_mag_hat=spectral_obj.curl_hat(w_hat.mag(),0);
	grad_temp_hat=spectral_obj.grad_hat(w_hat.temp(),0);
	
  //transform derivatives into real space
	spectral_obj.fft_ssc.inverse_transform(curl_vel,curl_vel_hat);
	spectral_obj.fft_ssc.inverse_transform(curl_mag,curl_mag_hat);
	spectral_obj.fft_ssc.inverse_transform(grad_temp,grad_temp_hat);
	
  //evaluate nonlinear term for velocity (first part)
	nonlinear=cross_product(basic.curl_vel(),vel)
		+cross_product(basic.mag(),curl_mag);
	
	RVF raux(nonlinear.shape());
	raux=basic.temp();
	raux*=grad_temp;
	nonlinear+=raux;
	
  //transform nonlinear into Fourier space (first part)
	spectral_obj.fft_ccs.direct_transform(nonlinear_hat,nonlinear);
	
  //evaluate velocity component (first part)
	
	CVF aux(out.vel().shape());
	aux=-w_hat.temp();
	aux*=cat::tvector<Real,3>(0,0,deltat);
	
	out.vel()=aux+visc*spectral_obj.lap_hat(w_hat.vel())+nonlinear_hat;
	
  //evaluate nonlinear term for velocity (second part)
	nonlinear=cross_product(basic.vel(),vel);
	
  //FFT nonlinear into Fourier space (second part)
	spectral_obj.fft_ssc.direct_transform(nonlinear_hat,nonlinear);
	
  //evaluate velocity component (second part)
	out.vel()+=(-spectral_obj.curl_hat(nonlinear_hat,1));
	
  //evaluate nonlinear term for magnetic field (first part)
	nonlinear=-cross_product(basic.curl_mag(),vel)
		-cross_product(basic.vel(),curl_mag);
  //Joule term
  //raux=-cross_product(basic.curl_mag(),grad_temp);
  //raux*=econd;
  //nonlinear+=raux;
  //raux=basic.curl_curl_mag();
  //raux*=temp;
  //raux*=econd;
  //nonlinear+=raux;
	
  //FFT nonlinear into Fourier space (first part)
	spectral_obj.fft_ccs.direct_transform(nonlinear_hat,nonlinear);
	
  //magnetic field component (first part)
	out.mag()=diff*spectral_obj.lap_hat(w_hat.mag())+nonlinear_hat;
	
  //evaluate nonlinear term for magnetic field (second part)
	nonlinear=cross_product(basic.mag(),vel);
	
  //FFT nonlinear into Fourier space (second part)
	spectral_obj.fft_ssc.direct_transform(nonlinear_hat,nonlinear);
	
  //magnetic field component (second part)
	out.mag()+=spectral_obj.curl_hat(nonlinear_hat,1);
	
	
	
  //evaluate nonlinear term (s_nonlinear) for temperature
	s_nonlinear=dot_product(basic.vel(),grad_temp);
	
  //FFT s_nonlinear into Fourier space
	spectral_obj.sfft_s.direct_transform(s_nonlinear_hat,s_nonlinear);
	
  //temperature component
	out.temp()=compress*g*w_hat.vel()[2]+tcond*spectral_obj.lap_hat(w_hat.temp())+s_nonlinear_hat;
	
  //dealiasing
	spectral_obj.dealias(out);
	
  //remove gradient
	spectral_obj.remove_gradient(out,0);
	
	return out;
	
}


//A_ONE

//Ctor
a_one::a_one(input & input_obj__,
             Spectral & spectral_obj__,
             Basic & basic__):
linops_base(input_obj__,spectral_obj__,basic__)
{
}

//Dtor
a_one::~a_one()
{
}


//a_one(index)
CBVF a_one::operator()(const CBVF & field,const int & index)
// const
{
	
	CBVF out(field.shape());
	
  //transform fields into real space
	spectral_obj.fft_ccs.inverse_transform(vel,field.vel());
	spectral_obj.fft_ccs.inverse_transform(mag,field.mag());
	spectral_obj.sfft_s.inverse_transform(temp,field.temp());
	
  //evaluate nonlinear term for velocity
  //nonlinear=unit_vector(index);
	nonlinear=0;
	nonlinear[index]=1;
	RSF rsaux(nonlinear.shape());
	rsaux=dot_product(basic.vel(),vel);
	rsaux+=-dot_product(basic.mag(),mag);
	nonlinear[index]*=rsaux;
	RVF raux(nonlinear.shape());
	raux=-vel;
	raux*=basic.vel()[index];
	nonlinear+=raux;
	raux=mag;
	raux*=basic.mag()[index];
	nonlinear+=raux;
	
  //FFT nonlinear into Fourier space
	spectral_obj.fft_ccs.direct_transform(nonlinear_hat,nonlinear);
	
  //evaluate velocity component
	out.vel()=spectral_obj.d_dx_index_hat(field.vel(),index);
	out.vel()*=visc;
	out.vel()*=2;
	out.vel()+=nonlinear_hat;
	
	
  //evaluate nonlinear term for magnetic field
	nonlinear=0;
	nonlinear=-basic.mag();
	nonlinear*=vel[index];
	raux=vel;
	raux*=basic.mag()[index];
	nonlinear+=raux;
	raux=basic.vel();
	raux*=mag[index];
	nonlinear+=raux;
	raux=-mag;
	raux*=basic.vel()[index];
	nonlinear+=raux;
	
  //FFT nonlinear into Fourier space
	spectral_obj.fft_ccs.direct_transform(nonlinear_hat,nonlinear);
	
  //evaluate magnetic component
	out.mag()=spectral_obj.d_dx_index_hat(field.mag(),index);
	out.mag()*=diff;
	out.mag()*=2;
	out.mag()+=nonlinear_hat;
	
	
  //evaluate nonlinear term (s_nonlinear) for temperature
  //Joule term
  //for(int o=0;o<3;++o)
  //  for(int k=0;k<3;++k)
  //    {
  //      rsaux=basic.curl_mag()[k];
  //      rsaux*=mag[o];
  //      rsaux*=levi_civita(index,o,k);
  //      s_nonlinear+=rsaux;
  //    }
  //s_nonlinear*=econd;
	rsaux=-temp;
	rsaux*=basic.vel()[index];
	s_nonlinear+=rsaux;
	
  //FFT s_nonlinear into Fourier space
	spectral_obj.sfft_s.direct_transform(s_nonlinear_hat,s_nonlinear);
	
  //evaluate temperature component
	out.temp()=spectral_obj.d_dx_index_hat(field.temp(),index);
	out.temp()*=tcond;
	out.temp()*=2;
	out.temp()+=s_nonlinear_hat;
	
  //dealiasing
	spectral_obj.dealias(out);
	
  //remove gradient
	spectral_obj.remove_gradient(out,0);
	
	return out;
	
}



////////////////////////////////
// Preconditioner
////////////////////////////////

precond::precond(Spectral & spectral_obj__,float qq__):
spectral_obj(spectral_obj__),
qq_(qq__)
{
}

precond::~precond()
{
}

CBVF precond::operator()
(const CBVF & x) const
{
	CBVF out(x.shape());
	array<CV,3>::iterator outvel_iterator(out.vel());
	array<CV,3>::iterator outmag_iterator(out.mag());
	array<CS,3>::iterator outtemp_iterator(out.temp());
	array<CV,3>::const_iterator xvel_iterator(x.vel());
	array<CV,3>::const_iterator xmag_iterator(x.mag());
	array<CS,3>::const_iterator xtemp_iterator(x.temp());
	array<RS,3>::const_iterator wv2_iterator(spectral_obj.wv2);
	for (outvel_iterator=out.vel().begin(),
	     outmag_iterator=out.mag().begin(),
	     outtemp_iterator=out.temp().begin(),
	     xvel_iterator=x.vel().begin(),
	     xmag_iterator=x.mag().begin(),
	     xtemp_iterator=x.temp().begin(),
	     wv2_iterator=spectral_obj.wv2.begin();
	     outvel_iterator!=out.vel().end();
	     ++outvel_iterator,
	     ++outmag_iterator,
	     ++outtemp_iterator,
	     ++xvel_iterator,
	     ++xmag_iterator,
	     ++xtemp_iterator,
	     ++wv2_iterator)
	{
		Real aux=pow((*wv2_iterator),-qq_);
		(*outvel_iterator)=(*xvel_iterator)*aux;
		(*outmag_iterator)=(*xmag_iterator)*aux;
		(*outtemp_iterator)=(*xtemp_iterator)*aux;
	}
	out.vel()(0,0,0)=x.vel()(0,0,0);
	out.mag()(0,0,0)=x.mag()(0,0,0);
	out.temp()(0,0,0)=x.temp()(0,0,0);
	return out;
}

