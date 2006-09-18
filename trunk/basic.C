
#include "input.h"

#include "globals.h"

#include "basic.h"


#include "block_vector.h"
#include "linops.h"

#include "gen_random.h"

#include <string>
#include <fstream>

#include <cat.h>
#include "spectral.h"

#include "vtkio.h"


using namespace std;
using namespace cat;


////////////////////////////////
// Class Basic definition
////////////////////////////////

//Accessors
RVF & Basic::vel()
{
	return vel_;
}
RVF & Basic::mag()
{	
	return mag_;
}
RSF & Basic::temp()
{
	return temp_;
}
RVF & Basic::curl_vel()
{
	return curl_vel_;
}
RVF & Basic::curl_mag()
{
	return curl_mag_;
}
RVF & Basic::curl_curl_mag()
{
	return curl_curl_mag_;
}
RVF & Basic::grad_temp()
{
	return grad_temp_;
}

//Constructors
Basic::Basic(const input & input_obj__,Spectral & spectral_obj__):
input_obj(input_obj__),
spectral_obj(spectral_obj__),
n1(input_obj__.n1),
n2(input_obj__.n2),
n3(input_obj__.n3),
vel_(n1,n2,n3),
mag_(n1,n2,n3),
temp_(n1,n2,n3),
curl_vel_(n1,n2,n3),
curl_mag_(n1,n2,n3),
curl_curl_mag_(n1,n2,n3),
grad_temp_(n1,n2,n3)
{
  //Load from files
	if (input_obj.basic_mode=="load")
	{
		//cout << "mode=load" << endl;
		load(input_obj.basic_fname);
	}
  //Randomly generated
	else if (input_obj.basic_mode=="random")
	{
		gen_random * gen_random_obj;
		int seed=input_obj.br_seed;
      //instantiate random generator
		if (seed>0)
			gen_random_obj=new gen_random(input_obj,spectral_obj,seed);
		else
			gen_random_obj=new gen_random(input_obj,spectral_obj);
      //generate fields
		//cout << input_obj.br_sym << endl;
		gen_random_obj->gen_random_field(vel_,input_obj.br_ki,input_obj.br_kf,input_obj.br_alpha,input_obj.br_rms_norm,input_obj.br_kind,input_obj.br_sym);
		gen_random_obj->gen_random_field(mag_,input_obj.br_ki,input_obj.br_kf,input_obj.br_alpha,input_obj.br_rms_norm,input_obj.br_kind,input_obj.br_sym);
		gen_random_obj->gen_random_field(temp_,input_obj.br_ki,input_obj.br_kf,input_obj.br_alpha,input_obj.br_rms_norm,input_obj.br_kind,input_obj.br_sym);
	cout	<< "pnvh vel" << endl;
		spectral_obj.pnvh(this->vel(),0);
		cout	<< "pnvh mag" << endl;
		spectral_obj.pnvh(this->mag(),0);
		cout	<< "pnvh temp" << endl;
		spectral_obj.pnvh(this->temp(),0);
		
		delete gen_random_obj;
		eval_derivatives();
	}
  //Plan forms
	else if (input_obj.basic_mode=="plan")
	{
		Real l1(input_obj.l1);
		Real l2(input_obj.l2);
		Real l3(input_obj.l3);
		
		vel_=0;
		temp_=0;
		mag_=0;
		
		double alpha1=1.;
		double alpha2=0.;
		double ell2=2*M_PI/l2;
		double ell1=2*M_PI/l1;
		double gamma=2;
		double gm=gamma*gamma+1;
		double pp=sqrt(gm);
		
		RVF pot(n1,n2,n3);
		for(int i=0;i<n1;++i)
			for(int j=0;j<n2;++j)
			{
				for(int k=0;k<n3;++k)
				{
					pot(i,j,k)=cat::tvector<double,3>
						(
						  0,
						  0,
						  (
						    alpha1*cos(ell1*i*l1/n1)*cos(ell2*j*l2/n2)
						    +alpha2*cos(pp*ell2*j*l2/n2)
						  ) * (
						        sin(k*l3/(n3-1))+sin(2*k*l3/(n3-1))
						      )
						);
				}
			}
		
		CVF pot_hat(n1,n2/2+1,n3);
		spectral_obj.fft_ccs.direct_transform(pot_hat,pot);
		
      //   spectral_obj.pnvh(pot_hat);
      //   exit(0);
		
		
		CVF aux_hat(n1,n2/2+1,n3);
		aux_hat=spectral_obj.curl_hat(pot_hat,0);
		CVF vv_hat(n1,n2/2+1,n3);
		vv_hat=spectral_obj.curl_hat(aux_hat,1);
		
		
      //  double mw=2;
      // double vn=sqrt(8./((2.*beta*beta+pow((alpha*mw),2))*(1.+pow((ell2*mw),2))));
		
		double mz=2;
		double vn=sqrt(8./(gm*(1.+gm*ell2*ell2+alpha1*alpha1*(mz*mz+gm*ell2*ell2))));
		
		vv_hat*=(vn/8.);
		vv_hat*=2.; //My basis is different from Vlad's
		
      //vv_hat*=(-1.);
		
		spectral_obj.fft_ccs.inverse_transform(vel_,vv_hat);
		
      //cout << spectral_obj.scalar_prod(vv_hat,vv_hat) << endl;
		
		cout << "Non-vanishing harmonics in Basic velocity" << endl;
	    //      spectral_obj.pnvh(vel_,0);
		
	cout << "Energy spectrum of Basic velocity: " << endl;
		
	    //cout << spectral_obj.eval_energ_spec(vv_hat) << endl;
		
	cout << "s_prod_hat: " << sqrt(spectral_obj.scalar_prod(vv_hat,vv_hat)) << endl;
		
      //Evaluate derivatives
		eval_derivatives();
	}
  //Expression
	else if (input_obj.basic_mode=="expression")
	{
		int n1(input_obj.n1);
		int n2(input_obj.n2);
		int n3(input_obj.n3);
		Real l1(input_obj.l1);
		Real l2(input_obj.l2);
		Real l3(input_obj.l3);
		for(int i=0;i<n1;++i)
			for(int j=0;j<n2;++j)
				for(int k=0;k<n3;++k)
				{
					double x=i*l1/n1;
					double y=j*l2/n2;
					double z=k*l3/(n3-1);
// 	      vel_(i,j,k)=cat::tvector<double,3>(sin(x)*cos(z),
// 						 sin(y)*cos(z),
// 						 -(cos(x)+cos(y))*sin(z));
// 	      mag_(i,j,k)=cat::tvector<double,3>(sin(x+y)*cos(z),
// 						 sin(x+y)*cos(z),
// 						 (-2.*cos(x+y))*sin(z));
// 	      temp_(i,j,k)=sin(z);
					vel_(i,j,k)=cat::tvector<double,3>(sin(10*x)*cos(10*z),
					                                   sin(10*y)*cos(10*z),
					                                   -(cos(10*x)+cos(10*y))*sin(10*z));
					mag_(i,j,k)=cat::tvector<double,3>(sin(10*x+10*y)*cos(10*z),
					                                   sin(10*x+10*y)*cos(10*z),
					                                   (-2.*cos(10*x+10*y))*sin(10*z));
					temp_(i,j,k)=sin(10*z);
				}
      //Evaluate derivatives
		eval_derivatives();
	}
}

//Destructor
Basic::~Basic()
{
}

void Basic::load(const string & vel_fname,
                 const string & mag_fname,
                 const string & temp_fname)
{
  //Basic velocity
	ifstream ifs(vel_fname.c_str());
	ifs >> vel_;
	ifs.close();
  //Basic magnetic field
	ifs.open(mag_fname.c_str());
	ifs >> mag_;
	ifs.close();
  //Basic temperature
	ifs.open(temp_fname.c_str());
	ifs >> temp_;
	ifs.close();
  //evaluate derivatives
	eval_derivatives();
}

void Basic::save(const string & vel_fname,
                 const string & mag_fname,
                 const string & temp_fname)
{
  //Basic velocity
	ofstream ofs(vel_fname.c_str());
	ofs << vel_ << endl;
	ofs.close();
  //Basic magnetic field
	ofs.open(mag_fname.c_str());
	ofs << mag_ << endl;
	ofs.close();
  //Basic temperature
	ofs.open(temp_fname.c_str());
	ofs << temp_ << endl;
	ofs.close();
}

void Basic::eval_derivatives()
{
  //curl_vel
	CVF aux(n1,n2/2+1,n3);
	CVF aux2(n1,n2/2+1,n3);
	aux=0;
	aux2=0;
	spectral_obj.fft_ccs.direct_transform(aux,vel_);
	aux2=spectral_obj.curl_hat(aux,0);
	spectral_obj.fft_ssc.inverse_transform(curl_vel_,aux2);
  //curl_mag
	aux=0;
	aux2=0;
	spectral_obj.fft_ccs.direct_transform(aux,mag_);
	aux2=spectral_obj.curl_hat(aux,0);
	spectral_obj.fft_ssc.inverse_transform(curl_mag_,aux2);
  //curl_curl_mag
	aux=0;
	aux=spectral_obj.curl_hat(aux2,1);
	spectral_obj.fft_ccs.inverse_transform(curl_curl_mag_,aux);
  //grad_temp
	CSF saux(n1,n2/2+1,n3);
	saux=0;
	aux=0;
	spectral_obj.sfft_s.direct_transform(saux,temp_);
	aux=spectral_obj.grad_hat(saux,0);
	spectral_obj.fft_ssc.inverse_transform(grad_temp_,aux);
}

void Basic::load(const string & filename)
{
	Spectral so(32,32,16,2*M_PI,2*M_PI,M_PI);
	RVF v32(32,32,16);
	vtkFileLoad(filename+"_basic_vel",v32);
	CVF v32_hat(32,32/2+1,16);
	v32_hat=0;
	so.fft_ccs.direct_transform(v32_hat,v32);
	so.pnvh_hat(v32_hat);
	CVF v64_hat(n1,n2/2+1,n3);
	v64_hat=0;
	for(int i=0;i<32/2+1;++i)
		for(int j=0;j<32/2+1;++j)
			for(int k=0;k<16;++k)
				v64_hat(i,j,k)=v32_hat(i,j,k);
	for(int i=32/2+1;i<32;++i)
		for(int j=0;j<32/2+1;++j)
			for(int k=0;k<16;++k)
				v64_hat(i+n1-32,j,k)=v32_hat(i,j,k);
	spectral_obj.pnvh_hat(v64_hat);
	spectral_obj.fft_ccs.inverse_transform(this->vel(),v64_hat);
	v32=0;
	vtkFileLoad(filename+"_basic_mag",v32);
	v32_hat=0;
	so.fft_ccs.direct_transform(v32_hat,v32);
	v64_hat=0;
	for(int i=0;i<32/2+1;++i)
		for(int j=0;j<32/2+1;++j)
			for(int k=0;k<16;++k)
				v64_hat(i,j,k)=v32_hat(i,j,k);
	for(int i=32/2+1;i<32;++i)
		for(int j=0;j<32/2+1;++j)
			for(int k=0;k<16;++k)
				v64_hat(i+n1-32,j,k)=v32_hat(i,j,k);
	spectral_obj.pnvh_hat(v64_hat);
	spectral_obj.fft_ccs.inverse_transform(this->mag(),v64_hat);
	RSF s32(32,32,16);
	vtkFileLoad(filename+"_basic_temp",s32);
	CSF s32_hat(32,32/2+1,16);
	s32_hat=0;
	so.sfft_s.direct_transform(s32_hat,s32);
	CSF s64_hat(n1,n2/2+1,n3);
	s64_hat=0;
	for(int i=0;i<32/2+1;++i)
		for(int j=0;j<32/2+1;++j)
			for(int k=0;k<16;++k)
				s64_hat(i,j,k)=s32_hat(i,j,k);
	for(int i=32/2+1;i<32;++i)
		for(int j=0;j<32/2+1;++j)
			for(int k=0;k<16;++k)
				s64_hat(i+n1-32,j,k)=s32_hat(i,j,k);
	spectral_obj.pnvh_hat(s64_hat);
	spectral_obj.sfft_s.inverse_transform(this->temp(),s64_hat);
	eval_derivatives();
}

void Basic::save(const string & filename)
{

	const int & n1 = input_obj.n1;
	const int & n2 = input_obj.n2;
	const int & n3 = input_obj.n3;
	cat::tvector<int,3> dims(n1,n2,n3);
	cat::tvector<double,3> ori(0,0,0);
	cat::tvector<double,3> sp(input_obj.l1/input_obj.n1,input_obj.l2/input_obj.n2,input_obj.l3/input_obj.n3);
	vtkFileSave(filename+"_basic_vel",this->vel(),ori,sp,"VelocityField");
	vtkFileSave(filename+"_basic_mag",this->mag(),ori,sp,"MagneticField");
	vtkFileSave(filename+"_basic_temp",this->temp(),ori,sp,"TemperatureField");
};

void Basic::save_energ_spec(const string & filename)
{
	CBVF aux(n1,n2/2+1,n3);
	spectral_obj.fft_ccs.direct_transform(aux.vel(),this->vel_);
	spectral_obj.fft_ccs.direct_transform(aux.mag(),this->mag_);
	spectral_obj.sfft_s.direct_transform(aux.temp(),this->temp_);
	//Evaluate and save energy spectra
	cout << "Energy spectrum (K Velocity Magnetic Temperature)" << endl;
	cout << spectral_obj.eval_energ_spec(aux,0) << endl;
	ofstream ofs(filename.c_str());
	ofs << spectral_obj.eval_energ_spec(aux,0) << endl;
	ofs.close();
}
