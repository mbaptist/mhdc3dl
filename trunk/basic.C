
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
Basic::Basic(const input & input_obj__,spectral & spectral_obj__):
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
      load(input_obj.basic_vel_fname,
	   input_obj.basic_mag_fname,
	   input_obj.basic_temp_fname);
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
	    gen_random_obj->gen_random_field(vel_,-1.,1.,input_obj.br_ki,input_obj.br_kf,input_obj.br_alpha,input_obj.br_rms_norm,input_obj.br_kind,input_obj.br_sym);
	    gen_random_obj->gen_random_field(mag_,-1.,1.,input_obj.br_ki,input_obj.br_kf,input_obj.br_alpha,input_obj.br_rms_norm,input_obj.br_kind,input_obj.br_sym);
	    gen_random_obj->gen_random_field(temp_,-1.,1.,input_obj.br_ki,input_obj.br_kf,input_obj.br_alpha,input_obj.br_rms_norm,input_obj.br_kind,input_obj.br_sym);
	    
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
      spectral_obj.pnvh(vel_);
      
      cout << "Energy spectrum of Basic velocity: " << endl;
      
      cout << spectral_obj.eval_energ_spec(vv_hat) << endl;
      
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


