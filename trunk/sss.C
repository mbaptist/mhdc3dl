
#include "sss.h"

#include "globals.h"

#include <iostream>
#include <fstream>

#include <cat.h>

#include "gen_random.h"

#include "spectral.h"

using namespace std;
using namespace cat;

void external_prodx_(int * nit,double * vi, double * vo, int * m)
{
	return ExternalProdx::instance()(nit,vi,vo,m);
}


///////////////////////
// Class sss
//////////////////////

//Ctor
sss::sss(input & input_obj_):
input_obj(input_obj_),
spectral_obj(input_obj.n1,input_obj.n2,input_obj.n3,
             input_obj.l1,input_obj.l2,input_obj.l3),
basic(input_obj,spectral_obj),
a_nought_obj(input_obj,spectral_obj,basic)
{
	ExternalProdx::instance().switch_ptr(this);
}

//Dtor() 
sss::~sss()
{
}

void sss::run()
{

	
  int n1=input_obj.n1;
  int n2=input_obj.n2;
  int n3=input_obj.n3;

  int m=2*7*n1*(n2/2+1)*n3;

  double xp=0;
  double eim=0;
  double ep=input_obj.ep;
  double thr=input_obj.thr;
  int mp=input_obj.mp;
  double sc=input_obj.sc;
  int nseq=input_obj.nseq;

  stringstream ofname;

  ofname << input_obj.sss_int_ofbname;

  double * v1 = new double[m];
  double * v2 = new double[m];
  double * v3 = new double[m];
  double * v4 = new double[m];
  double * v5 = new double[m];
  double * v6 = new double[m];
  double * v7 = new double[m];
  
  int basic_start_seed;
  int basic_stop_seed;
  int br_seed_step;

  if(input_obj.sss_ifname!="")
    {
      cout << "Loading initial condition from " 
	   << input_obj.sss_ifname << endl;
      vzdeigen_load_ffile_(v1,&m,(input_obj.sss_ifname).c_str());
      basic_start_seed=input_obj.sss_basic_restore_seed;
      basic_stop_seed=input_obj.sss_basic_restore_seed;
      br_seed_step=1;
    }
  else
    {
      int sym=input_obj.sym_sub;
      gen_random_v(v1,sym,input_obj.sss_seed);
      cout << "Initial flow with seed=" << input_obj.sss_seed
	   << " and symmetry: " << sym << endl;

    }


      vzdeigen_(v1,&xp,&eim,&ep,&thr,&m,
		v2,v3,v4,v5,v6,v7,&mp,&sc,
		&nseq,ofname.str().c_str());
      cout << "xp=" << xp << "  eim=" << eim << endl;
      

  delete[] v1;
  delete[] v2;
  delete[] v3;
  delete[] v4;
  delete[] v5;
  delete[] v6;
  delete[] v7;

}


void sss::gen_random_v(double * v,const int & sym,const int & seed)
{  

  int n1=input_obj.n1;
  int n2=input_obj.n2;
  int n3=input_obj.n3;

  CBVF bv(n1,n2/2+1,n3);

  gen_random b(input_obj,spectral_obj,seed);

  b.gen_random_field_hat(bv.vel(),0,5,4,1,0,sym);
  b.gen_random_field_hat(bv.mag(),0,5,4,1,0,sym);
  b.gen_random_field_hat(bv.temp(),0,5,4,1,0,sym);

  bv2v(v,bv);

}

void sss::external_prodx(int * nit,double * vi, double * vo, int * m)
{
  int n1=input_obj.n1;
  int n2=input_obj.n2;
  int n3=input_obj.n3;

  //Allocate block vectors
  CBVF bv_in(n1,n2/2+1,n3);
  CBVF bv_out(n1,n2/2+1,n3);
	
  //Copy fortran input array into block vector
  v2bv(bv_in,vi); 

  //Apply linear operator
	bv_out=a_nought_obj(bv_in);

  //Just laplacian
	//bv_out=spectral_obj.lap_hat(bv_in);
	
  //Remove laplacian
  //bv_out-=spectral_obj.lap_hat(bv_in);

  //Copy output from block vector to fortran output array
  bv2v(vo,bv_out);

  //Increase the number of iterations
  ++(*nit);
  
}
  



//Copy output from block vector to fortran output array
void sss::bv2v(double * v,const CBVF & bv)
{
  int n1=input_obj.n1;
  int n2=input_obj.n2;
  int n3=input_obj.n3;
  int pos=0;
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2/2+1;++j)
      for(int k=0;k<n3;++k)
	{
 	  for(int comp=0;comp<3;++comp)
 	    {
 	      v[pos]=((bv.vel())(i,j,k)[comp]).real();
 	      ++pos;
 	      v[pos]=((bv.vel())(i,j,k)[comp]).imag();
 	      ++pos;
 	    }
	  for(int comp=0;comp<3;++comp)
	    {
	      v[pos]=((bv.mag())(i,j,k)[comp]).real();
	      ++pos;
	      v[pos]=((bv.mag())(i,j,k)[comp]).imag();
	      ++pos;
	    }
 	  v[pos]=((bv.temp())(i,j,k)).real();
 	  ++pos;
 	  v[pos]=((bv.temp())(i,j,k)).imag();
 	  ++pos;
	}
}

//Copy fortran input array into block vector
void sss::v2bv(CBVF & bv,const double * v)
{
  int n1=input_obj.n1;
  int n2=input_obj.n2;
  int n3=input_obj.n3;
  int pos=0;
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2/2+1;++j)
      for(int k=0;k<n3;++k)
	{
 	  for(int comp=0;comp<3;++comp)
 	    {
 	      //cout << "Pos:" << pos << endl;
 	      ((bv.vel())(i,j,k)[comp]).real()=v[pos];
 	      ++pos;
 	      ((bv.vel())(i,j,k)[comp]).imag()=v[pos];
 	      ++pos;
 	    }
	  for(int comp=0;comp<3;++comp)
	    {
	      ((bv.mag())(i,j,k)[comp]).real()=v[pos];
	      ++pos;
	      ((bv.mag())(i,j,k)[comp]).imag()=v[pos];
	      ++pos;
 	    }
 	  ((bv.temp())(i,j,k)).real()=v[pos];
 	  ++pos;
 	  ((bv.temp())(i,j,k)).imag()=v[pos];
 	  ++pos;
	}
}
