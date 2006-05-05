
#include "mhdc3dl_sss.h"

#include <iostream>
#include <fstream>

#include <cat.h>

#include "gen_random.h"

#include "spectral.h"

using namespace std;
using namespace cat;

int main()
{
  input_obj=new input("default.cfg");

  int n1=input_obj->n1;
  int n2=input_obj->n2;
  int n3=input_obj->n3;

  int m=2*7*n1*(n2/2+1)*n3;
  spectral_obj= new spectral(input_obj->n1,
			     input_obj->n2,
			     input_obj->n3,
			     input_obj->l1,
			     input_obj->l2,
			     input_obj->l3);

  double xp=0;
  double eim=0;
  double ep=input_obj->ep;
  double thr=input_obj->thr;
  int mp=input_obj->mp;
  double sc=input_obj->sc;
  int nseq=input_obj->nseq;

  stringstream ofname;

  ofname << input_obj->sss_int_ofbname;

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

  if(input_obj->sss_ifname!="")
    {
      cout << "Loading initial condition from " 
	   << input_obj->sss_ifname << endl;
      vzdeigen_load_ffile_(v1,&m,(input_obj->sss_ifname).c_str());
      basic_start_seed=input_obj->sss_basic_restore_seed;
      basic_stop_seed=input_obj->sss_basic_restore_seed;
      br_seed_step=1;
    }
  else
    {
      int sym=input_obj->sym_sub;
      gen_random_v(v1,sym,input_obj->sss_seed);
      cout << "Initial flow with seed=" << input_obj->sss_seed
	   << " and symmetry: " << sym << endl;
      //basic_start_seed=input_obj->basic_start_seed;
      //basic_stop_seed=input_obj->basic_stop_seed;
      //br_seed_step=input_obj->br_seed_step;
    }

  for(int seed=basic_start_seed;
      seed<=basic_stop_seed;
      seed+=br_seed_step)
    {
      
      basic=new basic_fields(*input_obj,*spectral_obj);//,seed);
      cout << "Basic flow with seed=" << seed << endl;

      a_nought_obj=new a_nought(*input_obj,*spectral_obj,*basic);

      a_nought_obj->sym_sub()=input_obj->sym_sub;
      cout << "Symmetry in linear operator: " << a_nought_obj->sym_sub() << endl;

      vzdeigen_(v1,&xp,&eim,&ep,&thr,&m,
		v2,v3,v4,v5,v6,v7,&mp,&sc,
		&nseq,ofname.str().c_str());
      cout << "xp=" << xp << "  eim=" << eim << endl;
      
      delete a_nought_obj;
      delete basic;

    }

  delete[] v1;
  delete[] v2;
  delete[] v3;
  delete[] v4;
  delete[] v5;
  delete[] v6;
  delete[] v7;

  delete spectral_obj; 
  delete input_obj;
     
  return 0;

}


void gen_random_v(double * v,const int & sym,const int & seed)
{  

  int n1=input_obj->n1;
  int n2=input_obj->n2;
  int n3=input_obj->n3;

  CBVF * bv=
    new CBVF(n1,n2/2+1,n3);

  gen_random b(*input_obj,*spectral_obj,seed);

  b.gen_random_field_hat(bv->vel(),1,10,4,1,0,sym);
  b.gen_random_field_hat(bv->mag(),1,10,4,1,0,sym);
  b.gen_random_field_hat(bv->temp(),1,10,4,1,0,sym);

  
  bv->mag()=1+I;

  bv->vel()=0;
  bv->temp()=0;

  bv2v(v,bv);

  delete bv;

}

void external_prodx_(int * nit,double * vi, double * vo, int * m)
{
  int n1=input_obj->n1;
  int n2=input_obj->n2;
  int n3=input_obj->n3;

  //Allocate block vectors
  CBVF * bv_in=
    new CBVF(n1,n2/2+1,n3);
  CBVF * bv_out=
    new CBVF(n1,n2/2+1,n3);
	

  //cout << linops_obj->sym_sub() << endl;

  //Copy fortran input array into block vector
  v2bv(bv_in,vi); 

//   cout << "basic vel: " << endl;  
//   spectral_obj->pnvh((basic->vel()));
//   cout << "basic mag: " << endl;  
//   spectral_obj->pnvh((basic->mag()));
//   cout << "basic temp: " << endl;  
//   spectral_obj->pnvh((basic->temp()));

  //cout << "in: " << endl;
  //spectral_obj->pnvh(*bv_in);

  //Print energy spectrum
  //cout << spectral_obj->eval_energ_spec((*bv_in).mag(),1,20) << endl;

  //Apply linear operator
  //(*bv_out)=linops_obj->eval_a_nought((*bv_in));

  //Just laplacian
  (*bv_out)=spectral_obj->lap_hat((*bv_in));
  (*bv_out).vel()=0;
  (*bv_out).temp()=0;

  //cout << sum(norm(bv_in->mag())) << endl;

  //Remove laplacian
  //(*bv_out)-=spectral_obj->lap_hat((*bv_in));

  //cout << "out: " << endl;
  //spectral_obj->pnvh(*bv_out);
  

  //Print energy spectrum
  //cout << spectral_obj->eval_energ_spec((*bv_out).mag(),1,20) << endl;
  

  //exit(0);
  

  //Copy output from block vector to fortran output array
  bv2v(vo,bv_out);

  //Increase the number of iterations
  ++(*nit);
  
  delete bv_in;
  delete bv_out;	
}
  



//Copy output from block vector to fortran output array
void bv2v(double * v,const CBVF * bv)
{
  int n1=input_obj->n1;
  int n2=input_obj->n2;
  int n3=input_obj->n3;
  int pos=0;
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2/2+1;++j)
      for(int k=0;k<n3;++k)
	{
 	  for(int comp=0;comp<3;++comp)
 	    {
 	      v[pos]=((bv->vel())(i,j,k)[comp]).real();
 	      ++pos;
 	      v[pos]=((bv->vel())(i,j,k)[comp]).imag();
 	      ++pos;
 	    }
	  for(int comp=0;comp<3;++comp)
	    {
	      v[pos]=((bv->mag())(i,j,k)[comp]).real();
	      ++pos;
	      v[pos]=((bv->mag())(i,j,k)[comp]).imag();
	      ++pos;
	    }
 	  v[pos]=((bv->temp())(i,j,k)).real();
 	  ++pos;
 	  v[pos]=((bv->temp())(i,j,k)).imag();
 	  ++pos;
	}
}

//Copy fortran input array into block vector
void v2bv(CBVF * bv,const double * v)
{
  int n1=input_obj->n1;
  int n2=input_obj->n2;
  int n3=input_obj->n3;
  int pos=0;
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2/2+1;++j)
      for(int k=0;k<n3;++k)
	{
 	  for(int comp=0;comp<3;++comp)
 	    {
 	      //cout << "Pos:" << pos << endl;
 	      ((bv->vel())(i,j,k)[comp]).real()=v[pos];
 	      ++pos;
 	      ((bv->vel())(i,j,k)[comp]).imag()=v[pos];
 	      ++pos;
 	    }
	  for(int comp=0;comp<3;++comp)
	    {
	      ((bv->mag())(i,j,k)[comp]).real()=v[pos];
	      ++pos;
	      ((bv->mag())(i,j,k)[comp]).imag()=v[pos];
	      ++pos;
 	    }
 	  ((bv->temp())(i,j,k)).real()=v[pos];
 	  ++pos;
 	  ((bv->temp())(i,j,k)).imag()=v[pos];
 	  ++pos;
	}
}
