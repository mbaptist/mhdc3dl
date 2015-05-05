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
sss::sss(Input & input_obj_):
    input_obj(input_obj_),
    spectral_obj(input_obj.n1,input_obj.n2,input_obj.n3,
                 input_obj.l1,input_obj.l2,input_obj.l3),
    basic(input_obj,spectral_obj),
    ANought_obj(input_obj,spectral_obj,basic)
{
  ExternalProdx::instance().switch_ptr(this);
}

//Dtor()
sss::~sss()
{}

void sss::run(double & xp,double & eim)
{
  int n1=input_obj.n1;
  int n2=input_obj.n2;
  int n3=input_obj.n3;

  int m=2*7*n1*(n2/2+1)*n3;

  //double xp=0;
  //double eim=0;
  xp=0;
  eim=0;
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
	
	//Create complex block vectors from pointers to double
	CBVF bv(cat::Tvector<int,3>(n1,n2/2+1,n3),v);
	
  gen_random b(input_obj,spectral_obj,seed);

  b.gen_random_field_hat(bv.vel(),0,5,4,1,0,sym,"power");
  b.gen_random_field_hat(bv.mag(),0,5,4,1,0,sym,"power");
  b.gen_random_field_hat(bv.temp(),0,5,4,1,0,sym,"power");
	
}

void sss::external_prodx(int * nit,double * vi, double * vo, int * m)
{
  int n1=input_obj.n1;
  int n2=input_obj.n2;
  int n3=input_obj.n3;

	//Create complex block vectors from pointers to double
	//Data is not copied
	CBVF bv_in(cat::Tvector<int,3>(n1,n2/2+1,n3),vi);
	CBVF bv_out(cat::Tvector<int,3>(n1,n2/2+1,n3),vo);
	
  //Apply linear operator
  bv_out=ANought_obj(bv_in);

  //Increase the number of iterations
  ++(*nit);

}

