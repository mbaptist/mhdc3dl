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


#include "input.h"

#include "globals.h"

#include "spectral.h"

#include "gen_random.h"

#include "basic.h"

#include "linops.h"

#include <cat.h>
using namespace cat;

/////////////
// TESTING
////////////


Input * input_obj;

Spectral * spectral_obj;

int n1,n2,n3;

double l1,l2,l3;

void test_adjoint();

void test_cross_product();

void test_remove_gradient();


int main()
{

  input_obj= new Input("mhdc3dl_run_default.py");

  n1=input_obj->n1;
  n2=input_obj->n2;
  n3=input_obj->n3;
  l1=input_obj->l1;
  l2=input_obj->l2;
  l3=input_obj->l3;

  spectral_obj=new Spectral(n1,n2,n3,l1,l2,l3);

  test_adjoint();

  //test_cross_product();

  //test_remove_gradient();

  return 0;

}


void test_adjoint()
{

  Basic basic(*input_obj,*spectral_obj);

  ANought ANought_obj(*input_obj,*spectral_obj,basic);
  ANoughtAdjoint ANought_adj_obj(*input_obj,*spectral_obj,basic);

  CBVF a(n1,n2/2+1,n3),b(a),c(a),d(a);

  gen_random gr(*input_obj,*spectral_obj,555);

  const int skmin=0;
  const int skmax=10;
  const double decay=3;
  const double nnn=1;
  const bool kkk=0;
  const bool sss=1;

  gr.gen_random_field_hat(a.vel(),skmin,skmax,decay,nnn,kkk,sss,"power");
  gr.gen_random_field_hat(a.mag(),skmin,skmax,decay,nnn,kkk,sss,"power");
  gr.gen_random_field_hat(a.temp(),skmin,skmax,decay,nnn,kkk,sss,"power");

  gr.gen_random_field_hat(b.vel(),skmin,skmax,decay,nnn,kkk,sss,"power");
  gr.gen_random_field_hat(b.mag(),skmin,skmax,decay,nnn,kkk,sss,"power");
  gr.gen_random_field_hat(b.temp(),skmin,skmax,decay,nnn,kkk,sss,"power");


  //   a=0;

  //   a.vel()(1,1,1)=cat::Tvector<complex<double>,3>(1.,2.,3.);
  //   a.mag()(1,1,1)=cat::Tvector<complex<double>,3>(2.,2.,3.);
  //   a.temp()(1,1,1)=4;

  // 	a.vel()=0;
  //
  // 	a.temp()=0;

  //   b=0;

  //   b.vel()(1,1,1)=cat::Tvector<complex<double>,3>(1.,1.,8.);
  //   b.mag()(1,1,1)=cat::Tvector<complex<double>,3>(2.,3.,3.);
  //   b.temp()(1,1,1)=5;

  // 	b.vel()=0;
  //
  // 	b.temp()=0;

  c=ANought_obj(b);
  d=ANought_adj_obj(a);


  //   cout << "c=Ab" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*a" << endl;
  //   spectral_obj.pnvh(d);

  double sp1,sp2;

  sp1=spectral_obj->scalar_prod(a,c,kkk);
  sp2=spectral_obj->scalar_prod(d,b,kkk);

  cout << "<a,ANought(b)>= " << sp1 << endl;
  cout << "<ANought_star(a),b>= " << sp2 << endl;
  cout << "||<a,ANought(b)>-<ANought_star(a),b>||= " << abs(sp1-sp2) << endl;

  cout << endl;

  c=0;
  d=0;

  c=ANought_obj(a);
  d=ANought_adj_obj(b);


  //   cout << "c=Aa" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*b" << endl;
  //   spectral_obj.pnvh(d);

  sp1=spectral_obj->scalar_prod(b,c,kkk);
  sp2=spectral_obj->scalar_prod(d,a,kkk);

  cout << "<b,ANought(a)>= " << sp1 << endl;
  cout << "<ANought_star(b),a>= " << sp2 << endl;
  cout << "||<b,ANought(a)>-<ANought_star(b),a>||= " << abs(sp1-sp2) << endl;

  cout << endl;

}

#if 0
void test_cross_product()
{
  RVF a(n1,n2,n3);
  RVF b(n1,n2,n3);
  RVF cp(n1,n2,n3);
  CVF a_hat(n1,n2/2+1,n3);
  CVF b_hat(n1,n2/2+1,n3);
  CVF cp_hat(n1,n2/2+1,n3);
  a_hat=0;
  b_hat=0;

  a_hat(0,1,0)=cat::Tvector<complex<double>,3>(1,0,0);
  b_hat(0,1,0)=cat::Tvector<complex<double>,3>(0,1,0);

  spectral_obj->fft_ccs.inverse_transform(a,a_hat);
  spectral_obj->fft_ccs.inverse_transform(b,b_hat);
  cp=cross_product(a,b);
  spectral_obj->fft_ssc.direct_transform(cp_hat,cp);

  cout << "pnvh a: " << endl;
  spectral_obj->pnvh(a);
  cout << "pnvh b: " << endl;
  spectral_obj->pnvh(b);
  cout << "pnvh cp: " << endl;
  spectral_obj->pnvh(cp_hat);


#if 0

  a_hat=0;
  a_hat(0,0,0)=1;
  spectral_obj->fft_ssc.inverse_transform(a,a_hat);
  a_hat=0;
  spectral_obj->fft_ssc.direct_transform(a_hat,a);
  spectral_obj->pnvh(a_hat);
#endif

#if 0

  for(int i=0;i<n1;++i)
    for(int j=0;j<n2;++j)
      for(int k=0;k<n3;++k)

#endif

      }
#endif



void test_remove_gradient()
{
  RVF v(n1,n2,n3);
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2;++j)
      for(int k=0;k<n3;++k)
        v(i,j,k)=cat::Tvector<double,3>
                 (cos(l1/n1*i)*cos(l3/(n3-1)*k),
                  cos(l2/n2*j)*cos(l3/(n3-1)*k),
                  sin(l1/n1*i+l2/n2*j)*sin(l3/(n3-1)*k));

  CVF v_hat(n1,n2/2+1,n3);
  spectral_obj->fft_ccs.direct_transform(v_hat,v);

  cout << "v_hat before: " << endl;
  spectral_obj->pnvh_hat(v_hat);

  cout << "removing gradient... " << endl;
  CVF gv_hat(n1,n2/2+1,n3);
  gv_hat=spectral_obj->remove_gradient(v_hat,0);
  cout << "...done!" << endl;

  cout << "v_hat after: " << endl;
  spectral_obj->pnvh_hat(v_hat);

  cout << "gradient: " << endl;
  spectral_obj->pnvh_hat(gv_hat);



}




#if 0

void linops::test_adjoint(CBVF & a,CBVF & b)
{

  //   cout << "a" << endl;
  //   spectral_obj.pnvh(a);
  //   cout << "b" << endl;
  //   spectral_obj.pnvh(b);

  CBVF c(a),d(a);

  c=eval_ANought(b);
  d=eval_ANought_star(a);


  //   cout << "c=Ab" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*a" << endl;
  //   spectral_obj.pnvh(d);

  Real sp1,sp2;

  sp1=spectral_obj.scalar_prod(a,c);
  sp2=spectral_obj.scalar_prod(d,b);

  cout << "<a,ANought(b)>= " << sp1 << endl;
  cout << "<ANought_star(a),b>= " << sp2 << endl;
  cout << "||<a,ANought(b)>-<ANought_star(a),b>||= " << abs(sp1-sp2) << endl;

  cout << endl;


  c=eval_ANought(a);
  d=eval_ANought_star(b);


  //   cout << "c=Aa" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*b" << endl;
  //   spectral_obj.pnvh(d);

  sp1=spectral_obj.scalar_prod(b,c);
  sp2=spectral_obj.scalar_prod(d,a);

  cout << "<b,ANought(a)>= " << sp1 << endl;
  cout << "<ANought_star(b),a>= " << sp2 << endl;
  cout << "||<b,ANought(a)>-<ANought_star(b),a>||= " << abs(sp1-sp2) << endl;

  cout << endl;

}

#endif

#if 0
cout << "Basic Velocity:" << endl;
cout << "energy in real space: " << spectral_obj.energy(basic.vel()) << endl;
CVF vttt(n1,n2/2+1,n3);
spectral_obj.fft_ccs.direct_transform(vttt,basic.vel());
cout <<  "energy in fourier space (using the scalar product): " << .5*(l1*l2*l3)*spectral_obj.scalar_prod(vttt,vttt,0) << endl;
cat::Array<double,1> ves(spectral_obj.eval_energ_spec(vttt,0));
cout << "energy in fourier space (as sum of the enery spectrum): " << (l1*l2*l3)*sum(ves)*spectral_obj.wnstep << endl;

cout << "Basic Magnetic Field:" << endl;
cout << "energy in real space: " << spectral_obj.energy(basic.mag()) << endl;
CVF httt(n1,n2/2+1,n3);
spectral_obj.fft_ccs.direct_transform(httt,basic.mag());
cout <<  "energy in fourier space (using the scalar product): " << .5*(l1*l2*l3)*spectral_obj.scalar_prod(httt,httt,0) << endl;
cat::Array<double,1> hes(spectral_obj.eval_energ_spec(httt,0));
cout << "energy in fourier space (as sum of the enery spectrum): " << (l1*l2*l3)*sum(hes)*spectral_obj.wnstep << endl;

cout << "Basic Temperature:" << endl;
cout << "energy in real space: " << spectral_obj.energy(basic.temp()) << endl;
CSF tttt(n1,n2/2+1,n3);
spectral_obj.sfft_s.direct_transform(tttt,basic.temp());
cout <<  "energy in fourier space (using the scalar product): " << .5*(l1*l2*l3)*spectral_obj.scalar_prod(tttt,tttt,0) << endl;
cat::Array<double,1> tes(spectral_obj.eval_energ_spec(tttt,0));
cout << "energy in fourier space (as sum of the enery spectrum): " << (l1*l2*l3)*sum(tes)*spectral_obj.wnstep << endl;
#endif
