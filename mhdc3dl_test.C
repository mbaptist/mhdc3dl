
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


input * input_obj;

Spectral * spectral_obj;

int n1,n2,n3;

double l1,l2,l3;

void test_adjoint();

void test_cross_product();

void test_remove_gradient();


int main()
{

	input_obj= new input("mhdc3dl_run_default.py");

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

  a_nought a_nought_obj(*input_obj,*spectral_obj,basic);
  a_nought_adjoint a_nought_adj_obj(*input_obj,*spectral_obj,basic);

  CBVF a(n1,n2/2+1,n3),b(a),c(a),d(a);

  gen_random gr(*input_obj,*spectral_obj,555);

	const int skmin=0;
	const int skmax=10;
	const double decay=3;
	const double nnn=1;
	const bool kkk=0;
	const bool sss=1;
	
  gr.gen_random_field_hat(a.vel(),skmin,skmax,decay,nnn,kkk,sss);
	gr.gen_random_field_hat(a.mag(),skmin,skmax,decay,nnn,kkk,sss);
	gr.gen_random_field_hat(a.temp(),skmin,skmax,decay,nnn,kkk,sss);
  
	gr.gen_random_field_hat(b.vel(),skmin,skmax,decay,nnn,kkk,sss);
	gr.gen_random_field_hat(b.mag(),skmin,skmax,decay,nnn,kkk,sss);
	gr.gen_random_field_hat(b.temp(),skmin,skmax,decay,nnn,kkk,sss);


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

   c=a_nought_obj(b);
   d=a_nought_adj_obj(a);


  //   cout << "c=Ab" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*a" << endl;
  //   spectral_obj.pnvh(d);

  double sp1,sp2;

  sp1=spectral_obj->scalar_prod(a,c,kkk);
  sp2=spectral_obj->scalar_prod(d,b,kkk);

  cout << "<a,a_nought(b)>= " << sp1 << endl;
  cout << "<a_nought_star(a),b>= " << sp2 << endl;
  cout << "||<a,a_nought(b)>-<a_nought_star(a),b>||= " << abs(sp1-sp2) << endl;

  cout << endl;

	c=0;
	d=0;
  
  c=a_nought_obj(a);
  d=a_nought_adj_obj(b);


  //   cout << "c=Aa" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*b" << endl;
  //   spectral_obj.pnvh(d);

  sp1=spectral_obj->scalar_prod(b,c,kkk);
  sp2=spectral_obj->scalar_prod(d,a,kkk);

  cout << "<b,a_nought(a)>= " << sp1 << endl;
  cout << "<a_nought_star(b),a>= " << sp2 << endl;
  cout << "||<b,a_nought(a)>-<a_nought_star(b),a>||= " << abs(sp1-sp2) << endl;

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

  c=eval_a_nought(b);
  d=eval_a_nought_star(a);


  //   cout << "c=Ab" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*a" << endl;
  //   spectral_obj.pnvh(d);

  Real sp1,sp2;

  sp1=spectral_obj.scalar_prod(a,c);
  sp2=spectral_obj.scalar_prod(d,b);

  cout << "<a,a_nought(b)>= " << sp1 << endl;
  cout << "<a_nought_star(a),b>= " << sp2 << endl;
  cout << "||<a,a_nought(b)>-<a_nought_star(a),b>||= " << abs(sp1-sp2) << endl;

  cout << endl;


  c=eval_a_nought(a);
  d=eval_a_nought_star(b);


  //   cout << "c=Aa" << endl;
  //   spectral_obj.pnvh(c);
  //   cout << "d=A*b" << endl;
  //   spectral_obj.pnvh(d);

  sp1=spectral_obj.scalar_prod(b,c);
  sp2=spectral_obj.scalar_prod(d,a);

  cout << "<b,a_nought(a)>= " << sp1 << endl;
  cout << "<a_nought_star(b),a>= " << sp2 << endl;
  cout << "||<b,a_nought(a)>-<a_nought_star(b),a>||= " << abs(sp1-sp2) << endl;

  cout << endl;

}

#endif


