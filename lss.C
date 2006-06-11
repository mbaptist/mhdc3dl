
//class lss - imlementation
//LARGE SCALE STABILITY

#include "input.h"

#include "block_vector.h"
#include "spectral.h"
#include "basic.h"
#include "linops.h"

#include "lss.h"

#include <cat.h>
#include <lass.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;
//using namespace lass;
using namespace cat;

//Ctor
lss::lss(input & input_obj_):
input_obj(input_obj_),
spectral_obj(input_obj.n1,input_obj.n2,input_obj.n3,
             input_obj.l1,input_obj.l2,input_obj.l3),
basic(input_obj,spectral_obj),
a_nought_obj(input_obj,spectral_obj,basic),
a_nought_adjoint_obj(input_obj,spectral_obj,basic),
a_one_obj(input_obj,spectral_obj,basic),
precond_obj(spectral_obj,input_obj.qq),
precond_adjoint_obj(spectral_obj,input_obj.qq_adj)
{
}

//Dtor() 
lss::~lss()
{
}

//Evaluate eddy viscosity
void lss::run(double & lambda_minimal,
              double & lambda_maximal)
{
  //sizes
  int & n1=input_obj.n1;
	int & n2=input_obj.n2;
	int & n3=input_obj.n3;

	double & l1=input_obj.l1;
	double & l2=input_obj.l2;
	double & l3=input_obj.l3;

	cout << .5*sum(dot_product(basic.vel(),basic.vel()))*(l1*l2*l3)/(n1*n2*(n3-1)) << endl;
	CVF vttt(n1,n2/2+1,n3);
	spectral_obj.fft_ccs.direct_transform(vttt,basic.vel());
	//cout << .5*(l1*l2*l3)*spectral_obj.scalar_prod(vttt,vttt) << endl;
	cat::array<double,1> ves(spectral_obj.eval_energ_spec(vttt,0));
	cout << (l1*l2*l3)*sum(ves)*sqrt(max(spectral_obj.wv2))/(ves.size()-1) << endl;
	
	cout << .5*sum(dot_product(basic.mag(),basic.mag()))*(l1*l2*l3)/(n1*n2*(n3-1)) << endl;
	CVF httt(n1,n2/2+1,n3);
	spectral_obj.fft_ccs.direct_transform(httt,basic.mag());
	cat::array<double,1> hes(spectral_obj.eval_energ_spec(httt,0));
	cout << (l1*l2*l3)*sum(hes)*sqrt(max(spectral_obj.wv2))/(hes.size()-1) << endl;

 	double tmp=0;
	for(int i=0;i<n1;++i)
		for(int j=0;j<n2;++j)
			for(int k=0;k<n3;++k)
			{
				if(k==0||k==n3-1)
					tmp+=basic.temp()(i,j,k)*basic.temp()(i,j,k);
				else
					tmp+=2*basic.temp()(i,j,k)*basic.temp()(i,j,k);
			}
	cout << .5*tmp*(l1*l2*l3)/(n1*n2*(n3-1))/2. << endl;
	CSF tttt(n1,n2/2+1,n3);
	spectral_obj.sfft_c.direct_transform(tttt,basic.temp());
	cout << .5*(l1*l2*l3)*spectral_obj.scalar_prod(tttt,tttt,1) << endl;
	cat::array<double,1> tes(spectral_obj.eval_energ_spec(tttt,1));
	cout << (l1*l2*l3)*sum(tes)*sqrt(max(spectral_obj.wv2))/(tes.size()-1) << endl;
	
  //save basic fields
	basic.save(input_obj.basic_vel_fname,input_obj.basic_mag_fname,input_obj.basic_temp_fname);
	
  //define a block vector to contain the constants for
  //the rhs of auxiliary problems (AP)
	CBVF constant(n1,n2/2+1,n3);
	
  //Solving for auxiliary problems
  
  //Order zero
  //solutions of order 0 APs
  CBVF s[4]=
	{CBVF(n1,n2/2+1,n3),
			CBVF(n1,n2/2+1,n3),
			CBVF(n1,n2/2+1,n3),
			CBVF(n1,n2/2+1,n3)};
  //pressure part of order 0 APs
	CSF s_p[4]=
	{CSF(n1,n2/2+1,n3),
			CSF(n1,n2/2+1,n3),
			CSF(n1,n2/2+1,n3),
			CSF(n1,n2/2+1,n3)};
  //s0
	constant=0;
	constant.vel()(0,0,0)=cat::tvector<double,3>(1,0,0);
	solve_zero(s,s_p,constant,0);
  //s1
	constant=0;
	constant.vel()(0,0,0)=cat::tvector<double,3>(0,1,0);
	solve_zero(s,s_p,constant,1);
  //s2
	constant=0;
	constant.mag()(0,0,0)=cat::tvector<double,3>(1,0,0);
	solve_zero(s,s_p,constant,2);
  //s3
	constant=0;
	constant.mag()(0,0,0)=cat::tvector<double,3>(0,1,0);
	solve_zero(s,s_p,constant,3);
	
  //Order one
 
  //solutions of order 1 APs
	CBVF g[4][2]=
	{
		{CBVF(n1,n2/2+1,n3),
				CBVF(n1,n2/2+1,n3)},
		{CBVF(n1,n2/2+1,n3),
				CBVF(n1,n2/2+1,n3)},
		{CBVF(n1,n2/2+1,n3),
				CBVF(n1,n2/2+1,n3)},
		{CBVF(n1,n2/2+1,n3),
				CBVF(n1,n2/2+1,n3)}
	};
	for(int i=0;i<4;++i)
    {
	    stringstream ss;
	    ss << input_obj.apbfname
		    << "_s_"
		    << i
		    << ".dat";
	    save_block_vector(s[i],ss.str());
	    for(int j=0;j<2;++j)
	    {
		    solve_one(g[i][j],s,s_p,i,j);
		    stringstream sss;
		    sss << input_obj.apbfname
			    << "_gamma_"
			    << i
			    << "_"
			    << j
			    << ".dat";
		    save_block_vector(g[i][j],sss.str());
	    }
    }
  
//Evaluate the averages of B_k_Gamma_ij
  cout << "Evaluating averages of B_k_Gamma_ij... " << endl;
  for(int i=0;i<4;++i)
	  for(int j=0;j<2;++j)
		  for(int k=0;k<2;++k)
		  {
	  //Note that the average of the 3rd component is wrong
	  //but it is not used in further evaluataions
			  
	  //Velocity part
			  av_b_k_gamma_ij_vel[i][j][k]=
				  real((a_one_obj(g[i][j],k)).vel()(0,0,0));

			  //cout << 			  (a_one_obj(g[i][j],k)).vel()(0,0,0)<< endl;
			  
	  //	  cout << "   <d_vel(" << i << "," << j << "," << k << ")>="
	  //     << av_b_k_gamma_ij_vel[i][j][k] << endl;
			  
	  //Magnetic part
			  av_b_k_gamma_ij_mag[i][j][k]=
				  real((a_one_obj(g[i][j],k)).mag()(0,0,0));
			  
			  //cout << 			  (a_one_obj(g[i][j],k)).vel()(0,0,0)<< endl;
			  
	  // cout << "   <d_mag(" << i << "," << j << "," << k << ")>="
	 //     << av_b_k_gamma_ij_mag[i][j][k] << endl;
		  }
	cout << "... done! " << endl;



  //Build and diagonalise matrix
	cout << "Finding minimal and maximal growthrates..." << endl;
  
	cat::tvector<double,2> q(1.,0.);
	cat::array<double,2> ep=eval_ep(q);
  double lambda1,lambda2;
	diag(lambda1,lambda2,ep);
  double lambda_max,lambda_min;
	if (lambda1>lambda2)
	{
		lambda_max=lambda1;
		lambda_min=lambda2;
	}
	else
	{
		lambda_max=lambda2;
		lambda_min=lambda1;
	}
// 	int n_int=10000;
// 	double size_int=2*M_PI/n_int;
	double size_int=input_obj.ls_eps/2;
	int n_int=static_cast<int>(2*M_PI/size_int);
	for(int i=1;i<=n_int;++i)
	{
		double theta=i*size_int;
		q=tvector<double,2>(cos(theta),sin(theta));
		cat::array<double,2> ep=eval_ep(q);
		diag(lambda1,lambda2,ep);
		if(lambda1>lambda_max)
			lambda_max=lambda1;
		if(lambda2>lambda_max)
			lambda_max=lambda2;
		if(lambda1<lambda_min)
			lambda_min=lambda1;
		if(lambda2<lambda_min)
			lambda_min=lambda2;
	}
	cout << "...done!" << endl;

	lambda_maximal=lambda_max;
	lambda_minimal=lambda_min;

}
  



void lss::save_block_vector(CBVF & field,
		       const string & fname)
{
  ofstream ofs(fname.c_str());
  ofs << field << endl;
  ofs.close();
}


//Solving order 0 auxiliary problems
void lss::solve_zero(CBVF * s,
		CSF * s_p,
		const CBVF & rhs_constant_zero,
		const int & i)
{
  //Setting sym to anti-symmetric for zero order APs
  //linops_obj->sym_sub()=0;

  cout << "Solving for s(" << i << ")..." << endl;
  CBVF rhs_zero=-a_nought_obj(rhs_constant_zero);

	// cout << sum(rhs_zero.vel()) << endl;

  spectral_obj.remove_gradient(rhs_zero,0);
  s[i]=0;
  cgsolver(a_nought_obj,a_nought_adjoint_obj,
	   precond_obj,precond_adjoint_obj,
	   s[i],rhs_zero,
	   input_obj.ls_eps,
	   input_obj.kk,
	   input_obj.small,
	   input_obj.small_adj);
  s[i]+=rhs_constant_zero;
  cout << "... done!" << endl;
  cout << "Saving s(" << i << ")..." << endl;
  
  cout << "... done!" << endl;
  cout << "Solving for pressure s_p(" << i << ")..." << endl;
  s_p[i]=spectral_obj.
    poisson_hat(spectral_obj.div_hat( (a_nought_obj(s[i]) ).vel(),0));
  cout << "... done!" << endl;
  cout << "Printing energy spectrum" << endl;
  cout << spectral_obj.eval_energ_spec(s[i].mag(),0) << endl;
  cout << "... done!" << endl;
//   cout << "Printing non-vanishing harmonics in s" << endl;
//   spectral_obj.pnvh(s[i]);
	cout << "... done!" << endl;
}

//Solving order 1 auxiliary problems
//The rhs of the modified first order AP is
//M_ij = a_one_j S_i + gradient of artificial pressures + e_j S_p_i
//and is stored in b
//The solution is gij
void lss::solve_one(CBVF & gamma,
	       CBVF * s,
	       CSF * s_p,
	       const int & i,const int & j)
{  
  //sizes
  int & n1(input_obj.n1);
  int & n2(input_obj.n2);
  int & n3(input_obj.n3);  
 
  //Setting sym to sym for setting up first order APs
  //linops_obj.sym_sub()=1;
  
  cout << "Solving for gamma(" << i << "," << j << ") ..." << endl;

  //Evaluate artificial pressures
  CSF art_press_vel
    (spectral_obj.poisson_hat(-(((s[i]).vel())[j])));
  CSF art_press_mag
    (spectral_obj.poisson_hat(-(((s[i]).mag())[j])));
  //Evaluate the RHS of the modified Auxiliary Problem
  CBVF b(n1,n2/2+1,n3);
  b=0;
  (b.vel())[j]=s_p[i];
  b+=-a_one_obj(s[i],j);

  CBVF grad_art_press(n1,n2/2+1,n3);
  grad_art_press.vel()=spectral_obj.grad_hat(art_press_vel,1);
  grad_art_press.mag()=spectral_obj.grad_hat(art_press_mag,1);
  grad_art_press.temp()=0;

  b+=-a_nought_obj(grad_art_press);

  spectral_obj.remove_gradient(b,0);
	
  //Solve the modified problem
  gamma=0;
  cgsolver(a_nought_obj,a_nought_adjoint_obj,
	   precond_obj,precond_adjoint_obj,
	   gamma,b,
	   input_obj.ls_eps,
	   input_obj.kk,
	   input_obj.small,
	   input_obj.small_adj);
	 
  //Obtain the solution to the true auxiliary problem 
  //by adding the gradients
  gamma+=grad_art_press;
  
  cout << "Energy spectrum" << endl;
  cout << spectral_obj.eval_energ_spec(gamma.mag(),0) << endl; 

//   cout << "Printing non-vanishing harmonics in gamma" << endl;
//   spectral_obj.pnvh(gamma);

  cout << "... done!" << endl;
 
}

//Evaluate ep (2x2 matrix to be diagonalised)
cat::array<double,2> lss::eval_ep(const cat::tvector<double,2> & q)
{
  cat::array<double,2> out(2,2);
  cat::array<double,2> e=eval_e(q);
  out(0,0)=e(0,0)*q[1]*q[1]-e(0,1)*q[0]*q[1]-e(1,0)*q[1]*q[0]+e(1,1)*q[0]*q[0];
  out(0,1)=e(0,2)*q[1]*q[1]-e(0,3)*q[0]*q[1]-e(1,2)*q[1]*q[0]+e(1,3)*q[0]*q[0];
  out(1,0)=e(2,0)*q[1]*q[1]-e(2,1)*q[0]*q[1]-e(3,0)*q[1]*q[0]+e(3,1)*q[0]*q[0];
  out(1,1)=e(2,2)*q[1]*q[1]-e(2,3)*q[0]*q[1]-e(3,2)*q[1]*q[0]+e(3,3)*q[0]*q[0];
  return out;
}

//Evaluate e (4x4 matrix to be reduced to 2x2)
cat::array<double,2> lss::eval_e(const cat::tvector<double,2> & q)
{
	cat::array<double,2> out(4,4);
	for(int i=0;i<4;++i)
	{
		out(0,i)=0;
		out(1,i)=0;
		out(2,i)=0;
		out(3,i)=0;
		for(int k=0;k<2;++k)
			for(int j=0;j<2;++j)
			{
				out(0,i)+=(av_b_k_gamma_ij_vel[i][j][k])[0]*q[j]*q[k];
				out(1,i)+=(av_b_k_gamma_ij_vel[i][j][k])[1]*q[j]*q[k];
				out(2,i)+=(av_b_k_gamma_ij_mag[i][j][k])[0]*q[j]*q[k];
				out(3,i)+=(av_b_k_gamma_ij_mag[i][j][k])[1]*q[j]*q[k];
			}
		out(0,i)+=(input_obj.visc)*(0==i);
		out(1,i)+=(input_obj.visc)*(1==i);
		out(2,i)+=(input_obj.diff)*(2==i);
		out(3,i)+=(input_obj.diff)*(3==i);
	}
	out*=-1;
  return out;
}

void lss::diag(double & lambda1,double & lambda2,const cat::array<double,2> & matrix)
{
  double a=1.;
  double b=-matrix(0,0)-matrix(1,1);
  double c=matrix(0,0)*matrix(1,1)-matrix(0,1)*matrix(1,0);
  lambda1=-b/(2.*a)*(1.+sqrt(1.-4.*a*c/(b*b)));
  lambda2=c/(a*lambda1);
}

