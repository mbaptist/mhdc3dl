
//class lss - imlementation
//LARGE SCALE STABILITY

#include "input.h"

#include "block_vector.h"
#include "spectral.h"
#include "basic.h"
#include "linops.h"

#include "lss.h"

#include "io.h"

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
    n1(input_obj.n1),
    n2(input_obj.n2),
    n3(input_obj.n3),
    l1(input_obj.l1),
    l2(input_obj.l2),
    l3(input_obj.l3),
    spectral_obj(n1,n2,n3,l1,l2,l3),
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
void lss::run(double & theta_min,std::complex<double> & lambda_min,double & theta_max,std::complex<double> & lambda_max)
{	
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
    
    //define a block vector to contain the constants for
    //the rhs of order 0 auxiliary problems (AP)
    CBVF constant(n1,n2/2+1,n3);
    
    //Solving for auxiliary problems
    
    //Order zero
    //solutions of order 0 APs
    CBVF s[4]={CBVF(n1,n2/2+1,n3),CBVF(n1,n2/2+1,n3),CBVF(n1,n2/2+1,n3),CBVF(n1,n2/2+1,n3)};
    
    //Initial condition
    for(int i=0;i<4;++i)
    {
	if(input_obj.refine||input_obj.resume)
	{
	    stringstream ss,lr_ss;
	    lr_ss << input_obj.lr_runsname << "_s_" << i;
	    rawload_aux_field_hat(lr_ss.str(),s[i],input_obj.lr_n1,input_obj.lr_n2,input_obj.lr_n3);
	    ss << input_obj.runsname << "_s_" << i << "_read_out";
	    rawsave_aux_field_hat(ss.str(),s[i]);
	    vtksave_aux_field_real(ss.str(),s[i]);
	}
	else
	{
	    s[i]=0;
	}
    }
    
    //pressure part of order 0 APs
    CSF s_p[4]={CSF(n1,n2/2+1,n3),CSF(n1,n2/2+1,n3),CSF(n1,n2/2+1,n3),CSF(n1,n2/2+1,n3)};
    //s0
    constant=0;
    constant.vel()(0,0,0)=cat::Tvector<double,3>(1,0,0);
    solve_zero(s,s_p,constant,0);
    //s1
    constant=0;
    constant.vel()(0,0,0)=cat::Tvector<double,3>(0,1,0);
    solve_zero(s,s_p,constant,1);
    //s2
    constant=0;
    constant.mag()(0,0,0)=cat::Tvector<double,3>(1,0,0);
    solve_zero(s,s_p,constant,2);
    //s3
    constant=0;
    constant.mag()(0,0,0)=cat::Tvector<double,3>(0,1,0);
    solve_zero(s,s_p,constant,3);
    
    //Order one
    
    //solutions of order 1 APs
    CBVF g[4][2]=
    {
	{CBVF(n1,n2/2+1,n3),CBVF(n1,n2/2+1,n3)},
	{CBVF(n1,n2/2+1,n3),CBVF(n1,n2/2+1,n3)},
	{CBVF(n1,n2/2+1,n3),CBVF(n1,n2/2+1,n3)},
	{CBVF(n1,n2/2+1,n3),CBVF(n1,n2/2+1,n3)}
    };
    
    //initial condition
    for(int i=0;i<4;++i)
    {
    for(int j=0;j<2;++j)
    {
	//Initial condition
	if(input_obj.refine||input_obj.resume)
	{
	    stringstream lr_ss;
	    lr_ss << input_obj.lr_runsname << "_gamma_" << i << "_" << j;
	    rawload_aux_field_hat(lr_ss.str(),g[i][j],input_obj.lr_n1,input_obj.lr_n2,input_obj.lr_n3);
	}
	else
	{
	    g[i][j]=0;
	}
	// solve for gamma_i_j
	solve_one(g[i][j],s,s_p,i,j);
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
	    //Magnetic part
	    av_b_k_gamma_ij_mag[i][j][k]=
		real((a_one_obj(g[i][j],k)).mag()(0,0,0));
	    
	    cout << 				av_b_k_gamma_ij_vel[i][j][k] << "\n" << av_b_k_gamma_ij_mag[i][j][k] << endl;
	    
	}
    cout << "... done! " << endl;
    
    
    
    //Build and diagonalise matrix
    cout << "Finding minimal and maximal growthrates..." << endl;
    double theta=0.;
    cat::Tvector<double,2> q(cos(theta),sin(theta));
    cat::Array<double,2> ep=eval_ep(q);
    std::complex<double> lambda1,lambda2;
    diag(lambda1,lambda2,ep);
    //cout << lambda1 << " " << lambda2 << endl;
    //cout << ep << endl;
    theta_min=0.;
    lambda_min=(lambda1.real()<lambda2.real() ? lambda1 : lambda2);
    theta_max=0.;
    lambda_max=(lambda1.real()>lambda2.real() ? lambda1 : lambda2);
    double increment=.5e-6;
    while(theta<=2.*M_PI)
    {
	theta+=increment;
	q=cat::Tvector<double,2>(cos(theta),sin(theta));
	cat::Array<double,2> ep=eval_ep(q);
	diag(lambda1,lambda2,ep);
	//cout << theta << " " << q << " " << lambda1 << " " << lambda2 << endl;
	if(lambda1.real()<lambda_min.real())
	{
	    lambda_min=lambda1;
	    theta_min=theta;
	}
	if(lambda2.real()<lambda_min.real())
	{
	    lambda_min=lambda2;
	    theta_min=theta;
	}		
	if(lambda1.real()>lambda_max.real())
	{
	    lambda_max=lambda1;
	    theta_max=theta;
	}
	if(lambda2.real()>lambda_max.real())
	{
	    lambda_max=lambda2;
	    theta_max=theta;
	}
	
    }
    cout << "...done!" << endl;
}


//Solving order 0 auxiliary problems
void lss::solve_zero(CBVF * s,CSF * s_p,const CBVF & rhs_constant_zero,const int & i)
{
    //Solving for s
    cout << "Solving for s(" << i << ")..." << endl;
    CBVF rhs_zero=-a_nought_obj(rhs_constant_zero);//eval RHS
     spectral_obj.remove_gradient(rhs_zero,0);
    if(input_obj.refine||input_obj.resume)
	s[i]-=rhs_constant_zero;
    cgsolver(a_nought_obj,a_nought_adjoint_obj,
	     precond_obj,precond_adjoint_obj,
	     s[i],rhs_zero,
	     input_obj.ls_eps,
	     input_obj.kk,
	     input_obj.small,
	     input_obj.small_adj);
    s[i]+=rhs_constant_zero;
    cout << "... done!" << endl;
    //Solving for s_p
    cout << "Solving for pressure s_p(" << i << ")..." << endl;
    s_p[i]=spectral_obj.poisson_hat(spectral_obj.div_hat( (a_nought_obj(s[i]) ).vel(),0));
    cout << "... done!" << endl;
    //Save s
    cout << "Saving s(" << i << ")..." << endl;
    stringstream ss;
    ss << input_obj.runsname << "_s_" << i;
    rawsave_aux_field_hat(ss.str(),s[i]);
    vtksave_aux_field_real(ss.str(),s[i]);
    cout << "... done!" << endl;
    //Evaluate and save energy spectra
    cout << "Evaluating and saving energy spectrum for s(" << i << ")..." << endl;
    ss << "_spec.dat";
    ofstream ofs(ss.str().c_str());
    ofs << spectral_obj.eval_energ_spec(s[i],0) << endl;
    ofs.close();
    cout << "... done!" << endl;
    //Print energy spectra
    cout << "Energy spectrum (K Velocity Magnetic Temperature)" << endl;	
    cout << spectral_obj.eval_energ_spec(s[i],0) << endl;
    //   cout << "Printing non-vanishing harmonics in s" << endl;
    //   spectral_obj.pnvh(s[i]);
}

//Solving order 1 auxiliary problems
//The rhs of the modified first order AP is
//M_ij = a_one_j S_i + gradient of artificial pressures + e_j S_p_i
//and is stored in b
//The solution is gij
void lss::solve_one(CBVF & gamma,CBVF * s,CSF * s_p,const int & i,const int & j)
{  
    cout << "Solving for gamma(" << i << "," << j << ") ..." << endl;
    //Evaluate artificial pressures and gradients
    CSF art_press_vel
	(spectral_obj.poisson_hat( CSF(-(((s[i]).vel())[j]))));
    CSF art_press_mag
	(spectral_obj.poisson_hat( CSF(-(((s[i]).mag())[j]))));
    CBVF grad_art_press(n1,n2/2+1,n3);
    grad_art_press.vel()=spectral_obj.grad_hat(art_press_vel,1);
    grad_art_press.mag()=spectral_obj.grad_hat(art_press_mag,1);
    grad_art_press.temp()=0;
    //Evaluate the RHS of the original order 1 problem
    CBVF b(n1,n2/2+1,n3);
    b=0;
    (b.vel())[j]=s_p[i];
    b-=a_one_obj(s[i],j);
    //Evaluate the RHS of the modified Auxiliary Problem
    b-=a_nought_obj(grad_art_press);
    spectral_obj.remove_gradient(b,0);
    //Initial condition for the modified problem
    if(input_obj.refine||input_obj.resume)
	gamma-=grad_art_press;
    else
	gamma=0;	
    //Solve the modified problem
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
    cout << "... done!" << endl;
    //Save Gamma
    cout << "Saving gamma(" << i << "," << j << ") ..." << endl;
    stringstream ss;
    ss << input_obj.runsname << "_gamma_" << i << "_" << j;
    rawsave_aux_field_hat(ss.str(),gamma);
    vtksave_aux_field_real(ss.str(),gamma);
    cout << "... done!" << endl;
    //Evaluate and save energy spectra
    cout << "Evaluating and saving energy spectrum for gamma(" << i << "," << j << ") ..." << endl;
    ss << "_spec.dat";
    ofstream ofs(ss.str().c_str());
    ofs << spectral_obj.eval_energ_spec(gamma,0) << endl;
    ofs.close();
    cout << "... done!" << endl;
    //Print energy spectra
    cout << "Energy spectrum (K Velocity Magnetic Temperature)" << endl;	
    cout << spectral_obj.eval_energ_spec(gamma,0) << endl;
    //   cout << "Printing non-vanishing harmonics in gamma" << endl;
    //   spectral_obj.pnvh(gamma);
}

//Evaluate ep (2x2 matrix to be diagonalised)
cat::Array<double,2> lss::eval_ep(const cat::Tvector<double,2> & q)
{
    cat::Array<double,2> out(2,2);
    cat::Array<double,2> e=eval_e(q);
    out(0,0)=e(0,0)*q[1]*q[1]-e(0,1)*q[0]*q[1]-e(1,0)*q[1]*q[0]+e(1,1)*q[0]*q[0];
    out(0,1)=e(0,2)*q[1]*q[1]-e(0,3)*q[0]*q[1]-e(1,2)*q[1]*q[0]+e(1,3)*q[0]*q[0];
    out(1,0)=e(2,0)*q[1]*q[1]-e(2,1)*q[0]*q[1]-e(3,0)*q[1]*q[0]+e(3,1)*q[0]*q[0];
    out(1,1)=e(2,2)*q[1]*q[1]-e(2,3)*q[0]*q[1]-e(3,2)*q[1]*q[0]+e(3,3)*q[0]*q[0];
    return out;
}

//Evaluate e (4x4 matrix to be reduced to 2x2)
cat::Array<double,2> lss::eval_e(const cat::Tvector<double,2> & q)
{
    cat::Array<double,2> out(4,4);
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

void lss::diag(std::complex<double> & lambda1,std::complex<double> & lambda2,const cat::Array<double,2> & matrix)
{
    double a=1.;
    double b=-matrix(0,0)-matrix(1,1);
    double c=matrix(0,0)*matrix(1,1)-matrix(0,1)*matrix(1,0);
    double delta=1.-4.*a*c/(b*b);
    lambda1=-b/(2.*a)*(1.+sqrt(delta));
    if(delta==0)
	lambda2=lambda1;
    else
	lambda2=c/(a*lambda1);
}


void lss::rawsave_aux_field_hat(const string & filename,const CBVF & field)
{
    rawFileSave(filename+"_vel_hat",field.vel());
    rawFileSave(filename+"_mag_hat",field.mag());
    rawFileSave(filename+"_temp_hat",field.temp());
}

void lss::vtksave_aux_field_real	(const string & filename,const CBVF & field)
{
    const int & n1 = input_obj.n1;
    const int & n2 = input_obj.n2;
    const int & n3 = input_obj.n3;
    cat::Tvector<int,3> dims(n1,n2,n3);
    cat::Tvector<double,3> ori(0,0,0);
    cat::Tvector<double,3> sp(input_obj.l1/input_obj.n1,input_obj.l2/input_obj.n2,input_obj.l3/input_obj.n3);
    RVF vf(n1,n2,n3);
    vf=0;
    spectral_obj.fft_ccs.inverse_transform(vf,field.vel());
    vtkFileSave(filename+"_vel",vf,ori,sp,"VelocityField");
    vf=0;
    spectral_obj.fft_ccs.inverse_transform(vf,field.mag());
    vtkFileSave(filename+"_mag",vf,ori,sp,"MagneticField");
    RSF sf(n1,n2,n3);
    sf=0;
    spectral_obj.sfft_s.inverse_transform(sf,field.temp());
    vtkFileSave(filename+"_temp",sf,ori,sp,"TemperatureField");
}

void lss::rawload_aux_field_hat(const string & filename,CBVF & field,const int & lr_n1,const int & lr_n2,const int & lr_n3)
{
    rawFileLoad(filename+"_vel_hat",field.vel(),lr_n1,lr_n2/2+1,lr_n3);
    rawFileLoad(filename+"_mag_hat",field.mag(),lr_n1,lr_n2/2+1,lr_n3);
    rawFileLoad(filename+"_temp_hat",field.temp(),lr_n1,lr_n2/2+1,lr_n3);
}

