// -*- C++ -*-

#ifndef GLOBALS_H
#define GLOBALS_H
 
#include <complex>

#include <cat.h>

#include "block_vector.h"


//TYPE DEFINITIONS

//typedef float Real;//Real scalars
typedef double Real;//Real scalars
//typedef long double Real;//Real scalars
typedef std::complex<Real> Complex;//Complex scalars
typedef cat::tvector<Real,3> RV;//Real vectors
typedef cat::tvector<Complex,3> CV;//Complex Vectors
typedef cat::array<Real,3> RSF;//Real Scalar Fields
typedef cat::array<Complex,3> CSF;//Complex Scalar Fields
typedef cat::array<RV,3> RVF;//Real Vector Fields
typedef cat::array<CV,3> CVF;//Complex Vector Fields
typedef block_vector<Complex > CBVF;//Complex Block Vector

//GLOBAL CONSTANTS


//Symbolic Constants for Layer Problems
//representation
//scalars
const int S=0;
const int C=1;
//vectors
const int CCS=1;
const int SSC=0;
//symmetries
//scalars
const int SYM=1;
const int ASYM=0;
const int SSYM=SYM;
const int SASYM=ASYM;
//vectors
const int AAS=1;
const int SSA=0;
const int VSYM=AAS;
const int VASYM=SSA;

//GLOBAL FUNCTIONS

//Unit vector along direction i (i=0,1,2)
cat::tvector<int,3> unit_vector(const int & i);

//Element i,j of Levi-Civita tensor (i,j=0,1,2)
int levi_civita(const int & i,const int & j,const int & k);


#endif
