
#ifndef MHDCL_SSS_H
#define MHDCL_SSS_H

#include "input.h"

#include <cat.h>

#include "spectral.h"

#include "linops.h"
#include "basic_fields.h"

#include <iostream>

using namespace std;

extern"C"
{ 

  void external_prodx_(int * nit,double * vi, double * vo, int * m);

  void vzdeigen_(double * v1,
		 double * xp,
		 double * eim,
		 double * ep,
		 double * thr,
		 int * M,
		 double * v2,
		 double * v3,
		 double * v4,
		 double * v5,
		 double * v6,
		 double * v7,
		 int * mp,
		 double * sc,
		 int * nseq,
		 const char * name);

  void vzdeigen_load_ffile_(double * v,const int * m,const char * name);

  void vzdeigen_save_ffile_(double * v,const int * m,const char * name);

} 

input * input_obj;
spectral * spectral_obj;
basic_fields * basic;
a_nought * a_nought_obj;

void gen_random_v(double * v,const int & sym,const int & seed);
void bv2v(double * v,const CBVF * bv);
void v2bv(CBVF * bv,const double * v);



#endif
