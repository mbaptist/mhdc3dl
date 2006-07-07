
//C/C++ declaration of the interface of vzdeigen (implemented in vzdeigen.f)

#ifndef VZDEIGEN_H
#define VZDEIGEN_H

extern "C"
{

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

#endif

