 
#include "globals.h"

#include <cat.h>

using namespace cat;

//GLOBAL FUNCTIONS (IMPLEMENTATION)

//Unit vector along direction i (i=0,1,2)
cat::Tvector<int,3> unit_vector(const int & i)
{
  //assert(i>=0&&i<3)
  cat::Tvector<int,3> aux;
  aux=0;
  aux[i]=1;
  return aux;
}

//Element i,j of Levi-Civita tensor (i,j=0,1,2)
int levi_civita(const int & i,const int & j,const int & k)
{
  if (i==j || j==k || k==i)
    {
      return 0;
    }
  else
    {
      if ( (i==1 && j==2 && k==3) ||
	   (i==2 && j==3 && k==1) ||
	   (i==3 && j==1 && k==2) )
	{
	  return 1;
	}
      else
	{
	  return -1;
	}
    }
}
