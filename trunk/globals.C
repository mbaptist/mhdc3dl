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
