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

#include "mhdc3dl.h"
#include <iostream>
#include <string>
#include <cstring>
using namespace std;

int main(int argc,char * argv[])
{
  //get Input filename
  string inputfname("mhdc3dl_run_default.py");
  if (argc>1)
    if (strcmp(argv[1],""))
      inputfname=argv[1];
  cout << inputfname << endl;
  //load Input paarmeters
  Input input_obj(inputfname);
  //create lss object
  lss lss_obj(input_obj);
  //evaluate large-scale maximum and minimum eigenvalues
  double theta_min,theta_max;
  std::complex<double> lambda_min,lambda_max;
  lss_obj.run(theta_min,lambda_min,theta_max,lambda_max);
  //print maximum and minimum growth rates
  cout << "  Minimum Lambda = " << lambda_min << " for theta=" << theta_min << endl;
  cout << "  Maximum Lambda = " << lambda_max << " for theta=" << theta_max << endl;

  // 	sss sss_obj(input_obj);
  // 	double xp,eim;
  // 	sss_obj.run(xp,eim);

  return 0;
}
