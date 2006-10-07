#include "mhdc3dl.h"
#include <iostream>
#include <string>
#include <cstring>
using namespace std;

int main(int argc,char * argv[])
{ 
  //get input filename
	string inputfname("mhdc3dl_run_default.py");
	if (argc>1)
		if (strcmp(argv[1],""))
			inputfname=argv[1];
	cout << inputfname << endl;
	//load input paarmeters	
	input input_obj(inputfname);
	//create lss object
	lss lss_obj(input_obj);
	//evaluate large-scale maximum and minimum eigenvalues
	double theta_min,lambda_min,theta_max,lambda_max;
	lss_obj.run(theta_min,lambda_min,theta_max,lambda_max);
  //print maximum and minimum growth rates
	cout << "  Minimum Lambda = " << lambda_min << " for theta=" << theta_min << endl;
	cout << "  Maximum Lambda = " << lambda_max << " for theta=" << theta_max << endl;
		
// 	sss sss_obj(input_obj);
// 	double xp,eim;
// 	sss_obj.run(xp,eim);

	return 0;
}
