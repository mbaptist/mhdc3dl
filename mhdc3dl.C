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
	
	input input_obj(inputfname);


	lss lss_obj(input_obj);
	
	double lambda_min,lambda_max;
	
	lss_obj.run(lambda_min,lambda_max);
	
  //print maximum and minimum growth rates
	cout << "  Maximum Lambda = " << lambda_max
		<< "  Minimum Lambda = " << lambda_min
		<< endl;


// 	sss sss_obj(input_obj);
// 	double xp,eim;
// 	sss_obj.run(xp,eim);

	
	return 0;
}
