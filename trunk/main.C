#include "mhdc3dl.h"

int main(int argc,char * argv[])
{
 
  //get input filename
  string inputfname("default.cfg");
  if (argc>1)
    if (argv[1]!="")
      inputfname=argv[1];
  cout << inputfname << endl;

  double lambda_minimal,lambda_maximal;

  input input_obj(inputfname);

  lss lss_obj(input_obj);

  lss_obj.run(lambda_minimal,lambda_maximal);

  //Output maximum growth rate
  //To screen
  cout << "  Maximum Lambda = " << lambda_maximal
       << "  Minimum Lambda = " << lambda_minimal
       << endl;
  //Append to file
//   if (input_obj->grfname!="")
//     {
//       ofstream ofs;
//       ofs.open((input_obj->grfname).c_str(),ofstream::app);
//       ofs << lambda_maximal 
// 	  << " " << lambda_minimal
// 	  << endl;
//       ofs.close();
//     }  

  return 0;
}
