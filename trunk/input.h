// -*- C++ -*-


//input.h

//Class input declares the input parameters and is inherited
//by the problem class ns3d
//It is responsible for reading the configuration file and 
//seting the input parameters 

#ifndef INPUT_H
#define INPUT_H

#include <Python.h>

#include <iostream>
#include <string>

using namespace std;

//Forward declaration of class py_input_parser 
class py_input_parser;


///////////////////////////////////////////
// Declaration of class input
///////////////////////////////////////////
class input
{
public: 
  //constructor from filename
  input(std::string cfg_fname);
  //constructor from python object
  input(PyObject * module);
  //destructor
  ~input();
 public:
  string runsname;
  //input parameters
  int n1,n2,n3;//Spatial grid
  double l1,l2,l3,ell1,ell2;//Spatial extent
  //physical Parameters
  double visc,omegaz,compress,g,diff,deltat,econd,tcond;
  
  //basic fields
  std::string basic_mode;
  std::string basic_vel_fname;
  std::string basic_mag_fname;
  std::string basic_temp_fname;
	//Random Mode
	int br_seed,br_ki,br_kf;
	double 	br_alpha,br_rms_norm;
	bool	br_kind,br_sym;


  //large scale stability
  
  double ls_eps;
  double qq;
  double qq_adj;
  int kk;
  double small;
  double small_adj;

  std::string apbfname;

  std::string grfname;
  
  //short scale stability
  int sym_sub,mp,nseq,sss_seed,sss_basic_restore_seed;
  double ep,thr,sc;
  std::string sss_ifname,sss_int_ofbname;

 private:
  //load from config file
  void py_load_input(py_input_parser & parse);
};    


///////////////////////////////////////////
// Declaration of class py_input_parser
///////////////////////////////////////////
class py_input_parser
{
private:
  //Filename to parse
  std::string filename;
  //Python objects
  PyObject * module;
  PyObject * dict;
  PyObject * pval(const string & item);
  PyThreadState * inter;
  bool python_initialised;
public:
  //Constructor from filename
  py_input_parser(string filename_);
  //Constructor from dict
  py_input_parser(PyObject * module_);
  //Destructor
  ~py_input_parser();
  //Python to C++ conversion function
  void operator()(int & val,const std::string & item);
	void operator()(bool & val,const std::string & item);
  void operator()(double & val,const std::string & item);
  void operator()(std::string & val,const std::string & item);
private:
  //Forbidden constructors
  py_input_parser();//default
  py_input_parser(const py_input_parser &);//copy
};




#endif
