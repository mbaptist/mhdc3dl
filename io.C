//#include "io.h"
#include <string>
#include <iostream>
#include <fstream>

#include <cat.h>
using namespace std;
using namespace cat;

template <class T>
void vtkFileLoad(const std::string & vtkfilename,T & data)
{
	//temporary strings
	string stmp;
	//input stream
	ifstream ifs((vtkfilename+".vtk").c_str());
  //read file header
	ifs >> stmp;ifs >> stmp;ifs >> stmp;ifs >> stmp;ifs >> stmp;
	ifs >> stmp;
	ifs >> stmp;
	ifs >> stmp;ifs >> stmp;
	cat::tvector<int,3> dimensions;
	ifs >> stmp;
	ifs >> dimensions;
	ifs >> stmp;ifs >> stmp;ifs >> stmp;ifs >> stmp;
	cout << stmp << endl;
	ifs >> stmp;ifs >> stmp;ifs >> stmp;ifs >> stmp;
	cout << stmp << endl;
	ifs >> stmp;ifs >> stmp;
	cout << stmp << endl;
	string fieldtype;
	ifs >> fieldtype;
	cout << fieldtype << endl;
	ifs >> stmp;ifs >> stmp;
	cout << stmp << endl;
	if (fieldtype=="SCALARS")
	{
		ifs >> stmp;
		cout << stmp << endl;
		ifs >> stmp;
		cout << stmp << endl;
	}
	//read data
	int n1=dimensions[0];
	int n2=dimensions[1];
	int n3=dimensions[2];
	T data_read(n1,n2,n3);
	//cout << dimensions << endl;
	for(int k=0;k<n3;++k)
		for(int j=0;j<n2;++j)
			for(int i=0;i<n1;++i)
			{
				ifs >> data_read(i,j,k);
				//cout << data_read(i,j,k) << endl;
			}
	data=data_read;
}


template <class T>
void vtkFileSave(const std::string & vtkfilename,const T & data,cat::tvector<double,3> _vtkfile_origin_,cat::tvector<double,3> _vtkfile_spacing_,std::string _vtkfile_fieldname_)
{
   //output stream
	ofstream ofs((vtkfilename+".vtk").c_str());
	//Header data
	std::string vtkfile_version="2.0";
	std::string vtkfile_creator=".";
	std::string vtkfile_encoding="ASCII";
	std::string	vtkfile_datasettype="STRUCTURED_POINTS";
	cat::tvector<int,3>	vtkfile_dimensions=data.shape();
	int	vtkfile_size=data.size();
	cat::tvector<double,3> vtkfile_origin=_vtkfile_origin_;
	cat::tvector<double,3> vtkfile_spacing=_vtkfile_spacing_;
	std::string vtkfile_fieldtype;
	if (vtkFileTraits<T>::fieldtype==0)
			vtkfile_fieldtype="SCALARS";
		else
			vtkfile_fieldtype="VECTORS";
	std::string vtkfile_fieldname=_vtkfile_fieldname_;
	std::string vtkfile_numerictype="float";
   //write file header
	ofs << "# vtk DataFile Version " << vtkfile_version << "\n"
		<< vtkfile_creator << "\n"
		<< vtkfile_encoding	<< endl;
	//write dataset header
	ofs << "DATASET " << vtkfile_datasettype << "\n"
		<< "DIMENSIONS " << vtkfile_dimensions << "\n"
		<< "ORIGIN " << vtkfile_origin << "\n"
		<< "SPACING " << vtkfile_spacing << "\n"
		<< "POINT_DATA " << vtkfile_size << "\n"
		<< vtkfile_fieldtype << " " << vtkfile_fieldname << " " << vtkfile_numerictype << endl;
	if (vtkfile_fieldtype=="SCALARS")
		ofs << "LOOKUP_TABLE default" << endl;
	//write data
	int n1=vtkfile_dimensions[0];
	int n2=vtkfile_dimensions[1];
	int n3=vtkfile_dimensions[2];
	for(int k=0;k<n3;++k)
		for(int j=0;j<n2;++j)
			for(int i=0;i<n1;++i)
				ofs << data(i,j,k) << endl;
}

template <class T>
void rawFileLoad(const std::string & filename,T & data)
{
	cout << filename+".raw" << endl;
	  //open input stream
	ifstream ifs((filename+".raw").c_str());
	//output the field
	ifs >> data;
	//close output stream
	ifs.close();
}

template <class T>
void rawFileLoad(const std::string & filename,T & data,const int & lr_s1,const int & lr_s2,const int & lr_s3)
{
	cat::tvector<int,3> lr_shape=cat::tvector<int,3>(lr_s1,lr_s2,lr_s3);
	cout << lr_shape << endl;
	cout << data.shape() << endl;
	if(lr_shape==data.shape())
		rawFileLoad(filename,data);
	else
	{
		T aux(lr_shape);
		rawFileLoad(filename,aux);
		for(int i=0;i<lr_s1/2+1;++i)
			for(int j=0;j<lr_s2;++j)
				for(int k=0;k<lr_s3;++k)
					data(i,j,k)=aux(i,j,k);
		for(int i=lr_s1/2+1;i<lr_s1;++i)
			for(int j=0;j<lr_s2;++j)
				for(int k=0;k<lr_s3;++k)
					data(i+data.shape()[0]-lr_s1,j,k)=aux(i,j,k);
	}
}

template <class T>
void rawFileSave(const std::string & filename,const T & data)
{
  //open output stream
	ofstream ofs((filename+".raw").c_str());
	//output the field
	ofs << data << endl;
	//close output stream
	ofs.close();
}

