
// VTK File Format Header
// # vtk DataFile Version 2.0
// .
// ASCII

// VTK Scalar DATASET
// DATASET STRUCTURED_POINTS
// DIMENSIONS 32 32 16
// ORIGIN 0 0 0
// SPACING 0.19635 0.19635 0.19635
// POINT_DATA 16384
// SCALARS temperature float
// LOOKUP_TABLE default

// VTK Vector DATASET
// DATASET STRUCTURED_POINTS
// DIMENSIONS 32 32 16
// ORIGIN 0 0 0
// SPACING 0.19635 0.19635 0.19635
// POINT_DATA 16384
// VECTORS FieldName float

#ifndef VTKIO_H
#define VTKIO_H

#include <string>
#include <cat.h>
#include "block_vector.h"

template <class T>
struct vtkFileTraits
{
	enum {fieldtype};
};

template <>
struct vtkFileTraits<cat::Array<double,3> >
{
	enum {fieldtype=0};
};

template <>
struct vtkFileTraits<cat::Array<cat::Tvector<double,3>,3> >
{
	enum {fieldtype=1};
};
	                                     

template <class T>
void vtkFileLoad(const std::string & vtkfilename,T & data);

template <class T>
void vtkFileSave(const std::string & vtkfilename,const T & data,cat::Tvector<double,3> _vtkfile_origin_,cat::Tvector<double,3> _vtkfile_spacing_,std::string _vtkfile_fieldname_);

template <class T>
void rawFileLoad(const std::string & filename,T & data);

template <class T>
void rawFileLoad(const std::string & filename,T & data,const int & lr_n1,const int & lr_n2,const int & lr_n3);

template <class T>
void rawFileSave(const std::string & filename,const T & data);

#include "io.C"

#endif


