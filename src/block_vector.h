// -*- C++ -*-
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



//block_vector.h

#ifndef BLOCK_VECTOR_H
					#define BLOCK_VECTOR_H

#include <cat.h>
using namespace cat;


template <class T>
class BlockVector
  {
    //Members
  private:
    cat::Array<cat::Tvector<T,3>,3> velocity_;
    cat::Array<cat::Tvector<T,3>,3> magnetic_;
    cat::Array<T,3> temperature_;
    int sym_;

    //Accessors
  public:

    cat::Array<cat::Tvector<T,3>,3> & vel()
    {
      return velocity_;
    }
    cat::Array<cat::Tvector<T,3>,3> & mag()
    {
      return magnetic_;
    }
    cat::Array<T,3> & temp()
    {
      return temperature_;
    }

    const cat::Array<cat::Tvector<T,3>,3> & vel() const
      {
        return velocity_;
      }
    const cat::Array<cat::Tvector<T,3>,3> & mag() const
      {
        return magnetic_;
      }
    const cat::Array<T,3> & temp() const
      {
        return temperature_;
      }

    int & sym()
    {
      return sym_;
    }
    const int & sym() const
      {
        return sym_;
      }

    Tvector<int,3> & shape()
    {
      return this->temperature_.shape();
    }
    const Tvector<int,3> & shape() const
      {
        return this->temperature_.shape();
      }

    int & size()
    {
      return this->temperature_.size();
    }
    const int & size() const
      {
        return this->temperature_.size();
      }

    //Public Methods

    // 	void save(const string & filename)
    // 	{
    // 		ofstream ofs(filename);
    // 		ofs << "# vtk DataFile Version 2.0\n"
    // 			<< ".\n"
    // 			<< "ASCII"	<< endl;
    // 			ofs << "DATASET STRUCTURED_POINTS\n"
    // 			<< "DIMENSIONS " << n1 << " " << n2 << " " << n3 << "\n"
    // 			<< "ORIGIN 0 0 0"
    // 			<< "SPACING " << l1/n1 << " " << l2/n2 << " " << l3/n3 << "\n"
    // 			<< "POINT_DATA " << this->velocity_.size() << "\n"
    // 			<< "VECTORS velocity float" << endl;
    // 			for(int k=0;k<n3;++k)
    // 				for(int j=0;j<n2;++j)
    // 					for(int i=0;i<n1;++i)
    // 						ofs << this->velocity_(i,j,k) << endl;
    // 		ofs << "DATASET STRUCTURED_POINTS\n"
    // 			<< "DIMENSIONS " << n1 << " " << n2 << " " << n3 << "\n"
    // 			<< "ORIGIN 0 0 0"
    // 			<< "SPACING " << l1/n1 << " " << l2/n2 << " " << l3/n3 << "\n"
    // 			<< "POINT_DATA " << this->magnetic_.size() << "\n"
    // 			<< "VECTORS magnetic float" << endl;
    // 		for(int k=0;k<n3;++k)
    // 			for(int j=0;j<n2;++j)
    // 				for(int i=0;i<n1;++i)
    // 					ofs << this->magnetic_(i,j,k) << endl;
    // 		ofs << "DATASET STRUCTURED_POINTS\n"
    // 			<< "DIMENSIONS " << n1 << " " << n2 << " " << n3 << "\n"
    // 			<< "ORIGIN 0 0 0"
    // 			<< "SPACING " << l1/n1 << " " << l2/n2 << " " << l3/n3 << "\n"
    // 			<< "POINT_DATA " << this->temperature_.size() << "\n"
    // 			<< "SCALARS temperature float\n"
    // 			<< "LOOKUP_TABLE default" << endl;
    // 		for(int k=0;k<n3;++k)
    // 			for(int j=0;j<n2;++j)
    // 				for(int i=0;i<n1;++i)
    // 					ofs << this->temperature_(i,j,k) << endl;
    // 	};



    void enforce_sym()
    {
      if (sym>1)
        return;
      enforce_sym_();
    }

    void enforce_sym(const int & sym__)
    {
      sym_=sym__;
      if (sym_>1)
        return;
      enforce_sym_();
    }

  private:
    void enforce_sym_()
    {
      typename cat::Array<cat::Tvector<T,3>,3>::iterator vel_it(this->vel());
      typename cat::Array<cat::Tvector<T,3>,3>::iterator mag_it(this->mag());
      typename cat::Array<T,3>::_iterator temp_it(this->temp());
      for(vel_it=(this->vel()).begin(),
          mag_it=(this->mag()).begin(),
          temp_it=(this->temp()).begin();
          vel_it!=(this->vel()).end(),
          mag_it!=(this->mag()).end(),
          temp_it!=(this->temp()).end();
          ++vel_it,
          ++mag_it,
          ++temp_it)
        {
          if(this->sym_==0)
            {
              for(int m=0;m<2;++m)
                {
                  ((*vel_it)[m]).imag()=0;
                  ((*mag_it)[m]).imag()=0;
                }
              ((*vel_it)[2]).real()=0;
              ((*mag_it)[2]).real()=0;
              (*temp_it).real()=0;
            }
          else if (this->sym_==1)
            {
              for(int m=0;m<2;++m)
                {
                  ((*vel_it)[m]).real()=0;
                  ((*mag_it)[m]).real()=0;
                }
              ((*vel_it)[2]).imag()=0;
              ((*mag_it)[2]).imag()=0;
              (*temp_it).imag()=0;

            }
        }

      //     //symmetry about z axis
      //     int s1=this->vel().shape()[0];
      //     int s3=this->vel().shape()[2];
      //     for(int i=1;i<s1/2+1;++i)
      //       for(int k=0;k<s3;++k)
      // 	{
      // 	  for(int m=0;m<2;++m)
      // 	      {
      // 		(this->vel())(s1-i,0,k)[m]=-(sym_?1:-1)*(this->vel())(i,0,k)[m];
      // 		this->mag()(s1-i,0,k)[m]=-(sym_?1:-1)*this->mag()(i,0,k)[m];
      // 	      }
      // 	  this->vel()(s1-i,0,k)[2]=(sym_?1:-1)*this->vel()(i,0,k)[2];
      // 	  this->mag()(s1-i,0,k)[2]=(sym_?1:-1)*this->mag()(i,0,k)[2];
      // 	  this->temp()(s1-i,0,k)=(sym_?1:-1)*this->temp()(i,0,k);
      // 	}


    }

    //Forbidden constructors
  private:
    //default constructor
    BlockVector()
    {}
    ;

  public:
    //constructor from size
    BlockVector(int s1,int s2,int s3):
        velocity_(s1,s2,s3),
        magnetic_(s1,s2,s3),
        temperature_(s1,s2,s3),
        sym_(2)
    {
      //velocity_=0;
      //magnetic_=0;
      //temperature_=0;
    };

    //constructor from shape
    BlockVector(const Tvector<int,3> & shape__):
        velocity_(shape__),
        magnetic_(shape__),
        temperature_(shape__),
        sym_(2)
    {
      //velocity_=0;
      //magnetic_=0;
      //temperature_=0;
    };


    //constructor from shape and pointer
    BlockVector(const Tvector<int,3> & shape__, double * data__):
				velocity_(shape__,reinterpret_cast<cat::Tvector<T,3> *>(data__)),
        magnetic_(shape__,reinterpret_cast<cat::Tvector<T,3> *>(data__+3*2*velocity_.size())),
        temperature_(shape__,reinterpret_cast<T *>(data__+2*3*2*velocity_.size())),
        sym_(2)

    {
      //cout << "BV" << endl;
      //cout << shape__ << endl;
			//cout << "VS " << velocity_.size() << endl;
			//exit(0);

    }
    ;




    //contuctor from existing fields and symmetry (duplicates data)
    BlockVector(const cat::Array<cat::Tvector<T,3>,3> & velocity__,const cat::Array<cat::Tvector<T,3>,3> & magnetic__,const cat::Array<T,3> & temperature__,const int & sym__):
        velocity_(velocity__),
        magnetic_(magnetic__),
        temperature_(temperature__),
        sym_(sym__)
    {}
    ;


    //copy constructor
    BlockVector(const BlockVector & rhs):
        velocity_(rhs.vel()),
        magnetic_(rhs.mag()),
        temperature_(rhs.temp()),
        sym_(rhs.sym())
    {}
    ;

    //destructor
    ~BlockVector()
    {}
    ;


    //Public methods
  public:

    //IO
    //redefiniton of cout
    friend std::ostream& operator<<(std::ostream& output,const BlockVector& ovector)
    {
      output << ovector.vel() << ovector.mag() << ovector.temp();
      return output;
    };
    friend std::ostream& operator<<(std::ostream& output,BlockVector& ovector)
    {
      output << ovector.vel() << ovector.mag() << ovector.temp();
      return output;
    };
    //redefiniton of cin
    friend std::istream& operator>>(std::istream& input,BlockVector& ivector)
    {
      input >> ivector.vel() >> ivector.mag() >> ivector.temp();
      return input;
    };


    //operators

    //assignment
    BlockVector & operator=(const BlockVector & rhs)
    {
      velocity_=rhs.vel();
      magnetic_=rhs.mag();
      temperature_=rhs.temp();
      sym_=rhs.sym();
      return *this;
    }

    BlockVector & operator=(const T & rhs)
    {
      velocity_=rhs;
      magnetic_=rhs;
      temperature_=rhs;
      sym_=2;
      return *this;
    }


#define SL_BLOCK_VECTOR_UPDATE(op)		\
							inline BlockVector &				\
						operator op(const BlockVector & rhs)		\
						{						\
												this->velocity_ op rhs.vel();		\
												this->magnetic_ op rhs.mag();		\
												this->temperature_ op rhs.temp();		\
												return *this;				\
													}						\
					template <class AT>				\
					inline BlockVector &			\
					operator op(const AT & rhs)			\
					{						\
											this->velocity_ op rhs;			\
											this->magnetic_ op rhs;			\
											this->temperature_ op rhs;			\
											return *this;				\
												}

    SL_BLOCK_VECTOR_UPDATE(+=);
    SL_BLOCK_VECTOR_UPDATE(-=);
    SL_BLOCK_VECTOR_UPDATE(*=);
    SL_BLOCK_VECTOR_UPDATE(/=);

  };


//Declaration & definition of unary -
template <class T>
inline BlockVector<T>
operator -(const BlockVector<T>& rhs)
{
  BlockVector<T> aux(rhs);
  aux.vel()=-rhs.vel();
  aux.mag()=-rhs.mag();
  aux.temp()=-rhs.temp();
  return aux;
};


// Binary operators for BlockVectors
//Needs working


#define BLOCK_VECTOR_BINARY(op)						\
					template <class T1>							\
							inline BlockVector<T1>						\
							operator op(const BlockVector<T1> & lhs,const BlockVector<T1> & rhs) \
							{									\
														BlockVector<T1> out(lhs.shape());					\
														out.vel() = lhs.vel() op rhs.vel();					\
																	out.mag() = lhs.mag() op rhs.mag();					\
																	out.temp() = lhs.temp() op rhs.temp();				\
																				 return out;								\
																		 }									\
					template <class T1,class T2>						\
					inline BlockVector<T1>						\
					operator op(const BlockVector<T1> & lhs,const T2 & rhs)		\
					{									\
														BlockVector<T1> out(lhs.shape());					\
														out.vel() = lhs.vel() op rhs;					\
																	out.mag() = lhs.mag() op rhs;					\
																	out.temp() = lhs.temp() op rhs;					\
																				 return out;								\
																		 }									\
					template <class T1,class T2>						\
					inline BlockVector<T2>						\
					operator op(const T1 & lhs,const BlockVector<T2> & rhs)		\
					{									\
														BlockVector<T2> out(rhs.shape());					\
														out.vel() = lhs op rhs.vel();					\
																	out.mag() = lhs op rhs.mag();					\
																	out.temp() = lhs op rhs.temp();					\
																				 return out;								\
																		 }

BLOCK_VECTOR_BINARY(+);
BLOCK_VECTOR_BINARY(-);
BLOCK_VECTOR_BINARY(*);
BLOCK_VECTOR_BINARY(/);




#endif
