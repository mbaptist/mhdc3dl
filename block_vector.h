// -*- C++ -*-

//block_vector.h

#ifndef BLOCK_VECTOR_H
#define BLOCK_VECTOR_H

#include <cat.h>
using namespace cat;


template <class T>
class block_vector
{
  //Members
private:
  cat::array<cat::tvector<T,3>,3> velocity_;
  cat::array<cat::tvector<T,3>,3> magnetic_;
  cat::array<T,3> temperature_;
  int sym_;

  //Accessors
public:

  cat::array<cat::tvector<T,3>,3> & vel(){return velocity_;}
  cat::array<cat::tvector<T,3>,3> & mag(){return magnetic_;}
  cat::array<T,3> & temp(){return temperature_;}

  const cat::array<cat::tvector<T,3>,3> & vel() const {return velocity_;}
  const cat::array<cat::tvector<T,3>,3> & mag() const {return magnetic_;}
  const cat::array<T,3> & temp() const {return temperature_;}

  int & sym(){return sym_;}
  const int & sym() const {return sym_;}

  tvector<int,3> & shape(){return this->temperature_.shape();}
  const tvector<int,3> & shape() const {return this->temperature_.shape();}

  //Public Methods

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
    cat::array_iterator<cat::tvector<T,3>,3> vel_it(this->vel());
    cat::array_iterator<cat::tvector<T,3>,3> mag_it(this->mag());
    cat::array_iterator<T,3> temp_it(this->temp());
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
  block_vector(){};

public:
  //constructor from size
  block_vector(int s1,int s2,int s3):
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
  block_vector(const tvector<int,3> & shape__):
    velocity_(shape__),
    magnetic_(shape__),
    temperature_(shape__),
    sym_(2)
  {
    //velocity_=0;
    //magnetic_=0;
    //temperature_=0;
  };

  //contuctor from existing fields and symmetry (duplicates data)
  block_vector(const cat::array<cat::tvector<T,3>,3> & velocity__,const cat::array<cat::tvector<T,3>,3> & magnetic__,const cat::array<T,3> & temperature__,const int & sym__):
    velocity_(velocity__),
    magnetic_(magnetic__),
    temperature_(temperature__),
    sym_(sym__)
  {};

  //copy constructor
  block_vector(const block_vector & rhs):
    velocity_(rhs.vel()),
    magnetic_(rhs.mag()),
    temperature_(rhs.temp()),
    sym_(rhs.sym())
  {};

  //destructor
  ~block_vector(){};
  

  //Public methods
public:

  //IO
  //redefiniton of cout  
  friend std::ostream& operator<<(std::ostream& output,const block_vector& ovector)
  {
    output << ovector.vel() << ovector.mag() << ovector.temp();
    return output;
  };
  friend std::ostream& operator<<(std::ostream& output,block_vector& ovector)
  {
    output << ovector.vel() << ovector.mag() << ovector.temp();
    return output;
  };
  //redefiniton of cin
  friend std::istream& operator>>(std::istream& input,block_vector& ivector)
  {
    input >> ivector.vel() >> ivector.mag() >> ivector.temp();
    return input;
  }; 


  //operators

  //assignment
  block_vector & operator=(const block_vector & rhs)
  {
    velocity_=rhs.vel();
    magnetic_=rhs.mag();
    temperature_=rhs.temp();
    sym_=rhs.sym();
    return *this;
  }

  block_vector & operator=(const T & rhs)
  {
    velocity_=rhs;
    magnetic_=rhs;
    temperature_=rhs;
    sym_=2;
    return *this;
  }


#define SL_BLOCK_VECTOR_UPDATE(op)		\
  inline block_vector &				\
  operator op(const block_vector & rhs)		\
  {						\
    this->velocity_ op rhs.vel();		\
    this->magnetic_ op rhs.mag();		\
    this->temperature_ op rhs.temp();		\
    return *this;				\
  }						\
    template <class AT>				\
    inline block_vector &			\
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
inline block_vector<T>
operator -(const block_vector<T>& rhs)
{   
  block_vector<T> aux(rhs);
  aux.vel()=-rhs.vel();
  aux.mag()=-rhs.mag();
  aux.temp()=-rhs.temp();
  return aux;
};


// Binary operators for block_vectors
//Needs working


#define BLOCK_VECTOR_BINARY(op)						\
  template <class T1>							\
  inline block_vector<T1>						\
  operator op(const block_vector<T1> & lhs,const block_vector<T1> & rhs) \
  {									\
    block_vector<T1> out(lhs.shape());					\
    out.vel() = lhs.vel() op rhs.vel();					\
    out.mag() = lhs.mag() op rhs.mag();					\
    out.temp() = lhs.temp() op rhs.temp();				\
    return out;								\
  }									\
  template <class T1,class T2>						\
  inline block_vector<T1>						\
  operator op(const block_vector<T1> & lhs,const T2 & rhs)		\
  {									\
    block_vector<T1> out(lhs.shape());					\
    out.vel() = lhs.vel() op rhs;					\
    out.mag() = lhs.mag() op rhs;					\
    out.temp() = lhs.temp() op rhs;					\
    return out;								\
  }									\
  template <class T1,class T2>						\
  inline block_vector<T2>						\
  operator op(const T1 & lhs,const block_vector<T2> & rhs)		\
  {									\
    block_vector<T2> out(rhs.shape());					\
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
