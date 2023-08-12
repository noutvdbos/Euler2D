#include <bits/stdc++.h>

#include <stdio.h>
#include "Vector.h"


//-----------------------------------------------------------------------
//   Constructors
//-----------------------------------------------------------------------

template <class T>
Vector<T>::Vector( )
{
    rows_ = 0;

    data_ = new T[0]; 
}

template <class T> 
Vector<T>::Vector( int rows )
{
    rows_ = rows;

    data_ = new T[rows]; 
}


template <class T> 
Vector<T>::Vector( int rows, T value )
{
    rows_ = rows;

    data_ = new T[rows]; 

    for(int i = 0; i < rows; i++) 
    {    
        data_[i] = value;   
    }
}


template <class T> 
Vector<T>::Vector( const Vector<T>& vec)
{
    rows_ = vec.rowSize();

    data_ = new T[rows_]; 

    for(int i = 0; i < rows_; i++) 
    {    
        data_[i] = vec[i];
    }
}


template <class T> 
Vector<T>::Vector( const std::list<T>& list)
{
    rows_ = list.size();

    data_ = new T[rows_]; 
    int idx = 0;

    for(auto num = list.begin(); num != list.end(); ++num) 
    {    
        data_[idx] = *num;
        idx++;
    }
}


template <class T> 
Vector<T>::Vector( const std::initializer_list<T>& list)
{
    rows_ = list.size();

    data_ = new T[rows_]; 
    int idx = 0;

    for(auto num = list.begin(); num != list.end(); ++num) 
    {    
        data_[idx] = *num;
        idx++;
    }
}


//-----------------------------------------------------------------------
//   Destructors
//-----------------------------------------------------------------------

template <class T> 
Vector<T>::~Vector( )
{
    delete[] data_;
}


//-----------------------------------------------------------------------
//  Member functions
//-----------------------------------------------------------------------


template <class T>
void Vector<T>::print () const
{

    for ( int i=0; i < rows_; i++ )  
    {
        std::cout << data_[i] << "\n";
    }
    std::cout<<std::endl;
}


template<class T> 
T& Vector<T>::operator[](int row) const
{
    return data_[row];
}


template<class T> 
Vector<T> Vector<T>::operator()(std::tuple<int,int> slice) const
{
    int start = std::get<0>(slice);
    int end   = std::get<1>(slice);

    assert(end>start);

    Vector<T> vec(end-start);

    int idx = 0;
    for ( int i = start; i< end; i++)
    {
        vec[idx] = data_[i];
        idx++;
    }
    
    return vec;
}


template<class T> 
Vector<T> Vector<T>::operator()(const Vector<int>& slice) const
{
    assert( slice.max() <= rows_ );

    Vector<T> vec(slice.rowSize());

    for ( int i = 0; i< slice.rowSize() ; i++)
    {
        vec[i] = data_[slice[i]];
    }
    
    return vec;
}


template<class T>
Vector<T>& Vector<T>::operator= ( const Vector<T>& vec)
{
    rows_ = vec.rowSize();
    
    delete[] data_;
    data_ = new T[rows_]; 

    for(int i = 0; i < rows_; i++) 
    {    
        data_[i] = vec[i];
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator= ( T val )
{
    for(int i = 0; i < rows_; i++) 
    {    
        data_[i] = val;
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator*= ( const Vector<T>& vec )
{
    assert ( rows_ == vec.rowSize() );

    for(int i = 0; i < rows_; i++) 
    {    
       this->data_[i] *=  vec[i];
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator*= ( T val )
{

    for(int i = 0; i < rows_; i++) 
    {    
        this->data_[i] *= val;
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator+= ( const Vector<T>& vec )
{
    for(int i = 0; i < vec.rowSize(); i++) 
    {    
        this->data_[i] += vec[i];
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator+= ( T val )
{
    for(int i = 0; i < rows_; i++) 
    {    
        this->data_[i] += val;
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator-= ( const Vector<T>& vec )
{
    
    for(int i = 0; i < vec.rowSize(); i++) 
    {    
        this->data_[i] -= vec[i];
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator-= ( T val )
{
    for(int i = 0; i < rows_; i++) 
    {    
        this->data_[i] -= val;
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator/= ( const Vector<T>& vec )
{
    
    for(int i = 0; i < vec.rowSize(); i++) 
    {    
        this->data_[i] /= vec[i];
    }

    return *this;
}


template<class T>
Vector<T>& Vector<T>::operator/= ( T val )
{

    for(int i = 0; i < rows_; i++) 
    {    
       this->data_[i] /= val;
    }

    return *this;
}


template <class T>
void Vector<T>::reinit( int rows )
{
    reinit( rows, (T) 0 );
}


template <class T>
void Vector<T>::reinit (int rows, T value)
{
    rows_ = rows;

    delete[] data_;
    data_ = new T[rows]; 

    for (int i = 0; i < rows; i++) 
    {
        data_[i] = value; 
    }

}


template <class T>
void Vector<T>::resize (int rows)
{
    if (rows < rows_)
    {
        
        rows_ = rows;
        return;
    }

    T* newData = new T[rows];

    for (int i = 0; i < rows_; i++) 
    {
        newData[i] = data_[i];
    }
    for (int i = rows_; i < rows; i++)
    {
        newData[i] = (T) 0;
    }
    delete[] data_;
    data_ = newData;
    rows_ = rows;
}


template <class T>
void Vector<T>::sort ()
{
    std::sort(data_, data_+rows_);
}


template <class T>
double Vector<T>::norm () const
{
    double norm = 0;

    for (int i = 0; i<rows_; i++)
    {
        norm += std::pow(data_[i],2);
    }
    norm = std::sqrt(norm);

    return norm;
}


template <class T>
T Vector<T>::sum () const
{
    T sum = (T) 0;

    for (int i = 0; i<rows_; i++)
    {
        sum += data_[i];
    }

    return sum;
}


template <class T>
T Vector<T>::min () const
{
    return *std::min_element(data_,data_+rows_);
}


template <class T>
T Vector<T>::max () const
{
    return *std::max_element(data_,data_+rows_);
}


template <class T>
T Vector<T>::mean () const
{
    return sum()/rows_;
}



// Operator definitions

template <class T>
inline Vector<T> operator*( Vector<T> left, Vector<T> const& right ) 
{
  return left*=right;
}
template <class T>
inline Vector<T> operator*(  Vector<T> left, T right ) 
{
  return left*=right;
}
template <class T>
inline Vector<T> operator*(  T left, Vector<T> right ) 
{
  return right*=left;
}


template <class T>
inline Vector<T> operator+( Vector<T> left, Vector<T> const& right ) 
{
  return left+=right;
}
template <class T>
inline Vector<T> operator+(  Vector<T> left, T right ) 
{
  return left+=right;
}
template <class T>
inline Vector<T> operator+(  T left, Vector<T> right ) 
{
  return right+=left;
}


template <class T>
inline Vector<T> operator-( Vector<T> left, Vector<T> const& right ) 
{
  return left-=right;
}
template <class T>
inline Vector<T> operator-(  Vector<T> left, T right ) 
{
  return left-=right;
}
template <class T>
inline Vector<T> operator-(  T left, Vector<T> right ) 
{
  return -right+=left;
}


template <class T>
inline Vector<T> operator/( Vector<T> left, Vector<T> const& right ) 
{
  return left/=right;
}
template <class T>
inline Vector<T> operator/(  Vector<T> left, T right ) 
{
  return left/=right;
}
template <class T>
Vector<T> operator/(  T left, Vector<T> right ) 
{
  Vector<T> vec(right.rowSize());
  for (int i = 0; i < right.rowSize(); i++)
  {
      vec[i] = left/right[i];
  }
  return vec;
}