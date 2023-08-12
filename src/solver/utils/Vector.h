#pragma once

#include <list>

//-----------------------------------------------------------------------
//   class Config
//-----------------------------------------------------------------------

template <class T>
class Vector
{

public:

    //constructor
    Vector              ();
    Vector              ( int rows );
    Vector              ( int rows, T value );
    Vector              ( const Vector<T>& vec );
    Vector              ( const std::list<T>& list); 
    Vector              ( const std::initializer_list<T>& list);

    //Destructor
    ~Vector();

    //member functions
    void            print          () const;
    T&              operator[]     (int row) const;
    Vector<T>       operator()     ( std::tuple<int,int> slice ) const;
    Vector<T>       operator()     ( const Vector<int>&  slice ) const;
    Vector<T>&      operator=      ( const Vector<T>& vec );
    Vector<T>&      operator=      ( T val );
    Vector<T>&      operator*=     ( const Vector<T>& vec );
    Vector<T>&      operator*=     ( T val );
    Vector<T>&      operator+=     ( const Vector<T>& vec );
    Vector<T>&      operator+=     ( T val );
    Vector<T>&      operator-=     ( const Vector<T>& vec );
    Vector<T>&      operator-=     ( T val );
    Vector<T>&      operator/=     ( const Vector<T>& vec );
    Vector<T>&      operator/=     ( T val );

    void            reinit         ( int rows );
    void            reinit         ( int rows, T value );
    void            resize         ( int rows );

    void            sort           ();
    int             adress         (int row);
    double          norm           () const;
    T               sum            () const;
    T               min            () const;
    T               max            () const;
    T               mean           () const;
    
    //get functions (1 liner so it is defined here)
    int rowSize() const {return rows_;}
    


private:


    // variables
    T  *data_;
    int rows_;

};

#include "Vector.tpp"