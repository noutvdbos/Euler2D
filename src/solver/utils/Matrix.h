#pragma once

#include "Vector.h"

//-----------------------------------------------------------------------
//   class Matrix
//-----------------------------------------------------------------------

template <class T>
class Matrix
{

public:
    //constructor
    Matrix              ();
    Matrix              ( int rows, int cols );
    Matrix              ( Vector<int> matSize );
    Matrix              ( const Matrix<T>& mat );
    Matrix              ( const std::list<Vector<T>>& list);
    Matrix              ( const std::initializer_list<Vector<T>>& list); 

    //Destructor
    ~Matrix();

    //member functions
    void            print          ();
    T&              operator()     ( int row, int col );
    Vector<T>&      operator[]     ( int row ) const;
    Vector<T>       operator()     ( int row, std::tuple<int,int> slice ) const;
    Matrix<T>&      operator=      ( const Matrix<T>& mat );
    Matrix<T>&      operator=      ( std::list<Vector<T>>& list );
    Matrix<T>&      operator=      ( std::initializer_list<Vector<T>> list); 
    Vector<T>       operator*      ( const Vector<T>& vec );
    Matrix<T>&      operator+=     ( const Matrix<T>& mat );

    Vector<T>       col            ( int col);

    void            reinit         ( int rows, int cols );
    void            reinit         ( int rows, int cols, T value );
    void            reinit         ( int rows, const std::list<T>& list );
    void            reinit         ( int rows, const Vector<T>& values );
    void            resize         ( int rows );
    Vector<int>     size           ();

    Matrix<T>       transpose      ();
    Matrix<T>&      elemwiseMultiply ( const Vector<T>& vec );

    // copy functions
    void            copy           (int idx, const Vector<T>& vec);
    void            copy           (int idxRow, int idxStart, const Vector<T>& vec);

    //get functions (1 liner so it is defined here)
    int rowSize() const {return rows_;}
    int colSize() const {return cols_;}

private:

    // variables
    Vector< Vector<T> > data_;
    int rows_,cols_;

};

#include "Matrix.tpp"