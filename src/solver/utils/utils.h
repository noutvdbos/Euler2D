#pragma once

#include "Matrix.h"
#include "Vector.h"

//-----------------------------------------------------------------------
//   static class utils
//-----------------------------------------------------------------------


class utils
{   
    public:

    template <class T> static Vector<T> find(Matrix<T> mat, int idx, T val, const char* type);
    template <class T> static Vector<T> find(Matrix<T> mat, T val, const char* type);
    template <class T> static Vector<T> find(Matrix<T> mat, Vector<T> vals, const char* type );

    template <class T> static Vector<T> intersect(Matrix<T> mat1, Matrix<T> mat2 );
    template <class T> static Vector<T> intersect(Vector<T> vec1, Vector<T> vec2 );

    
    template <class T> static Vector<T> flip(Vector<T> vec );
    template <class T> static Vector<T> pow (Vector<T> vec, double pow );
    template <class T> static Vector<T> sqrt(Vector<T> vec);
    template <class T> static void      transpose (const Matrix<T>& mat, Matrix<T>& transMat);

    template <class T> static int       binarySearch(Vector<T> vec, int l, int r, T x);

    private:
    //disallow creating an instance of this class

    //Constructor
    utils();

    //Destructor
    ~utils();

    // private member functions

};

#include "utils.tpp"