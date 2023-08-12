#include <iostream>
#include <string>
#include <fstream>
#include "utils.h"
#include "Matrix.h"
#include "Vector.h"
#include "utils/Types.h"


//=======================================================================
//   static class utils
//=======================================================================

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//   find()
//-----------------------------------------------------------------------


// Find value in a certain row and return the column indices

template <class T>
Vector<T> utils::find(Matrix<T> mat, int idx, T val, const char* type )
{
    std::list<int> list; 
    int i, end;

    if ( type == Types::CELL )
    {
        i = 1;
        end = mat(idx,0) + 1;
    }
    else
    {
        i = 0;
        end = mat.colSize();
    }

    for (; i < end; i++)
    {
        if( mat(i,idx) == val)
        {
            list.push_back(i);
        }
    }
    Vector<T> results(list);

    return results;
}


//Find a value and return the row indices

template <class T>
Vector<T> utils::find(Matrix<T> mat, T val, const char* type )
{
    std::list<int> list; 
    for (int i = 0; i < mat.rowSize(); i++)
    {
        int j, end;

        if ( type == Types::CELL )
        {
            j = 1;
            end = mat(i,0) + 1;
        }
        else
        {
            j = 0;
            end = mat.colSize();
        }

        for (; j < end; j++)
        {
            if( mat(i,j) == val)
            {
                list.push_back(i);
            }
        }
    }
    Vector<T> results(list);

    return results;
}


//Find the rows that has all the input values, return the indices

template <class T>
Vector<T> utils::find(Matrix<T> mat, Vector<T> vals, const char* type )
{
    std::list<int> list; 
    for (int i = 0; i < mat.rowSize(); i++)
    {
        int hasVals = 0;
        
        for (int ivals = 0; ivals < vals.rowSize(); ivals++)
        {
            int j, end;

            if ( type == Types::CELL )
            {
                j = 1;
                end = mat(i,0) + 1;
            }
            else
            {
                j = 0;
                end = mat.colSize();
            }

            for (; j < end; j++)
            {
                if( mat(i,j) == vals(ivals))
                {
                    hasVals++;
                    break;
                }
            }
        }
        if (hasVals == vals.rowSize())
        {
            list.push_back(i);
        }
    }
    Vector<T> results(list);

    return results;
}


//-----------------------------------------------------------------------
//   intersect()
//-----------------------------------------------------------------------


template <class T>
Vector<T> utils::intersect (Vector<T> vec1, Vector<T> vec2 )
{
    std::list<T> list; 

    // Before finding intersection, make sure vec1
    // is smaller
    if ( vec1.rowSize() > vec2.rowSize() ) {
        

        Vector<T> temp = vec1;
        vec2 = vec1;
        vec1 = temp;
    }

    // Sort smaller array vec1
    vec1.sort();
 
    // Search every element of bigger array in smaller
    // array and print the element if found
    for (int i = 0; i < vec2.rowSize(); i++)
    {
        if ( utils::binarySearch(vec1, 0, vec1.rowSize() - 1, vec2(i)) != -1 )
        {
            list.push_back(vec2(i));
        }
    }

    Vector<T> results(list);

    return results;    
}


//-----------------------------------------------------------------------
//   flip()
//-----------------------------------------------------------------------


template <class T>
Vector<T> utils::flip (Vector<T> vec )
{
    Vector<T> results(vec.rowSize());

    for (int i = 0; i < vec.rowSize(); i++)
    {
        results[i] = vec[vec.rowSize()-i-1];
    }
    return results;
}


//-----------------------------------------------------------------------
//   pow()
//-----------------------------------------------------------------------


template <class T>
Vector<T> utils::pow (Vector<T> vec, double pow )
{
    Vector<T> results(vec.rowSize());

    for (int i = 0; i < vec.rowSize(); i++)
    {
        results[i] = std::pow(vec[i],pow) ;
    }
    return results;
}


//-----------------------------------------------------------------------
//   sqrt()
//-----------------------------------------------------------------------


template <class T>
Vector<T> utils::sqrt (Vector<T> vec )
{
    Vector<T> results(vec.rowSize());

    for (int i = 0; i < vec.rowSize(); i++)
    {
        results[i] = std::sqrt(vec[i]);
    }
    return results;
}


//-----------------------------------------------------------------------
//   sqrt()
//-----------------------------------------------------------------------


template <class T>
void utils::transpose (const Matrix<T>& mat, Matrix<T>& transMat)
{
    assert(mat.rowSize() == transMat.colSize());
    assert(mat.colSize() == transMat.rowSize());

    for (int i = 0; i < mat.rowSize(); i++)
    {
        for ( int j = 0; j < mat.colSize(); j++)
        {
            transMat[j][i] = mat[i][j];
        }
    }
    
    return;
}


//-----------------------------------------------------------------------
//   binarySearch()
//-----------------------------------------------------------------------


template <class T>
int utils:: binarySearch(Vector<T> vec, int l, int r, T x)
{
    if (r >= l) {
        int mid = l + (r - l) / 2;
 
        // If the element is present at the middle itself
        if (vec(mid) == x)
            return mid;
 
        // If element is smaller than mid, then it can only
        // be present in left subarray
        if (vec(mid) > x)
            return binarySearch(vec, l, mid - 1, x);
 
        // Else the element can only be present in right
        // subarray
        return binarySearch(vec, mid + 1, r, x);
    }
 
    // We reach here when element is not present in array
    return -1;
}