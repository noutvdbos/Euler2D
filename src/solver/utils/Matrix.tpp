#include <stdio.h>
#include "Matrix.h"
#include "Vector.h"


//-----------------------------------------------------------------------
//   Constructors
//-----------------------------------------------------------------------

template <class T>
Matrix<T>::Matrix( )
{
    reinit(0, 0);
}


template <class T> 
Matrix<T>::Matrix( int rows, int cols )
{
    reinit(rows, cols);
}


template <class T> 
Matrix<T>::Matrix( Vector<int> matSize )
{
    reinit(matSize[0],matSize[1]);
}


template <class T> 
Matrix<T>::Matrix( const Matrix<T>& mat )
{
    rows_ = mat.rowSize();
    cols_ = mat.colSize();


    for (int i = 0; i < rows_; i++) 
    {
        data_[i] = Vector<T>(cols_); 
    }

    for(int i = 0; i < rows_; i++) 
    {
        for(int j = 0; j < cols_; j++)
        {
            data_[i][j] = mat[i][j];
        }
    }
}


template<class T>
Matrix<T>::Matrix ( const std::list<Vector<T>>& list)
{
    rows_ = list.size();
    cols_ = list.front()->rowSize();

    int idx = 0;

    data_.reinit(rows_);

    for(auto num = list.begin(); num != list.end(); ++num) 
    {
        data_[idx] = *num;
        idx++;
    }
}


template<class T>
Matrix<T>::Matrix ( const std::initializer_list<Vector<T>>& list)
{
    rows_ = list.size();
    cols_ = (*list.begin()).rowSize();
    
    int idx = 0;

    data_.reinit(rows_);

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
Matrix<T>::~Matrix( )
{
    
}


//-----------------------------------------------------------------------
//  Member functions
//-----------------------------------------------------------------------


template <class T>
void Matrix<T>::print ()
{
    int i,j;

    for (i=0;i < rows_;i++) 
    {
        for(j=0;j < cols_;j++) 
        {
            printf("%.4f    ",(float) data_[i][j]);
        }
        printf("\n");
    }
    std::cout<<std::endl;
}


template<class T> 
T& Matrix<T>::operator()(int row, int col)
{
    return data_[row][col];
}


template<class T> 
Vector<T>& Matrix<T>::operator[] (int row) const 
{
    return data_[row];
}


template<class T> 
Vector<T> Matrix<T>::operator()(int row, std::tuple<int,int> slice) const
{
    int start = std::get<0>(slice);
    int end   = std::get<1>(slice);

    assert(end>start);

    Vector<T> vec(end-start);

    int idx = 0;
    for ( int i = start; i< end; i++)
    {
        vec[idx] = data_[row][i];
        idx++;
    }
    
    return vec;
}


template<class T>
Matrix<T>& Matrix<T>::operator= ( const Matrix<T>& mat)
{
    rows_ = mat.rowSize();
    cols_ = mat.colSize();
    

    data_.reinit(rows_);

    for(int i = 0; i < rows_; i++) 
    {
        data_[i] = Vector<T>(cols_);
        
        for(int j = 0; j < cols_; j++)
        {
            data_[i][j] = mat[i][j];
        }
    }

    return *this;
}


template<class T>
Matrix<T>& Matrix<T>::operator= ( std::list<Vector<T>>& list)
{
    rows_ = list.size();
    cols_ = list.front()->rowSize();

    int idx = 0;

    data_.reinit(rows_);

    for(auto num = list.begin(); num != list.end(); ++num) 
    {
        data_[idx] = *num;
        idx++;
    }

    return *this;
}


template<class T>
Matrix<T>& Matrix<T>::operator= ( std::initializer_list<Vector<T>> list)
{
    rows_ = list.size();
    cols_ = (*list.begin()).rowSize();
    
    int idx = 0;

    data_.reinit(rows_);

    for(auto num = list.begin(); num != list.end(); ++num) 
    {
        data_[idx] = *num;
        idx++;
    }

    return *this;
}


template<class T>
Vector<T> Matrix<T>::operator* ( const Vector<T>& vec )
{
    assert ( vec.rowSize() == cols_ );

    Vector<T> result(vec.rowSize(), 0.0);
    for (int i = 0; i < rows_; i++ )
    {
        for (int j = 0; j < cols_; j++ )
        {
            result[i] += data_[i][j]*vec[j];
        }
    }

    return result;
}


template<class T>
Matrix<T>& Matrix<T>::operator+= (const Matrix<T>& mat)
{
    assert (mat.colSize() == this->cols_);
    assert (mat.rowSize() == this->rows_);

    for (int i = 0; i < this->rows_; i++)
    {
        for (int j = 0; j < this->cols_; j++)
        {
           this->data_[i][j] += mat[i][j];
        }
    }

    return *this;
}


template <class T>
Vector<T> Matrix<T>::col( int col )
{
    Vector<T> vec (cols_);

    for (int i = 0; i < rows_; i++ )
    {
        vec[i] = data_[i][col];
    }

    return vec;
}


template <class T>
void Matrix<T>::reinit( int rows, int cols )
{
    reinit( rows, cols, (T) 0 );
}


template <class T>
void Matrix<T>::reinit( int rows, int cols, T value )
{
    rows_ = rows;
    cols_ = cols;

    data_.reinit(rows);

    for (int i = 0; i < rows; i++) 
    {
        data_[i] = Vector<T>(cols); 
    }

    for(int i = 0; i < rows; i++) 
    {
        for(int j = 0; j < cols; j++) 
        {
            data_[i][j] = value;
        }
    }
}


template <class T>
void Matrix<T>::reinit( int rows, const Vector<T>& values )
{

    rows_ = rows;
    cols_ = values.rowSize();

    data_.reinit(rows);

    for (int i = 0; i < rows; i++) 
    {
       data_[i] = Vector<T>(cols_); 
    }

    for(int i = 0; i < rows; i++) 
    {
        for(int j = 0; j < values.rowSize(); j++) 
        {
            data_[i][j] = values[j];
        }
    }
}


template<class T>
void Matrix<T>::reinit( int rows, const std::list<T>& list )
{
    rows_ = rows;
    cols_ = list.size();

    data_.reinit(rows);

    for (int i = 0; i < rows; i++) 
    {
        data_[i] = Vector<T>(cols_);  
    }

    for(int i = 0; i < rows; i++) 
    {
        int j = 0;
        for(auto num = list.begin(); num != list.end(); ++num) 
        {
            data_[i][j] = *num;
            j++;
        }
    } 
}


template <class T>
void Matrix<T>::resize (int rows)
{
    if (rows < rows_)
    {
        rows_ = rows;
        return;
    }

    Vector< Vector<T> > temp = data_;

    data_.reinit(rows);

    for (int i = 0; i < rows; i++) 
    {
       data_[i] = Vector<T>(cols_); 
    }

    for(int i = 0; i < rows; i++) 
    {
        for(int j = 0; j < cols_; j++) 
        {
            data_[i][j] = temp[i][j];
        }
    }
    rows_ = rows;
}


template<class T>
Vector<int> Matrix<T>::size ()
{
    Vector<int> matSize( 2 );

    matSize[0] = rows_;
    matSize[1] = cols_;

    return matSize;
}


// Note, this returns a full matrix, thus it should only be used for small 
// matrices. For large matrices, use the transpose function of the utils class.

template<class T>
Matrix<T> Matrix<T>::transpose()
{
    int rows = cols_;
    int cols = rows_;
   
    Matrix<T> mat(rows,cols);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            mat(i,j) = data_[j][i];
        }
    }

    return mat;
}


template<class T>
Matrix<T>& Matrix<T>::elemwiseMultiply (const Vector<T>& vec)
{
    assert ( rows_ = vec.rowSize() );

    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
        {
            data_[i][j] *= vec[i];
        }
    }

    return *this;
}



//-----------------------------------------------------------------------
//   copy()
//-----------------------------------------------------------------------

// This function copies a vector to a row of the matrix. If the vector
// is larger than the row size, it is not copied. If the vector is 
// smaller or equal to the row size, it is copied. If the vector only has 
// one entry, it is copied to the whole row.

template<class T>
void Matrix<T>::copy
(
    int         idx,
    const Vector<T>&   vec
)
{
    if (vec.rowSize() > cols_)
    {
        std::cout<< "can not copy vector into matrix: Vector is too large\n";
        return;
    }

    if ( vec.rowSize() == 1 )
    {
        for (int i = 0; i < cols_; i++)
        {
            data_[idx][i] = vec[0];
        }
    }

    else
    {
        for (int i = 0; i < vec.rowSize(); i++)
        {
            data_[idx][i] = vec[i];
        }
    }
}


// Since a start idx is given, it is copied from here untill the vector length,
// unless the vector is too large for the matrix row.

template<class T>
void Matrix<T>::copy
(
    int         idxRow,
    int         idxStart,
    const Vector<T>&  vec
)
{
    if (vec.rowSize() + idxStart > cols_)
    {
        std::cout<< "can not copy vector into matrix: Vector is too large\n";
        return;
    }

    for (int i = idxStart; i < vec.rowSize() + idxStart; i++)
    {
        data_[idxRow][i] = vec[i-idxStart];
    }

}