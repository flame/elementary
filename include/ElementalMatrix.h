/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef ELEMENTAL_MATRIX_H
#define ELEMENTAL_MATRIX_H 1

#include "ElementalEnvironment.h"

namespace Elemental
{
    template<typename T>
    class Matrix
    {
        bool      _viewing;
        bool      _lockedView;
        int       _height;
        int       _width;
        int       _ldim;
        T*        _data;
        const T*  _lockedData;
        Memory<T> _memory;

    public:    

        Matrix();                

        Matrix( const int height, const int width );

        Matrix( const int height, const int width, const int ldim );

        Matrix( const Matrix<T>& A );

        ~Matrix();
        
        void Print() const;
        void Print(std::string msg) const;

        void SetToRandom();

        T& operator() ( const int i, const int j );
        T  operator() ( const int i, const int j ) const;
        
        int Height() const;
        int Width() const;
        int LDim() const;
        int MemorySize() const;

        T* Buffer();
        T* Buffer
        ( const int i, const int j );
        T* Buffer
        ( const int i, const int j, const int height, const int width );

        T* Pointer();
        T* Pointer( const int i, const int j );

        const T* LockedBuffer() const;
        const T* LockedBuffer
        ( const int i, const int j ) const;
        const T* LockedBuffer
        ( const int i, const int j, const int height, const int width ) const;

        const T* LockedPointer() const;
        const T* LockedPointer( const int i, const int j ) const;

        // Resize the matrix
        void ResizeTo( const int height, const int width );
        void ResizeTo( const int height, const int width, const int ldim );

        void View
        ( Matrix<T>& A);

        void View
        ( Matrix<T>& A, 
          const int i, const int j, const int height, const int width );

        void View1x2( Matrix<T>& AL, Matrix<T>& AR );

        void View2x1( Matrix<T>& AT, 
                      Matrix<T>& AB );

        void View2x2( Matrix<T>& ATL, Matrix<T>& ATR,
                      Matrix<T>& ABL, Matrix<T>& ABR );
        
        void LockedView( const Matrix<T>& A );

        void LockedView
        ( const Matrix<T>& A, 
          const int i, const int j, const int height, const int width );

        void LockedView1x2( const Matrix<T>& AL, const Matrix<T>& AR );

        void LockedView2x1( const Matrix<T>& AT, 
                            const Matrix<T>& AB );

        void LockedView2x2( const Matrix<T>& ATL, const Matrix<T>& ATR,
                            const Matrix<T>& ABL, const Matrix<T>& ABR );

        void SetToIdentity();

        void SetToZero();

        const Matrix<T>& operator=( const Matrix<T>& A );
    };
}

/*----------------------------------------------------------------------------*/

template<typename T>
inline
Elemental::Matrix<T>::Matrix()
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _ldim(0), _data(NULL), _lockedData(NULL),
  _memory()
{ }

template<typename T>
inline
Elemental::Matrix<T>::Matrix( const int height, const int width )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _ldim(std::max(height,1)), 
  _lockedData(NULL)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix(height,width)");
    if( height < 0 || width < 0 )
    {
        std::cerr << "Height and width must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    _memory.Require( _ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline 
Elemental::Matrix<T>::Matrix
( const int height, const int width, const int ldim )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _ldim(ldim), _lockedData(NULL)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix(height,width,ldim)");
    if( height < 0 || width < 0 )
    {
        std::cerr << "Height and width must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( ldim < height )
    {
        std::cerr << "Initialized with ldim(" << ldim << ") < "
                  << "height(" << height << ")." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( ldim == 0 )
    {
        std::cerr << "Leading dimensions cannot be zero (for BLAS compat.)."
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    _memory.Require( ldim*width );
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::Matrix<T>::Matrix
( const Matrix<T>& A )
: _viewing(false), _lockedView(false), _lockedData(false)
{
#ifndef RELEASE
    PushCallStack("Matrix::Matrix( const Matrix& )");
#endif
    if( &A != this )
    {
        *this = A;
    }
    else
    {
        std::cout << "You just tried to construct a Matrix with itself..." 
                  << std::endl;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::Matrix<T>::~Matrix()
{ }

template<typename T>
inline T& 
Elemental::Matrix<T>::operator()
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( i > Height() || j > Width() )
    {
        std::cerr << "Out of bounds: "
                  << "(" << i << "," << j << ") of " << Height() 
                  << " x " << Width() << " Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return _data[i+j*_ldim];
}

template<typename T>
inline T
Elemental::Matrix<T>::operator()
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::operator");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( i > Height() || j > Width() )
    {
        std::cerr << "Out of bounds: "
                  << "(" << i << "," << j << ") of " << Height()
                  << " x " << Width() << " Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    if( _lockedData )
        return _lockedData[i+j*_ldim];
    else
        return _data[i+j*_ldim];
}

template<typename T>
inline int
Elemental::Matrix<T>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::Matrix<T>::Width() const
{ return _width; }

template<typename T>
inline int
Elemental::Matrix<T>::LDim() const
{ return _ldim; }

template<typename T>
inline int
Elemental::Matrix<T>::MemorySize() const
{ return _memory.Size(); }

template<typename T>
inline T*
Elemental::Matrix<T>::Buffer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( _lockedView )
    {
        std::cerr << "Cannot return non-const buffer of locked Matrix." 
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedBuffer() const
{
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
Elemental::Matrix<T>::Buffer
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( _lockedView )
    {
        std::cerr << "Cannot return non-const buffer of locked Matrix."
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    // The height or width of the buffer could be zero, so we 
    // use strict inequalities for flexibility. Pointer() does not.
    if( i>_height || j>_width )
    {
        std::cerr << "Requested out-of-bounds buffer of Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedBuffer
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    // The height or width of the buffer could be zero, so we 
    // use strict inequalities for flexibility. LockedPointer() does not.
    if( i>_height || j>_width )
    {
        std::cerr << "Requested out-of-bounds buffer of Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
Elemental::Matrix<T>::Buffer
( const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::Buffer");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( _lockedView )
    {
        std::cerr << "Cannot return non-const buffer of locked Matrix."
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( (i+height)>_height || (j+width)>_width )
    {
        std::cerr << "Requested out-of-bounds buffer of Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedBuffer
( const int i, const int j, const int height, const int width ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedBuffer");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( height < 0 || width < 0 )
    {
        std::cerr << "Height and width must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( (i+height)>_height || (j+width)>_width )
    {
        std::cerr << "Requested out-of-bounds buffer of Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

template<typename T>
inline T*
Elemental::Matrix<T>::Pointer()
{
#ifndef RELEASE
    PushCallStack("Matrix::Pointer");
    if( _lockedView )
    {
        std::cerr << "Cannot return non-const pointer to locked Matrix."
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( _height == 0 || _width == 0 )
    {
        std::cerr << "Requested pointer to empty Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return _data;
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedPointer() const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedPointer");
    if( _height == 0 || _width == 0 )
    {
        std::cerr << "Requested locked pointer to empty Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    if( _lockedView )
        return _lockedData;
    else
        return _data;
}

template<typename T>
inline T*
Elemental::Matrix<T>::Pointer
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("Matrix::Pointer");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( _lockedView )
    {
        std::cerr << "Cannot return non-const pointer to locked Matrix." 
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( i >= _height || j >= _width )
    {
        std::cerr << "Requested out-of-bounds pointer to Matrix." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return &_data[i+j*_ldim];
}

template<typename T>
inline const T*
Elemental::Matrix<T>::LockedPointer
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedPointer");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( i >= _height || j >= _width )
    {
        std::cerr << "Requested out-of-bounds locked pointer to Matrix." 
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    if( _lockedView )
        return &_lockedData[i+j*_ldim];
    else
        return &_data[i+j*_ldim];
}

#endif /* ELEMENTAL_MATRIX_H */

