/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_DISTMATRIX_STAR_VR_H
#define ELEMENTAL_DISTMATRIX_STAR_VR_H 1

#include "ElementalDistMatrix.h"

namespace Elemental
{
    // Partial specialization to A[* ,VR]
    //
    // The rows of these distributed matrices are spread throughout the 
    // process grid in a row-major fashion, while the columns are not 
    // distributed.
    template<typename T>
    class DistMatrix<T,Star,VR> 
    {
        bool      _viewing;
        bool      _lockedView;
        int       _height;
        int       _width;
        Memory<T> _auxMemory;
        Matrix<T> _localMatrix;

        bool _constrainedRowDist;
        int  _rowAlignment;
        int  _rowShift;
        const Grid* _grid;

    public:

        DistMatrix
        ( const Grid& grid );

        DistMatrix
        ( const int height, const int width, const Grid& grid );

        DistMatrix
        ( const bool constrainedRowDist, const int rowAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const int height, const int width,
          const bool constrainedRowDist, const int rowAlignment,
          const Grid& grid                                      );

        DistMatrix
        ( const DistMatrix<T,Star,VR>& A );

        ~DistMatrix();

        //--------------------------------------------------------------------//
        // Operations that can be performed by individual processes           //
        //--------------------------------------------------------------------//

        const Grid& GetGrid() const;

        bool Viewing() const;

        // Matrix dimensions
        int Height() const;
        int Width() const;
        int LocalHeight() const;
        int LocalWidth() const;
        int LocalLDim() const;

        // Retrieve (a reference to) an entry from the local matrix
        T& LocalEntry( const int i, const int j );
        T  LocalEntry( const int i, const int j ) const;

        // Return an (immutable) reference to the local matrix
              Matrix<T>& LocalMatrix();
        const Matrix<T>& LockedLocalMatrix() const;

        // Generic distribution parameters
        bool ConstrainedColDist() const;
        bool ConstrainedRowDist() const;
        int  ColAlignment() const;
        int  RowAlignment() const;
        int  ColShift() const;
        int  RowShift() const;

        //--------------------------------------------------------------------//
        // Operations that must be collectively performed                     //
        //--------------------------------------------------------------------//

        // Get/Set an entry from the distributed matrix
        T    Get( const int i, const int j );
        void Set( const int i, const int j, const T u );
 
        // Zero out necessary entries to make distributed matrix trapezoidal:
        //
        //   If side equals 'Left', then the diagonal is chosen to pass through 
        //   the upper-left corner of the matrix.
        //
        //   If side equals 'Right', then the diagonal is chosen to pass through
        //   the lower-right corner of the matrix.
        //
        // Upper trapezoidal with offset = 0:
        //
        //    |x x x x x x x| <-- side = Left      |0 0 x x x x x|
        //    |0 x x x x x x|                      |0 0 0 x x x x|
        //    |0 0 x x x x x|     side = Right --> |0 0 0 0 x x x|
        //    |0 0 0 x x x x|                      |0 0 0 0 0 x x|
        //    |0 0 0 0 x x x|                      |0 0 0 0 0 0 x|
        //
        // Upper trapezoidal with offset = 1:
        //   
        //    |0 x x x x x x| <-- side = Left      |0 0 0 x x x x|
        //    |0 0 x x x x x|                      |0 0 0 0 x x x|
        //    |0 0 0 x x x x|     side = Right --> |0 0 0 0 0 x x|
        //    |0 0 0 0 x x x|                      |0 0 0 0 0 0 x|
        //    |0 0 0 0 0 x x|                      |0 0 0 0 0 0 0|
        //
        // Lower trapezoidal with offset = 1:
        //    
        //    |x x 0 0 0 0 0| <-- side = Left      |x x x x 0 0 0|
        //    |x x x 0 0 0 0|                      |x x x x x 0 0|
        //    |x x x x 0 0 0|     side = Right --> |x x x x x x 0|
        //    |x x x x x 0 0|                      |x x x x x x x|
        //    |x x x x x x 0|                      |x x x x x x x|
        //
        void MakeTrapezoidal
        ( const Side side, const Shape shape, const int offset = 0 );

        void Print( const std::string msg ) const;
        void ResizeTo( const int height, const int width );
        void SetToIdentity();
        void SetToRandom();
        void SetToRandomDiagDominant();
        void SetToZero();
       
        // For aligning the row and/or column distributions with another matrix.
        // Often useful when two distributed matrices are added together.
        //
        // The top part of this list contains the (valid) distributions that
        // contain 'VectorRow'.
        void AlignWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignWith( const DistMatrix<T,MR,  Star>& A );
        void AlignWith( const DistMatrix<T,Star,MR  >& A );
        void AlignWith( const DistMatrix<T,Star,VR  >& A );
        void AlignWith( const DistMatrix<T,VR,  Star>& A );
        void AlignRowsWith( const DistMatrix<T,MC,  MR  >& A );
        void AlignRowsWith( const DistMatrix<T,MR,  MC  >& A );
        void AlignRowsWith( const DistMatrix<T,MR,  Star>& A );
        void AlignRowsWith( const DistMatrix<T,Star,MR  >& A );
        void AlignRowsWith( const DistMatrix<T,Star,VR  >& A );
        void AlignRowsWith( const DistMatrix<T,VR,  Star>& A );
        // These are no-ops, but they exist for template flexibility
        void AlignWith( const DistMatrix<T,Star,MC  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,MD  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,VC  >& A ) {}
        void AlignWith( const DistMatrix<T,Star,Star>& A ) {}
        void AlignWith( const DistMatrix<T,MC,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,MD,  Star>& A ) {}
        void AlignWith( const DistMatrix<T,VC,  Star>& A ) {}
        void AlignColsWith( const DistMatrix<T,Star,MC  >& A ) {}
        void AlignColsWith( const DistMatrix<T,Star,MD  >& A ) {}
        void AlignColsWith( const DistMatrix<T,Star,MR  >& A ) {}
        void AlignColsWith( const DistMatrix<T,Star,VC  >& A ) {}
        void AlignColsWith( const DistMatrix<T,Star,VR  >& A ) {}
        void AlignColsWith( const DistMatrix<T,Star,Star>& A ) {}
        void AlignColsWith( const DistMatrix<T,MC,  Star>& A ) {}
        void AlignColsWith( const DistMatrix<T,MD,  Star>& A ) {}
        void AlignColsWith( const DistMatrix<T,MR,  Star>& A ) {}
        void AlignColsWith( const DistMatrix<T,VC,  Star>& A ) {}
        void AlignColsWith( const DistMatrix<T,VR,  Star>& A ) {}

        // So that matrix-multiplication will make sense, we force alignment
        // with a single distribution type that can be inferred.
        void ConformWith( const DistMatrix<T,Star,VR>& A );
        void ConformWith( const DistMatrix<T,VR,Star>& A );
        // This is a no-op, but it exists for template flexibility
        void ConformWith( const DistMatrix<T,Star,Star>& A ) {}

        // Clear the alignment constraints
        void FreeConstraints();

        // (Immutable) view of a distributed matrix
        void View( DistMatrix<T,Star,VR>& A );
        void LockedView( const DistMatrix<T,Star,VR>& A );

        // (Immutable) view of a portion of a distributed matrix
        void View
        ( DistMatrix<T,Star,VR>& A,
          const int i, const int j, const int height, const int width );

        void LockedView
        ( const DistMatrix<T,Star,VR>& A,
          const int i, const int j, const int height, const int width );

        // (Immutable) view of two horizontally contiguous partitions of a
        // distributed matrix
        void View1x2
        ( DistMatrix<T,Star,VR>& AL, DistMatrix<T,Star,VR>& AR );

        void LockedView1x2
        ( const DistMatrix<T,Star,VR>& AL, const DistMatrix<T,Star,VR>& AR );

        // (Immutable) view of two vertically contiguous partitions of a
        // distributed matrix
        void View2x1
        ( DistMatrix<T,Star,VR>& AT,
          DistMatrix<T,Star,VR>& AB );

        void LockedView2x1
        ( const DistMatrix<T,Star,VR>& AT,
          const DistMatrix<T,Star,VR>& AB );

        // (Immutable) view of a contiguous 2x2 set of partitions of a 
        // distributed matrix
        void View2x2
        ( DistMatrix<T,Star,VR>& ATL, DistMatrix<T,Star,VR>& ATR,
          DistMatrix<T,Star,VR>& ABL, DistMatrix<T,Star,VR>& ABR );

        void LockedView2x2
        ( const DistMatrix<T,Star,VR>& ATL, const DistMatrix<T,Star,VR>& ATR,
          const DistMatrix<T,Star,VR>& ABL, const DistMatrix<T,Star,VR>& ABR );
 
        // Bury communication behind '=' operator
        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,MC,MR>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,MC,Star>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,Star,MR>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,MD,Star>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,Star,MD>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,MR,MC>& A );
        
        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,MR,Star>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,Star,MC>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,VC,Star>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,Star,VC>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,VR,Star>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,Star,VR>& A );

        const DistMatrix<T,Star,VR>&
        operator=( const DistMatrix<T,Star,Star>& A );
    };
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::DistMatrix
( const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedRowDist(false), _rowAlignment(0), _rowShift(grid.VRRank()),
  _grid(&grid)
{ }

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::DistMatrix
( const int height, const int width, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedRowDist(true), _rowAlignment(0), _rowShift(grid.VRRank()),
  _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::DistMatrix(height,width)");
    if( height < 0 || width < 0 )
    {
        if( grid.VCRank() == 0 )
            std::cerr << "Height and width must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    _localMatrix.ResizeTo
    ( height, utilities::LocalLength( width, grid.VRRank(), grid.Size() ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::DistMatrix
( const bool constrainedRowDist, const int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(0), _width(0), _auxMemory(), _localMatrix(),
  _constrainedRowDist(constrainedRowDist),
  _rowAlignment(rowAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::DistMatrix(rowAlign)");
    if( rowAlignment < 0 || rowAlignment >= grid.Size() )
    {
        if( grid.VCRank() == 0 )
        {
            std::cerr << "rowAlignment for [* ,VR] must be in [0,p-1], p=rxc."
                      << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
#endif
    _rowShift = utilities::Shift( grid.VRRank(), rowAlignment, grid.Size() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::DistMatrix
( const int height, const int width,
  const bool constrainedRowDist, const int rowAlignment, const Grid& grid )
: _viewing(false), _lockedView(false),
  _height(height), _width(width), _auxMemory(),
  _constrainedRowDist(constrainedRowDist),
  _rowAlignment(rowAlignment), _grid(&grid)
{ 
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::DistMatrix(m,n,rowAlign)");
    if( height < 0 || width < 0 )
    {
        if( grid.VCRank() == 0 )
            std::cerr << "Height and width must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( rowAlignment < 0 || rowAlignment >= grid.Size() )
    {
        if( grid.VCRank() == 0 )
        {
            std::cerr << "rowAlignment for [* ,VR] must be in [0,p-1], p=rxc."
                      << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
#endif
    _rowShift = utilities::Shift( grid.VRRank(), _rowAlignment, grid.Size() );
    _localMatrix.ResizeTo
    ( height, utilities::LocalLength(width,_rowShift,grid.Size()) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::DistMatrix
( const DistMatrix<T,Elemental::Star,Elemental::VR>& A )
: _viewing(false), _lockedView(false),
  _constrainedRowDist(A.ConstrainedRowDist()),
  _rowAlignment(A.RowAlignment()), _rowShift(A.RowShift()),
  _grid( &( A.GetGrid() ) )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::DistMatrix( const DistMatrix[* ,VR]& )");
#endif
    if( &A != this )
    {
        *this = A;
    }
    else
    {
        std::cerr << "You just tried to construct a DistMatrix[* ,VR] with"
                  << " itself!" << std::endl;
#ifndef RELEASE
        DumpCallStack();
#endif
        throw std::exception();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::~DistMatrix()
{ }

template<typename T>
inline const Elemental::Grid&
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::GetGrid() const
{ return *_grid; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::Viewing() const
{ return _viewing; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::Height() const
{ return _height; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::Width() const
{ return _width; }

template<typename T>
inline T&
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::LocalEntry
( const int i, const int j )
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LocalEntry(i,j)");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( _viewing && _lockedView )
    {
        std::cerr << "Cannot alter data with locked view." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    T& value = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline T
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::LocalEntry
( const int i, const int j ) const
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LocalEntry(i,j)");
    if( i < 0 || j < 0 )
    {
        std::cerr << "Indices must be non-negative." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    T value = _localMatrix(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return value;
}

template<typename T>
inline Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::LocalMatrix()
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,VR]::LocalMatrix");
    if( _viewing && _lockedView )
    {
        std::cerr << "Cannot alter data with locked view." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    PopCallStack();
#endif
    return _localMatrix;
}

template<typename T>
inline const Elemental::Matrix<T>&
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::LockedLocalMatrix() 
const
{ return _localMatrix; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::LocalHeight() const
{ return _localMatrix.Height(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::LocalWidth() const
{ return _localMatrix.Width(); }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::LocalLDim() const
{ return _localMatrix.LDim(); }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::ConstrainedColDist() 
const
{ return false; }

template<typename T>
inline bool
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::ConstrainedRowDist() 
const
{ return _constrainedRowDist; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::ColAlignment() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::RowAlignment() const
{ return _rowAlignment; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::ColShift() const
{ return 0; }

template<typename T>
inline int
Elemental::DistMatrix<T,Elemental::Star,Elemental::VR>::RowShift() const
{ return _rowShift; }

#endif /* ELEMENTAL_DISTMATRIX_STAR_VR_H */
