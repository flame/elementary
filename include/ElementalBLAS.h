/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef ELEMENTAL_BLAS_H
#define ELEMENTAL_BLAS_H 1

#include <cmath>
#include "ElementalDistMatrix.h"
#include "ElementalPartitioning.h"
#include "wrappers/BLAS.h"

namespace Elemental 
{
    namespace BLAS 
    {
        //--------------------------------------------------------------------//
        // Local BLAS: Level 0 (extensions)                                   //
        //--------------------------------------------------------------------//

        template<typename R>
        R 
        Conj( const R& alpha );

        //--------------------------------------------------------------------//
        // Local BLAS: Level 1                                                //
        //--------------------------------------------------------------------//

        // AXPY: Y := Alpha X Plus Y 
        template<typename T>
        void
        Axpy( const T alpha, const Matrix<T>& X, Matrix<T>& Y );

        // COPY: Copy
        template<typename T>
        void
        Copy( const Matrix<T>& X, Matrix<T>& Y );


        // SCAL: SCALe X by alpha
        template<typename T>
        void
        Scal( const T alpha, Matrix<T>& X );
        
        //--------------------------------------------------------------------//
        // Local BLAS: Level 1 (extensions)                                   //
        //--------------------------------------------------------------------//

        // CONJ: CONJugate in-place
        //
        // There are two implementations because, for real datatypes, Conj is
        // a no-op. Partial specialization of function templates is not allowed,
        // so we must have two declarations.
        template<typename R>
        void
        Conj( Matrix<R>& A );


        // CONJ: CONJugated copy
        template<typename T>
        void
        Conj( const Matrix<T>& A, Matrix<T>& B );

        // CONJTRANS: CONJugated Transposed copy
        template<typename T>
        void
        ConjTrans( const Matrix<T>& A, Matrix<T>& B );

        // TRANS: TRANSposed copy
        template<typename T>
        void
        Trans( const Matrix<T>& A, Matrix<T>& B );

        //--------------------------------------------------------------------//
        // Local BLAS: Level 2                                                //
        //--------------------------------------------------------------------//

        // GEMV: GEneral Matrix-Vector multiply
        template<typename T>
        void
        Gemv( const Orientation orientation,
              const T alpha, const Matrix<T>& A, const Matrix<T>& x,
              const T beta,        Matrix<T>& y                     );

        // GER: GEneral Rank-one update
        //
        // For complex datatypes it routes to Gerc, as x (tensor product) y 
        // is x * conj(y)^T. That is, the dual of y is its conjugate transpose
        // thanks to the Riesz map. Thus our generalized Ger is equivalent to
        // our generalized Gerc.
        template<typename T>
        void
        Ger( const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                  Matrix<T>& A                     );

        // GERC: GEneral Rank-one Conjugated update
        template<typename T>
        void
        Gerc( const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                   Matrix<T>& A                     );

        // GERU: GEneral Rank-one Unconjugated update
        template<typename T>
        void
        Geru( const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                   Matrix<T>& A                     );

        //--------------------------------------------------------------------//
        // Local BLAS: Level 3                                                //
        //--------------------------------------------------------------------//

        // GEMM: GEneral Matrix-Matrix multiplication
        template<typename T>
        void
        Gemm
        ( const Orientation orientationOfA, const Orientation orientationOfB,
          const T alpha, const Matrix<T>& A, const Matrix<T>& B,
          const T beta,        Matrix<T>& C                                  );

        
        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 1                                          //
        //--------------------------------------------------------------------//

        // AXPY: Y := Alpha X Plus Y 
        template<typename T, Distribution U, Distribution V>
        void
        Axpy( const T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );

        // COPY: Copy
        //
        // In our case, it is just a wrapper around the '=' operator for those
        // that prefer BLAS/PLAPACK syntax.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        void
        Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

        // DOT: alpha := conj(x)^T * y
        // 
        // Though the standard BLAS interface only defines DOT for real 
        // datatypes, it is naturally generalized to an inner product over the
        // complex field. Recall that the conjugate symmetry of inner products 
        // requires that (x,y) = conj(y,x), so that (x,x) = conj( (x,x) ) => 
        // (x,x) is real. This requires that we choose (x,x) = conj(x)^T * x.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        T
        Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

        // DOTC: alpha := conj(x)^T * y
        //
        // This is the sister routine to DOT; while DOT is originally defined 
        // only over the reals, DOTC was defined only over the complex field. 
        // They are each others' extensions, and so, to us, they are 
        // identical.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        T
        Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

        // DOTU: alpha := x^T * y
        //
        // Standard BLAS defines DOTU for complex datatypes, but the operation
        // is perfectly valid over the reals (clearly), so we extend it.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        T
        Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

        // NRM2: NoRM 2 (Euclidean norm)
        template<typename R>
        R
        Nrm2( const DistMatrix<R,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        R
        Nrm2( const DistMatrix< std::complex<R>, MC, MR >& x );
#endif

        // SCAL: SCALe by a constant
        template<typename T, Distribution U, Distribution V>
        void
        Scal
        ( const T alpha, DistMatrix<T,U,V>& A );

        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 1 (extensions)                             //
        //--------------------------------------------------------------------//


        // TRANS: TRANSposed copy
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        void
        Trans( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 2                                          //
        //--------------------------------------------------------------------//

        // GEMV: GEneral Matrix-Vector multiplication
        template<typename T>
        void
        Gemv
        ( const Orientation orientationOfA,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& x,
          const T beta,        DistMatrix<T,MC,MR>& y );


        template<typename T>
        void
        Gemv
        ( const Orientation orientationOfA,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,VR,Star>& x,
          const T beta,        DistMatrix<T,VC,Star>& y );

        // GER: GEneral Rank-one update
        //
        // For complex datatypes it routes to Gerc, as x (tensor product) y 
        // is x * conj(y)^T. That is, the dual of y is its conjugate transpose
        // thanks to the Riesz map. Thus our generalized Ger is equivalent to
        // our generalized Gerc.
        template<typename T>
        void
        Ger
        ( const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );

        template<typename T>
        void
        Ger
        ( const T alpha, const DistMatrix<T,VR,Star>& x,
                         const DistMatrix<T,VC,Star>& y,
                               DistMatrix<T,MC,MR>& A );

        // GERC: GEneral Rank-one Conjugated update
        //
        // Since the extension of Ger to complex datatypes 
        template<typename T>
        void
        Gerc
        ( const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );
        
        // GERU: GEneral Rank-one Unconjugated update
        template<typename T>
        void
        Geru
        ( const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );

        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 3                                          //
        //--------------------------------------------------------------------//

        // GEMM: GEneral Matrix-Matrix multiplication
        template<typename T>
        void
        Gemm
        ( const Orientation orientationOfA, 
          const Orientation orientationOfB,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& B,
          const T beta,        DistMatrix<T,MC,MR>& C );

    }
}

/*----------------------------------------------------------------------------*/

//----------------------------------------------------------------------------//
// Local BLAS: Level 0 (extensions)                                           //
//----------------------------------------------------------------------------//


template<typename R>
inline R
Elemental::BLAS::Conj
( const R& alpha )
{ return alpha; }

//----------------------------------------------------------------------------//
// Local BLAS: Level 1                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Elemental::BLAS::Axpy
( const T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Axpy");
    if( X.Height() != Y.Height() || X.Width() != Y.Width() )
    {
        std::cerr << "Nonconformal Axpy." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    if( X.Width() <= X.Height() )
    {
        for( int j=0; j<X.Width(); ++j )
        {
            Elemental::wrappers::BLAS::Axpy
            ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
        }
    }
    else
    {
        for( int i=0; i<X.Height(); ++i )
        {
            Elemental::wrappers::BLAS::Axpy
            ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                Y.Buffer(i,0),       Y.LDim() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Copy
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Scal
( const T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Scal");
#endif
    if( alpha != (T)1 )
    {
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                Elemental::wrappers::BLAS::Scal
                ( X.Height(), alpha, X.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                Elemental::wrappers::BLAS::Scal
                ( X.Width(), alpha, X.Buffer(i,0), X.LDim() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 1 (extensions)                                           //
//----------------------------------------------------------------------------//

// Default case is for real datatypes
template<typename R>
inline void
Elemental::BLAS::Conj
( Matrix<R>& A )
{ }

#ifndef WITHOUT_COMPLEX
// Specialization is to complex datatypes
template<typename R>
inline void
Elemental::BLAS::Conj
( Matrix< std::complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Conj (in-place)");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A(i,j) = BLAS::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template<typename T>
inline void
Elemental::BLAS::Conj
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Conj");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(i,j) = BLAS::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::ConjTrans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::ConjTrans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( n, m );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(j,i) = BLAS::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Trans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( n, m );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(j,i) = A(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 2                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Elemental::BLAS::Gemv
( const Orientation orientation,
  const T alpha, const Matrix<T>& A, const Matrix<T>& x,
  const T beta,        Matrix<T>& y                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == Normal )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::cerr << "A must conform with x and y:" << 
            std::endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            std::endl << "  x ~ " << x.Height() << " x " << x.Width() <<
            std::endl << "  y ~ " << y.Height() << " x " << y.Width() << 
            std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::cerr << "A must conform with x and y:" << 
            std::endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            std::endl << "  x ~ " << x.Height() << " x " << x.Width() <<
            std::endl << "  y ~ " << y.Height() << " x " << y.Width() << 
            std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
#endif
    const char transChar = OrientationToChar( orientation );
    const int m = A.Height();
    const int n = A.Width();
    const int k = ( transChar == 'N' ? n : m );
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( k != 0 )
    {
        Elemental::wrappers::BLAS::Gemv
        ( transChar, m, n, 
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
          beta,  y.Buffer(), incy                                   );
    }
    else
    {
        Elemental::BLAS::Scal( beta, y );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Ger
( const T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Ger");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::cerr << "Nonconformal Ger: " << std::endl;
        std::cerr << "  x ~ " << x.Height() << " x " << x.Width() << std::endl;
        std::cerr << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
        std::cerr << "  A ~ " << A.Height() << " x " << A.Width() << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim()                           );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 3                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Elemental::BLAS::Gemm
( const Orientation orientationOfA, const Orientation orientationOfB,
  const T alpha, const Matrix<T>& A, const Matrix<T>& B,
  const T beta,        Matrix<T>& C                                  )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemm");
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height()    )
        {
            std::cerr << "Nonconformal GemmNN." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientationOfA == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width()     )
        {
            std::cerr << "Nonconformal GemmN(T/C)." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientationOfB == Normal )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height()    )
        {
            std::cerr << "Nonconformal Gemm(T/C)N." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width()     )
        {
            std::cerr << "Nonconformal Gemm(T/C)(T/C)." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == Normal ? A.Width() : A.Height() );
    if( k != 0 )
    {
        Elemental::wrappers::BLAS::Gemm
        ( transA, transB, m, n, k, 
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        Elemental::BLAS::Scal( beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}


//----------------------------------------------------------------------------//
// Distributed BLAS: Level 1                                                  //
//----------------------------------------------------------------------------//

template<typename T, Elemental::Distribution U, Elemental::Distribution V>
inline void
Elemental::BLAS::Axpy
( const T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Axpy");
    const Grid& grid = X.GetGrid();
    if( X.GetGrid() != Y.GetGrid() )
    {
        if( grid.VCRank() == 0 ) 
        {
            std::cerr << "X and Y must be distributed over the same grid."
                      << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
    if( X.ColAlignment() != Y.ColAlignment() ||
        X.RowAlignment() != Y.RowAlignment()   )
    {
        if( grid.VCRank() == 0 )
            std::cerr << "Axpy requires X and Y be aligned." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    BLAS::Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V,
                     Elemental::Distribution W, Elemental::Distribution Z >
inline void
Elemental::BLAS::Copy
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V>
inline void
Elemental::BLAS::Scal
( const T alpha, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Scal");
#endif
    BLAS::Scal( alpha, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V,
                     Elemental::Distribution W, Elemental::Distribution Z >
inline void
Elemental::BLAS::Trans
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trans");
#endif
    DistMatrix<T,Z,W> C( B.GetGrid() );
    if( B.ConstrainedColDist() )
        C.AlignRowsWith( B );
    if( B.ConstrainedRowDist() )
        C.AlignColsWith( B );

    C = A;

    if( !B.ConstrainedColDist() )
        B.AlignColsWith( C );
    if( !B.ConstrainedRowDist() )
        B.AlignRowsWith( C );

    B.ResizeTo( A.Width(), A.Height() );
    BLAS::Trans( C.LockedLocalMatrix(), B.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif /* ELEMENTAL_BLAS_H */

