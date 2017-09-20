/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "ElementalBLAS_Internal.h"
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Gemv
( const Orientation orientation,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& x,
  const T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemv");
#endif
    if( orientation == Normal )
        BLAS::Internal::GemvN( alpha, A, x, beta, y );
    //else
    //    BLAS::Internal::GemvT( orientation, alpha, A, x, beta, y );
    else {
        // Not supported
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Gemv
( const Orientation orientation,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& x,
  const float beta,        DistMatrix<float,MC,MR>& y );

template void Elemental::BLAS::Gemv
( const Orientation orientation,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& x,
  const double beta,        DistMatrix<double,MC,MR>& y );


template<typename T>
void
Elemental::BLAS::Gemv
( const Orientation orientation,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,VR,Star>& x,
  const T beta,        DistMatrix<T,VC,Star>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemv");
#endif
    if( orientation == Normal )
        BLAS::Internal::GemvN( alpha, A, x, beta, y );
    //else
    //    BLAS::Internal::GemvT( orientation, alpha, A, x, beta, y );
    else {
        // Not supported
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Gemv ( const Orientation orientation,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,VR,Star>& x,
  const float beta,        DistMatrix<float,VC,Star>& y );

template void Elemental::BLAS::Gemv
( const Orientation orientation,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,VR,Star>& x,
  const double beta,        DistMatrix<double,VC,Star>& y );

