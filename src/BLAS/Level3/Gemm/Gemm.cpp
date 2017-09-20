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
Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemm");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN( alpha, A, B, beta, C );
    }
    else
    {
        // Not supported
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Gemm_A");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN_A( alpha, A, B, beta, C );
    }
    else {
        // Not supported
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Gemm_B");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN_B( alpha, A, B, beta, C );
    }
    else {
        // Not supported
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Gemm_C");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN_C( alpha, A, B, beta, C );
    }
    else {
        // Not supported
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );


