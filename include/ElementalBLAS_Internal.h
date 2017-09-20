/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef ELEMENTAL_BLAS_INTERNAL_H
#define ELEMENTAL_BLAS_INTERNAL_H 1

#include "ElementalBLAS.h"

namespace Elemental
{
    namespace BLAS
    {
        namespace Internal
        {

            //----------------------------------------------------------------//
            // Distributed BLAS: Level 2                                      //
            //----------------------------------------------------------------//
            
            // Gemv where A is not transposed
            template<typename T>
            void
            GemvN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& x,
              const T beta,        DistMatrix<T,MC,MR>& y );

            template<typename T>
            void
            GemvN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,VR,Star>& x,
              const T beta,        DistMatrix<T,VC,Star>& y );

            //----------------------------------------------------------------//
            // Distributed BLAS: Level 3                                      //
            //----------------------------------------------------------------//

            // Gemm where we avoid redistributing A.
            template<typename T>
            void
            Gemm_A
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, 
              const DistMatrix<T,MC,MR>& A,
              const DistMatrix<T,MC,MR>& B,
              const T beta,
                    DistMatrix<T,MC,MR>& C    );

            // Gemm where we avoid redistributing B.
            template<typename T>
            void
            Gemm_B
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, 
              const DistMatrix<T,MC,MR>& A,
              const DistMatrix<T,MC,MR>& B,
              const T beta, 
                    DistMatrix<T,MC,MR>& C     );

            // Gemm where we avoid redistributing C.
            template<typename T>
            void
            Gemm_C
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, 
              const DistMatrix<T,MC,MR>& A,
              const DistMatrix<T,MC,MR>& B,
              const T beta,
                    DistMatrix<T,MC,MR>& C     );

            // Normal Normal Gemm.
            template<typename T>
            void
            GemmNN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal Normal Gemm where we avoid redistributing A.
            template<typename T>
            void
            GemmNN_A
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal Normal Gemm where we avoid redistributing B.
            template<typename T>
            void
            GemmNN_B
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal Normal Gemm where we avoid redistributing C.
            template<typename T>
            void
            GemmNN_C
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );


            //----------------------------------------------------------------//
            // Level 3 BLAS Utility Functions                                 //
            //----------------------------------------------------------------//
            template<typename T>
            double 
            GemmGFlops
            ( const int m, const int n, const int k, const double seconds );


        }
    }
}

/*----------------------------------------------------------------------------*/

namespace Elemental
{
    namespace BLAS
    {
        namespace Internal
        {
            // Level 3 Utility functions
            template<>
            inline double
            GemmGFlops<float>
            ( const int m, const int n, const int k, const double seconds )
            { return (2.*m*n*k)/(1.e9*seconds); }

            template<>
            inline double
            GemmGFlops<double>
            ( const int m, const int n, const int k, const double seconds )
            { return GemmGFlops<float>(m,n,k,seconds); }

            
        }
    }
}

#endif /* ELEMENTAL_BLAS_INTERNAL_H */

