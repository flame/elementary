/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "ElementalBLAS_Internal.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Internal::GemmNN
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNN");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Width();
    const float weightTowardsC = 0.5;

    if( m <= n && m <= weightTowardsC*k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNN routing to GemmNN_B." << endl;
#endif
        BLAS::Internal::GemmNN_B( alpha, A, B, beta, C );    
    }
    else if( n <= m && n <= weightTowardsC*k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNN routing to GemmNN_A." << endl;
#endif
        BLAS::Internal::GemmNN_A( alpha, A, B, beta, C );
    }
    else
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNN routing to GemmNN_C." << endl;
#endif
        BLAS::Internal::GemmNN_C( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Normal Gemm that avoids communicating the matrix A.
template<typename T>
void
Elemental::BLAS::Internal::GemmNN_A
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNN_A");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Nonconformal GemmNN_A: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() << 
            endl << "  B ~ " << B.Height() << " x " << B.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> BL(grid), BR(grid),
                        B0(grid), B1(grid), B2(grid);
    DistMatrix<T,MC,MR> CL(grid), CR(grid),
                        C0(grid), C1(grid), C2(grid);

    // Temporary distributions
    DistMatrix<T,MR,Star> B1_MR_Star(grid);
    DistMatrix<T,MC,Star> D1_MC_Star(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionRight( B, BL, BR );
    PartitionRight( C, CL, CR );
    while( BR.Width() > 0 )
    {
        LockedRepartitionRight( BL, /**/     BR,
                                B0, /**/ B1, B2 );

        RepartitionRight( CL, /**/     CR,
                          C0, /**/ C1, C2 );

        B1_MR_Star.ConformWith( A );
        D1_MC_Star.AlignWith( A );
        D1_MC_Star.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        B1_MR_Star = B1; // B1[MR,*] <- B1[MC,MR]

        // D1[MC,*] := alpha A[MC,MR] B1[MR,*]
        BLAS::Gemm( Normal, Normal, 
                    alpha, A.LockedLocalMatrix(),
                           B1_MR_Star.LockedLocalMatrix(),
                    (T)0,  D1_MC_Star.LocalMatrix()       );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.ReduceScatterUpdate( (T)1, D1_MC_Star );
        //--------------------------------------------------------------------//
        B1_MR_Star.FreeConstraints();
        D1_MC_Star.FreeConstraints();

        SlideLockedPartitionRight( BL,     /**/ BR,
                                   B0, B1, /**/ B2 );

        SlidePartitionRight( CL,     /**/ CR,
                             C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Normal Gemm that avoids communicating the matrix B.
template<typename T>
void 
Elemental::BLAS::Internal::GemmNN_B
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNN_B");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Nonconformal GemmNN_B: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() << 
            endl << "  B ~ " << B.Height() << " x " << B.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(grid),  A0(grid),
                        AB(grid),  A1(grid),
                                   A2(grid);
    DistMatrix<T,MC,MR> CT(grid),  C0(grid),
                        CB(grid),  C1(grid),
                                   C2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MC> A1_Star_MC(grid);
    DistMatrix<T,Star,MR> D1_Star_MR(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionDown( A, AT,
                            AB );
    PartitionDown( C, CT,
                      CB );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown( AT,  A0,
                              /**/ /**/
                                    A1,
                               AB,  A2 );

        RepartitionDown( CT,  C0,
                        /**/ /**/
                              C1,
                         CB,  C2 );

        A1_Star_MC.ConformWith( B );
        D1_Star_MR.AlignWith( B );
        D1_Star_MR.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        A1_Star_MC = A1; // A1[*,MC] <- A1[MC,MR]

        // D1[*,MR] := alpha A1[*,MC] B[MC,MR]
        BLAS::Gemm( Normal, Normal, 
                    alpha, A1_Star_MC.LockedLocalMatrix(),
                           B.LockedLocalMatrix(),
                    (T)0,  D1_Star_MR.LocalMatrix()       );

        // C1[MC,MR] += scattered result of D1[*,MR] summed over grid cols
        C1.ReduceScatterUpdate( (T)1, D1_Star_MR );
        //--------------------------------------------------------------------//
        A1_Star_MC.FreeConstraints();
        D1_Star_MR.FreeConstraints();

        SlideLockedPartitionDown( AT,  A0,
                                       A1,
                                 /**/ /**/
                                  AB,  A2 );
 
        SlidePartitionDown( CT,  C0,
                                 C1,
                           /**/ /**/
                            CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}                     

// Normal Normal Gemm that avoids communicating the matrix C.
template<typename T>
void 
Elemental::BLAS::Internal::GemmNN_C
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNN_C");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Nonconformal GemmNN_C: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() << 
            endl << "  B ~ " << B.Height() << " x " << B.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AL(grid), AR(grid),
                        A0(grid), A1(grid), A2(grid);         

    DistMatrix<T,MC,MR> BT(grid),  B0(grid),
                        BB(grid),  B1(grid),
                                   B2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(grid);
    DistMatrix<T,Star,MR> B1_Star_MR(grid); 

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionRight( A, AL, AR ); 
    LockedPartitionDown( B, BT, 
                            BB ); 
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1, 
                               BB,  B2 );

        A1_MC_Star.AlignWith( C );
        B1_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MC_Star = A1; // A1[MC,*] <- A1[MC,MR]
        B1_Star_MR = B1; // B1[*,MR] <- B1[MC,MR]

        // C[MC,MR] += alpha A1[MC,*] B1[*,MR]
        BLAS::Gemm( Normal, Normal, 
                    alpha, A1_MC_Star.LockedLocalMatrix(),
                           B1_Star_MR.LockedLocalMatrix(),
                    (T)1,  C.LocalMatrix()                );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeConstraints();
        B1_Star_MR.FreeConstraints();

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::GemmNN
( const float alpha, const DistMatrix<float,MC,MR>& A,     
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNN_A
( const float alpha, const DistMatrix<float,MC,MR>& A,     
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNN_B
( const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNN_C
( const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNN
( const double alpha, const DistMatrix<double,MC,MR>& A,         
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNN_A
( const double alpha, const DistMatrix<double,MC,MR>& A,         
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNN_B
( const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNN_C
( const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

