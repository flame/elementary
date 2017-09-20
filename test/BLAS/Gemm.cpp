/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include <cmath>
#include <ctime>
#include <sstream>
#include "Elemental.h"
#include "ElementalBLAS_Internal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

void Usage()
{
    cout << "GEneral Matrix Matrix multiplication." << endl << endl;;
    cout << "  Gemm <r> <c> <orient. of A?> <orient. of B?> <m> <n> <k> <nb> "  
         << endl
         << "       <serial gemm?> <parallel gemm?> <correctness?> <print?> " 
         << endl << endl;
    cout << "  r: number of process rows    " << endl;
    cout << "  c: number of process cols    " << endl;
    cout << "  orient. of A: {N,T,C}        " << endl;
    cout << "  orient. of B: {N,T,C}        " << endl;
    cout << "  m: height of C               " << endl;
    cout << "  n: width  of C               " << endl;
    cout << "  k: inner dimension of AB     " << endl;
    cout << "  nb: algorithmic blocksize    " << endl;
    cout << "  serial gemm? [0/1]           " << endl;
    cout << "  parallel gemm? [0/1]         " << endl;
    cout << "  correctness?: [0/1]          " << endl;
    cout << "  print?: [0/1]                " << endl;
    cout << endl;
}

template<typename T>
bool OKRelativeError( T truth, T computed );

template<>
bool OKRelativeError( float truth, float computed )
{
    return ( fabs(truth-computed) / max(fabs(truth),(float)1)  <= 1e-5 );
}

template<>
bool OKRelativeError( double truth, double computed )
{
    return ( fabs(truth-computed) / max(fabs(truth),(double)1)  <= 1e-13 );
}

//#ifndef WITHOUT_COMPLEX
//template<>
//bool OKRelativeError( scomplex truth, scomplex computed )
//{
//    return ( norm(truth-computed) / max(norm(truth),(float)1) <= 1e-5 );
//}
//
//template<>
//bool OKRelativeError( dcomplex truth, dcomplex computed )
//{
//    return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-13 );
//}
//#endif

template<typename T>
void ManualGemm
( const Orientation orientationOfA, const Orientation orientationOfB,
  const T alpha, const Matrix<T>& A, const Matrix<T>& B,
  const T beta,        Matrix<T>& C                                  )
{
    const int m = C.Height();
    const int n = C.Width();
    BLAS::Scal( beta, C );
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        const int k = A.Width();
        for( int j=0; j<n; ++j )
            for( int i=0; i<m; ++i )
                for( int l=0; l<k; ++l )
                    C(i,j) += alpha * A(i,l) * B(l,j);
    }
    //else if( orientationOfA == Normal )
    //{
    //    const int k = A.Width();
    //    if( orientationOfB == Transpose )
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * A(i,l) * B(j,l);
    //    }
    //    else
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * A(i,l) * BLAS::Conj(B(j,l));
    //    }
    //}
    //else if( orientationOfB == Normal )
    //{
    //    const int k = A.Height();
    //    if( orientationOfA == Transpose )
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * A(l,i) * B(l,j);
    //    }
    //    else
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * BLAS::Conj(A(l,i)) * B(l,j);
    //    }
    //}
    //else
    //{
    //    const int k = A.Height();
    //    if( orientationOfA == Transpose && orientationOfB == Transpose )
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * A(l,i) * B(j,l);
    //    }
    //    else if( orientationOfA == Transpose )
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * A(l,i) * BLAS::Conj(B(j,l));
    //    }
    //    else if( orientationOfB == Transpose )
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * BLAS::Conj(A(l,i)) * B(j,l);
    //    }
    //    else
    //    {
    //        for( int j=0; j<n; ++j )
    //            for( int i=0; i<m; ++i )
    //                for( int l=0; l<k; ++l )
    //                    C(i,j) += alpha * BLAS::Conj(A(l,i)) 
    //                                    * BLAS::Conj(B(j,l));
    //    }
    //}
}

template<typename T>
void TestSerialCorrectness
( const Orientation orientationOfA, const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,Star,Star>& A_ref,
                 const DistMatrix<T,Star,Star>& B_ref,
  const T beta,        DistMatrix<T,Star,Star>& C_ref,
  const DistMatrix<T,Star,Star>& C, const bool printMatrices )
{
    const Grid& grid = C.GetGrid();

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    ManualGemm( orientationOfA, orientationOfB,
                alpha, A_ref.LockedLocalMatrix(),
                       B_ref.LockedLocalMatrix(),
                beta,  C_ref.LocalMatrix()       );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        C_ref.Print("Truth");

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<C.Width(); ++j )
    {
        for( int i=0; i<C.Height(); ++i )
        {
            T truth = C_ref.LocalEntry(i,j);
            T computed = C.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                cout << "FAILED on process " << grid.VCRank() 
                     << " at index (" << i << "," << j << "): truth="
                     << truth << ", computed=" << computed << endl;
                throw exception();
            }
        }
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestParallelCorrectness
( const Orientation orientationOfA, const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,Star,Star>& A_ref,
                 const DistMatrix<T,Star,Star>& B_ref,
  const T beta,        DistMatrix<T,Star,Star>& C_ref,
  const DistMatrix<T,MC,MR>& C, const bool printMatrices )
{
    const Grid& grid = C.GetGrid();
    DistMatrix<T,Star,Star> C_copy(grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    C_copy = C;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    BLAS::Gemm( orientationOfA, orientationOfB,
                alpha, A_ref.LockedLocalMatrix(),
                       B_ref.LockedLocalMatrix(),
                beta,  C_ref.LocalMatrix()       );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        C_ref.Print("Truth");

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    for( int j=0; j<C.Width(); ++j )
    {
        for( int i=0; i<C.Height(); ++i )
        {
            T truth = C_ref.LocalEntry(i,j);
            T computed = C_copy.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                cout << "FAILED on process " << grid.VCRank() 
                     << " at index (" << i << "," << j << "): truth="
                     << truth << ", computed=" << computed << endl;
                throw exception();
            }
        }
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestSerialGemm
( const Orientation orientationOfA, const Orientation orientationOfB,
  const int m, const int n, const int k, const T alpha, const T beta,
  const bool testCorrectness, const bool printMatrices, const Grid& grid )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,Star,Star> A(grid);
    DistMatrix<T,Star,Star> B(grid);
    DistMatrix<T,Star,Star> C(grid);
    DistMatrix<T,Star,Star> A_ref(grid);
    DistMatrix<T,Star,Star> B_ref(grid);
    DistMatrix<T,Star,Star> C_ref(grid);

    if( orientationOfA == Normal )
        A.ResizeTo( m, k );
    else
        A.ResizeTo( k, m );

    if( orientationOfB == Normal )
        B.ResizeTo( k, n );
    else
        B.ResizeTo( n, k );

    C.ResizeTo( m, n );

    if( grid.VCRank() == 0 )
        cout << "Serial Gemm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        A_ref = A;
        B_ref = B;
        C_ref = C;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( grid.VCRank() == 0 )
    {
        cout << "  Starting Gemm...";
        cout.flush();
    }
    Barrier( grid.VCComm() );
    startTime = Time();
    BLAS::Gemm
    ( orientationOfA, orientationOfB, 
      alpha, A.LockedLocalMatrix(), B.LockedLocalMatrix(), 
      beta,  C.LocalMatrix()                              );
    Barrier( grid.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = BLAS::Internal::GemmGFlops<T>(m,n,k,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestSerialCorrectness
        ( orientationOfA, orientationOfB, 
          alpha, A_ref, B_ref, beta, C_ref, C, printMatrices );
    }
}

template<typename T>
void TestParallelGemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const int m, const int n, const int k, 
  const T alpha, const T beta, 
  const bool testCorrectness, 
  const bool printMatrices, 
  const Grid& grid                      )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,  MR  > A(grid);
    DistMatrix<T,MC,  MR  > B(grid);
    DistMatrix<T,MC,  MR  > C(grid);
    DistMatrix<T,Star,Star> A_ref(grid);
    DistMatrix<T,Star,Star> B_ref(grid);
    DistMatrix<T,Star,Star> C_ref(grid);

    if( orientationOfA == Normal )
        A.ResizeTo( m, k );
    else
        A.ResizeTo( k, m );

    if( orientationOfB == Normal )
        B.ResizeTo( k, n );
    else
        B.ResizeTo( n, k );

    C.ResizeTo( m, n );

    // Test the variant of Gemm that keeps A stationary
    if( grid.VCRank() == 0 )
        cout << "Stationary A Algorithm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        A_ref = A;
        B_ref = B;
        C_ref = C;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( grid.VCRank() == 0 )
    {
        cout << "  Starting Parallel Gemm...";
        cout.flush();
    }
    Barrier( grid.VCComm() );
    startTime = Time();
    BLAS::Internal::Gemm_A
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    Barrier( grid.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = BLAS::Internal::GemmGFlops<T>(m,n,k,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestParallelCorrectness
        ( orientationOfA, orientationOfB, 
          alpha, A_ref, B_ref, beta, C_ref, C, printMatrices );
    }

    // Test the variant of Gemm that keeps B stationary
    if( grid.VCRank() == 0 )
        cout << endl << "Stationary B Algorithm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        A_ref = A;
        B_ref = B;
        C_ref = C;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( grid.VCRank() == 0 )
    {
        cout << "  Starting Parallel Gemm...";
        cout.flush();
    }
    Barrier( grid.VCComm() );
    startTime = Time();
    BLAS::Internal::Gemm_B
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    Barrier( grid.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = BLAS::Internal::GemmGFlops<T>(m,n,k,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestParallelCorrectness
        ( orientationOfA, orientationOfB,
          alpha, A_ref, B_ref, beta, C_ref, C, printMatrices );
    }

    // Test the variant of Gemm that keeps C stationary
    if( grid.VCRank() == 0 )
        cout << endl << "Stationary C Algorithm:" << endl;
    A.SetToRandom();
    B.SetToRandom();
    C.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        A_ref = A;
        B_ref = B;
        C_ref = C;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
        B.Print("B");
        C.Print("C");
    }
    if( grid.VCRank() == 0 )
    {
        cout << "  Starting Parallel Gemm...";
        cout.flush();
    }
    Barrier( grid.VCComm() );
    startTime = Time();
    BLAS::Internal::Gemm_C
    ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    Barrier( grid.VCComm() );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = BLAS::Internal::GemmGFlops<T>(m,n,k,runTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        ostringstream msg;
        msg << "C := " << alpha << " A B + " << beta << " C";
        C.Print( msg.str() );
    }
    if( testCorrectness )
    {
        TestParallelCorrectness
        ( orientationOfA, orientationOfB,
          alpha, A_ref, B_ref, beta, C_ref, C, printMatrices );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 13 )
    {
        if( rank == 0 )
            Usage();
        Elemental::Finalize();
        return 0;
    }
    try
    {
        const int         r = atoi( argv[1] );
        const int         c = atoi( argv[2] );
        const Orientation orientationOfA = CharToOrientation( *argv[3] );
        const Orientation orientationOfB = CharToOrientation( *argv[4] );
        const int         m = atoi( argv[5] );
        const int         n = atoi( argv[6] );
        const int         k = atoi( argv[7] );
        const int         nb = atoi( argv[8] );
        const bool        testSerial = atoi( argv[9] );
        const bool        testParallel = atoi( argv[10] );
        const bool        testCorrectness = atoi( argv[11] );
        const bool        printMatrices = atoi( argv[12] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 ) {
            cout << "---------------------------------------------" << endl;
            cout << "Constructing grid..." << endl;
        }
        const Grid grid( MPI_COMM_WORLD, r, c );
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 )
            cout << "Done building grid." << endl;
        SetBlocksize( nb );
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 )
            cout << "Done setting blocksize." << endl;

        if( rank == 0 )
        {
            cout << "Will test Gemm" << OrientationToChar(orientationOfA) 
                                     << OrientationToChar(orientationOfB) 
                                     << endl;
        }

        if( rank == 0 )
        {
            cout << "--------------------" << endl;
            cout << "Testing with floats:" << endl;
            cout << "--------------------" << endl;
        }
        if( testSerial )
        {
            TestSerialGemm<float>
            ( orientationOfA, orientationOfB,
              m, n, k, (float)3, (float)4, 
              testCorrectness, printMatrices, grid );
        }
        if( testParallel )
        {
            TestParallelGemm<float>
            ( orientationOfA, orientationOfB,
              m, n, k, (float)3, (float)4,
              testCorrectness, printMatrices, grid );
        }
        if( rank == 0 )
            cout << endl;

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        if( testSerial )
        {
            TestSerialGemm<double>
            ( orientationOfA, orientationOfB,
              m, n, k, (double)3, (double)4,
              testCorrectness, printMatrices, grid );
        }
        if( testParallel )
        {
            TestParallelGemm<double>
            ( orientationOfA, orientationOfB,
              m, n, k, (double)3, (double)4,
              testCorrectness, printMatrices, grid );
        }
        if( rank == 0 )
            cout << endl;
    }
    catch( exception e )
    {
        cerr << "Caught exception on process " << rank << endl;
    }
    Elemental::Finalize();
    return 0;
}

