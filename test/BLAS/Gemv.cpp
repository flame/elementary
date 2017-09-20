/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "Elemental.h"
#include "ElementalBLAS_Internal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

template<typename T, Distribution ColDist, Distribution RowDist>
void Uniform( DistMatrix<T, ColDist, RowDist> &A, int m, int n) {
    A.ResizeTo( m, n );
    A.SetToRandom();
}

void Usage()
{
    cout << "GEneral Matrix Vector multiplication." << endl << endl;;
    cout << "  Gemv <r> <c> <m> <n> <print?> "  
         << endl << endl;
    cout << "  r: number of process rows    " << endl;
    cout << "  c: number of process cols    " << endl;
    cout << "  m: height of C               " << endl;
    cout << "  n: width  of C               " << endl;
    cout << "  print?: [0/1]                " << endl;
    cout << endl;
}

int
main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if ( argc != 6 )
    {
        if ( rank == 0 )
            Usage();
        Elemental::Finalize;
        return 0;
    }
    try 
    {
        const int  r = atoi( argv[1] );
        const int  c = atoi( argv[2] );
        const int  m = atoi( argv[3] );
        const int  n = atoi( argv[4] );
        //const int  nb = atoi( argv[8] );
        const bool print = atoi( argv[5] );

        Barrier( MPI_COMM_WORLD );
        if( rank == 0 ) {
            cout << "---------------------------------------------" << endl;
            cout << "Constructing grid..." << endl;
        }
        const Grid grid( MPI_COMM_WORLD, r, c );
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 )
            cout << "Done building grid." << endl;
        //SetBlocksize( nb );
        Barrier( MPI_COMM_WORLD );
        if( rank == 0 )
            cout << "Done setting blocksize." << endl;

        //const Orientation orientation = ( adjoint ? ADJOINT : NORMAL );
        const Orientation orientation = Normal;

        if( rank == 0 )
        {
            cout << "Will test Gemv" << OrientationToChar(orientation) 
                                     << endl;
        }

        if( rank == 0 )
        {
            cout << "--------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "--------------------" << endl;
        }

        DistMatrix<double, MC, MR> A(grid);
        Uniform( A, m, n );

        // Draw the entries of the original x and y from uniform distributions 
        // over the complex unit ball
        //DistMatrix<double,MC,MR> x(grid);
        //DistMatrix<double,MC,MR> y(grid);
        DistMatrix<double,VR,Star> x(grid);
        DistMatrix<double,VC,Star> y(grid);
 
        if( orientation == Normal )
        {
            Uniform( x, n, 1 );
            Uniform( y, m, 1 );
        }
        else
        {
            Uniform( x, m, 1 );
            Uniform( y, n, 1 );
        }
        if( print )
        {
            A.Print( "A" );
            x.Print( "x" );
            y.Print( "y" );
        }

        // Run the matrix-vector product
        BLAS::Gemv( orientation, (double)(1.0), A, x, (double)(1.0), y );

        if( print )
        {
            if( orientation == Normal )
                y.Print( "y := 1 A x + 1 y" );
            else
                y.Print( "y := 1 A^H x + 1 y" );
        }
    }
    catch( exception& e ) {
        cerr << "Caught exception on process " << rank << endl;
    }

    Elemental::Finalize();
    return 0;
}
