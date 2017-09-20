/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "ElementalEnvironment.h"
#include "wrappers/MPI.h"

#include <cmath>

using namespace Elemental;
using namespace std;

Elemental::Grid::Grid
( MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid(comm)");
#endif
    // Pull in data from comm
    MPI_Comm_dup  (  comm, &_comm  );
    MPI_Comm_group( _comm, &_group );
    MPI_Comm_size ( _comm, &_p     );
    MPI_Comm_rank ( _comm, &_rank  );

    // Factor p
    int r = static_cast<int>(sqrt(static_cast<double>(_p)));
    while( _p % r != 0 )
        ++r;
    int c = _p / r;

    Init( r, c );

#ifndef RELEASE
    PopCallStack();
#endif
}

// Currently forces a columnMajor absolute rank on the grid
Elemental::Grid::Grid
( MPI_Comm comm, int r, int c )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid( comm, r, c )");
#endif
    
    // Pull in data from comm
    MPI_Comm_dup  (  comm, &_comm  );
    MPI_Comm_group( _comm, &_group );
    MPI_Comm_size ( _comm, &_p     );
    MPI_Comm_rank ( _comm, &_rank  );

    Init( r, c );

#ifndef RELEASE
    PopCallStack();
#endif
}

void
Elemental::Grid::Init
( int r, int c )
{
#ifndef RELEASE
    PushCallStack("Grid::Init(r,c)");
    if( r <= 0 || c <= 0 )
    {
        cerr << "r and c must be positive." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    if( _p != r*c )
    {
        if( _rank == 0 )
        {
            cerr << "Number of processes must match grid size:" << endl;
            cerr << "  p=" << _p << ", (r,c)=(" << r << "," << c << ")" << endl;
        }
#ifndef RELEASE
        DumpCallStack();
#endif
        throw exception();
    }
    _r = r;
    _c = c;

    _gcd = utilities::GCD( r, c );
    int lcm = _p / _gcd;

#ifndef RELEASE
    if( _rank == 0 )
    {
        cout << "Building process grid with:" << endl;   
        cout << "  p=" << _p << ", (r,c)=(" << r << "," << c << ")" << endl;
        cout << "  gcd=" << _gcd << endl;
    }
#endif

    // Set up the MatrixCol communicator
    int* matrixColRanks = new int[r];
    for( int i=0; i<r; ++i )
        matrixColRanks[i] = (_rank / r)*r + i;
    MPI_Group_incl( _group, r, matrixColRanks, &_matrixColGroup );
    MPI_Comm_create( _comm, _matrixColGroup, &_matrixColComm );
    MPI_Comm_rank( _matrixColComm, &_matrixColRank );
    delete[] matrixColRanks;

    // Set up the MatrixRow communicator
    int* matrixRowRanks = new int[c];
    for( int j=0; j<c; ++j )
        matrixRowRanks[j] = (_rank % r) + j*r;
    MPI_Group_incl( _group, c, matrixRowRanks, &_matrixRowGroup );
    MPI_Comm_create( _comm, _matrixRowGroup, &_matrixRowComm );
    MPI_Comm_rank( _matrixRowComm, &_matrixRowRank );
    delete[] matrixRowRanks;

    // Set up the VectorCol communicator
    int* vectorColRanks = new int[_p];
    for( int i=0; i<_p; ++i )
        vectorColRanks[i] = i;
    MPI_Group_incl( _group, _p, vectorColRanks, &_vectorColGroup );
    MPI_Comm_create( _comm, _vectorColGroup, &_vectorColComm );
    MPI_Comm_rank( _vectorColComm, &_vectorColRank );
    delete[] vectorColRanks;

    // Set up the VectorRow communicator
    int* vectorRowRanks = new int[_p];
    for( int i=0; i<r; ++i )
        for( int j=0; j<c; ++j )
            vectorRowRanks[j+i*c] = i+j*r;
    MPI_Group_incl( _group, _p, vectorRowRanks, &_vectorRowGroup );
    MPI_Comm_create( _comm, _vectorRowGroup, &_vectorRowComm );
    MPI_Comm_rank( _vectorRowComm, &_vectorRowRank );
    delete[] vectorRowRanks;

    // Compute which diagonal 'path' we're in, and what our rank is, then
    // perform AllGather world to store everyone's info
    _diagPathsAndRanks = new int[2*_p];
    int* myDiagPathAndRank = new int[2];
    myDiagPathAndRank[0] = (_matrixRowRank+r-_matrixColRank) % _gcd;
    int diagPathRank = 0;
    int row = 0;
    int col = myDiagPathAndRank[0];
    for( int j=0; j<lcm; ++j )
    {
        if( row == _matrixColRank && col == _matrixRowRank )
        {
            myDiagPathAndRank[1] = diagPathRank;
            break;
        }
        else
        {
            row = (row + 1) % r;
            col = (col + 1) % c;
            ++diagPathRank;
        }
    }
    wrappers::MPI::AllGather
    ( myDiagPathAndRank, 2, _diagPathsAndRanks, 2, _vectorColComm );
    delete myDiagPathAndRank;

#ifndef RELEASE
    MPI_Errhandler_set( _matrixColComm, MPI_ERRORS_RETURN );
    MPI_Errhandler_set( _matrixRowComm, MPI_ERRORS_RETURN );
    MPI_Errhandler_set( _vectorColComm, MPI_ERRORS_RETURN );
    MPI_Errhandler_set( _vectorRowComm, MPI_ERRORS_RETURN );
#endif

#ifndef RELEASE
    PopCallStack();
#endif
}

