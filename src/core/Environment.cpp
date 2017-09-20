/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "ElementalEnvironment.h"
using namespace std;
using namespace Elemental;

static bool elementalInitializedMPI;
static stack<int> blocksizeStack;

void
Elemental::Init
( int* argc, char** argv[] )
{
#ifndef RELEASE
    PushCallStack("Init");
#endif
    int initialized;
    MPI_Initialized( &initialized );
    if( initialized == 0 )
    {
        int finalized; 
        MPI_Finalized( &finalized );
        if( finalized != 0 )
        {
            cerr << "Cannot initialize Elemental after MPI_Finalize." << endl;
#ifndef RELEASE
            DumpCallStack();
#endif
            throw exception();
        }
        MPI_Init( argc, argv );
        ::elementalInitializedMPI = true;
    }
    else
    {
        ::elementalInitializedMPI = false;
    }

    // Seed the random number generator with out rank
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    srand( rank );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
Elemental::Finalize()
{
#ifndef RELEASE
    PushCallStack("Finalize");
#endif
    int finalized;
    MPI_Finalized( &finalized );
    if( finalized != 0 )
        cerr << "Warning: MPI was finalized before Elemental." << endl;

    if( ::elementalInitializedMPI && finalized == 0 )
        MPI_Finalize();
#ifndef RELEASE
    PopCallStack();
#endif
}

int 
Elemental::Blocksize()
{ 
    if( ::blocksizeStack.size() == 0 )
        ::blocksizeStack.push( 192 );
    return ::blocksizeStack.top(); 
}

void
Elemental::SetBlocksize( int blocksize )
{ 
    if( ::blocksizeStack.size() == 0 )
        ::blocksizeStack.push( blocksize );
    else
        ::blocksizeStack.top() = blocksize; 
}

void
Elemental::PushBlocksizeStack( int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void
Elemental::PopBlocksizeStack()
{ ::blocksizeStack.pop(); }

// If we are not in RELEASE mode, then implement wrappers for a CallStack
#ifndef RELEASE
static stack<string> callStack;

void
Elemental::PushCallStack( string name )
{ ::callStack.push( name ); }

void
Elemental::PopCallStack()
{ ::callStack.pop(); }

void
Elemental::DumpCallStack()
{
    while( ! ::callStack.empty() )
    {
        cerr << "[" << ::callStack.size() << "]: " << ::callStack.top() << endl;
        ::callStack.pop();
    }
}
#endif

