/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef ELEMENTAL_ENVIRONMENT_H 
#define ELEMENTAL_ENVIRONMENT_H 1

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stack>
#include <vector>
#include "mpi.h"

namespace Elemental
{
#ifndef RELEASE
    void PushCallStack( std::string );
    void PopCallStack();
    void DumpCallStack();
#endif
}

#include "ElementalMemory.h"
#include "ElementalGrid.h"
#include "ElementalRandom.h"
#include "ElementalTypes.h"
#include "ElementalUtilities.h"

namespace Elemental
{
    void Init( int* argc, char** argv[] );
    void Finalize();

    int Blocksize();
    void SetBlocksize( int blocksize );

    void PushBlocksizeStack( int blocksize );
    void PopBlocksizeStack();
}

#endif /* ELEMENTAL_ENVIRONMENT_H */

