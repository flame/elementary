/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef ELEMENTAL_DISTMATRIX_H
#define ELEMENTAL_DISTMATRIX_H 1

#include "ElementalMatrix.h"
#include "wrappers/MPI.h"

namespace Elemental
{
    // We will partially specialize for each valid distribution 
    template<typename T, Distribution ColDist, Distribution RowDist> 
    class DistMatrix;
}

#include "ElementalDistMatrix_MC_MR.h"
#include "ElementalDistMatrix_MC_Star.h"
#include "ElementalDistMatrix_MR_MC.h"
#include "ElementalDistMatrix_MR_Star.h"
#include "ElementalDistMatrix_Star_MC.h"
#include "ElementalDistMatrix_Star_MR.h"
#include "ElementalDistMatrix_Star_Star.h"
#include "ElementalDistMatrix_Star_VC.h"
#include "ElementalDistMatrix_Star_VR.h"
#include "ElementalDistMatrix_VC_Star.h"
#include "ElementalDistMatrix_VR_Star.h"

#endif /* ELEMENTAL_DISTMATRIX_H */
