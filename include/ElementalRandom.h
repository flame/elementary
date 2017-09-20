/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef ELEMENTAL_RANDOM_H
#define ELEMENTAL_RANDOM_H 1

namespace Elemental
{
    // Generate a sample from a uniform PDF over the unit ball about the origin 
    // of the vector space implied by the type T
    template<typename T> T Random();
}

/*----------------------------------------------------------------------------*/

namespace Elemental
{
    template<>
    inline int
    Random<int>()
    {
        int sample = rand();
        if( sample <= RAND_MAX/3 )
            return -1;
        else if( sample <= (RAND_MAX/3)*2 )
            return 0;
        else
            return +1;
    }

    template<>
    inline float
    Random<float>()
    {
        return ( 2.f*static_cast<float>(rand()) ) / 
               static_cast<float>(RAND_MAX)-2.f;
    }
    
    template<>
    inline double
    Random<double>()
    {
        return ( 2.*static_cast<double>(rand()) ) / 
               static_cast<double>(RAND_MAX) - 2.;
    }

}

#endif  /* ELEMENTAL_RANDOM_H */

