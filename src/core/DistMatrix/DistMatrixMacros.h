/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#pragma once

#ifndef RELEASE

#define REPORT_UNIMPLEMENTED_FEATURE \
{ \
    std::cerr << "Sorry, feature not yet implemented." << std::endl; \
    DumpCallStack(); \
    throw std::exception(); \
}
#else
#define REPORT_UNIMPLEMENTED_FEATURE \
{ \
    std::cerr << "Sorry, feature not yet implemented." << std::endl; \
    std::cerr << "Run debug version to see which routine." << std::endl; \
    throw std::exception(); \
}
#endif

#ifndef RELEASE

#define CHECK_IF_ALIGNING_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
    { \
        std::cerr << "Tried to align with a matrix from diff. grid." \
                  << std::endl; \
        DumpCallStack(); \
        throw exception(); \
    } \
}

#define CHECK_IF_CONFORMING_1x2( AL, AR ) \
{ \
    if( AL.Height() != AR.Height() ) \
    { \
        std::cerr << "1x2 not conformant. Left is " << AL.Height() \
                  << " x " << AL.Width() << ", right is " \
                  << AR.Height() << " x " << AR.Width() << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}
#define CHECK_IF_CONFORMING_2x1( AT, AB ) \
{ \
    if( AT.Width() != AB.Width() ) \
    { \
        std::cerr << "2x1 is not conformant. Top is " << AT.Height() \
                  << " x " << AT.Width() << ", bottom is " \
                  << AB.Height() << " x " << AB.Width() << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}
#define CHECK_IF_CONFORMING_2x2( ATL, ATR, ABL, ABR ) \
{ \
    if( ATL.Width() != ABL.Width()   ||   \
        ATR.Width() != ABR.Width()   ||   \
        ATL.Height() != ATR.Height() ||   \
        ABL.Height() != ABR.Height()    ) \
    { \
        std::cerr << "2x2 is not conformant." << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_CONFORMING_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
    { \
        std::cerr << "Tried to conform with a matrix from diff. grid." \
                  << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_LOCKED_VIEW \
{ \
    if( _viewing && _lockedView ) \
    { \
        std::cerr << "Cannot alter data with locked view." << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_UNFREED_COL_CONSTRAINT \
{ \
    if( ConstrainedColDist() ) \
    { \
        std::cerr << "Forgot to free column constraint before reconstraining." \
                  << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_UNFREED_ROW_CONSTRAINT \
{ \
    if( ConstrainedRowDist() ) \
    { \
        std::cerr << "Forgot to free row constraint before reconstraining." \
                  << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_NONEXISTANT_CONSTRAINTS \
{ \
    if( !ConstrainedColDist() && !ConstrainedRowDist() ) \
    { \
        std::cerr << "Tried to free nonexistant constraints." << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_OUT_OF_BOUNDS( A, i, j, height, width ) \
{ \
    if( i < 0 || j < 0 ) \
    { \
        std::cerr << "Indices must be non-negative." << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
    if( height < 0 || width < 0 ) \
    { \
        std::cerr << "Height and width must be non-negative." << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
    if( (i+height) > A.Height() || (j+width) > A.Width() ) \
    { \
        std::cerr << "Out of bounds of distributed matrix: up to (" \
                  << i+height-1 << "," << j+width-1 << ") of " \
                  << A.Height() << "x" << A.Width() << " matrix." << std::endl;\
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_REDIST_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
    { \
        std::cerr << \
        "Tried to redist. btw. matrices from different grids." \
        << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_DIFF_SIZE( A ) \
{ \
    if( Height() != A.Height() || Width() != A.Width() ) \
    { \
        std::cerr << "Distributed matrices are different sizes." << std::endl; \
        DumpCallStack(); \
        throw exception(); \
    } \
}

#define CHECK_IF_VIEWING_DIFF_SIZE( A ) \
{ \
    if( _viewing && (Height() != A.Height() || Width() != A.Width()) ) \
    { \
        std::cerr << "Cannot assign to different sized view." \
                  << "Views cannot be resized." << std::endl; \
        DumpCallStack(); \
        throw exception(); \
    } \
}

#define CHECK_IF_VIEWING_AND_STORING \
{ \
    if( LockedLocalMatrix().MemorySize() > 0 ) \
    { \
        std::cerr << "Attempting to view and store with distributed matrix." \
                  << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#define CHECK_IF_VIEWING_DIFF_GRID( A ) \
{ \
    if( GetGrid() != A.GetGrid() ) \
    { \
        std::cerr << "Tried to view a matrix from a different grid." \
                  << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}

#endif

