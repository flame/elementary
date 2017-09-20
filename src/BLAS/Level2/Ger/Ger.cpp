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
Elemental::BLAS::Ger
( const T alpha, const DistMatrix<T,MC,MR>& x,
                 const DistMatrix<T,MC,MR>& y,
                       DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Ger");
#endif
    const Grid& grid = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "{A,x,y} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y are assumed to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal Ger: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            endl << "  x ~ " << x.Height() << " x " << x.Width() << 
            endl << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    if( x.Width() == 1 && y.Width() == 1 )
    {
        //cout << " x.Width() == 1 && y.Width() == 1 " << endl;
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(grid);
        DistMatrix<T,MR,Star> y_MR_Star(grid);

        // Begin the algoritm
        x_MC_Star.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        y_MR_Star = y;
        BLAS::Ger( alpha, x_MC_Star.LockedLocalMatrix(),
                          y_MR_Star.LockedLocalMatrix(),
                          A.LocalMatrix()               );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        y_MR_Star.FreeConstraints();
    }
    else
    {
        // Not supported
    }
    
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Ger
( const float alpha, const DistMatrix<float,MC,MR>& x,
                     const DistMatrix<float,MC,MR>& y,
                           DistMatrix<float,MC,MR>& A );

template void Elemental::BLAS::Ger
( const double alpha, const DistMatrix<double,MC,MR>& x,
                      const DistMatrix<double,MC,MR>& y,
                            DistMatrix<double,MC,MR>& A );


template<typename T>
void
Elemental::BLAS::Ger
( const T alpha, const DistMatrix<T,VR,Star>& x,
                 const DistMatrix<T,VC,Star>& y,
                       DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Ger");
#endif
    const Grid& grid = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "{A,x,y} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
    {
        if( grid.VCRank() == 0 )
            cerr << "x and y are assumed to be vectors." << endl;
        DumpCallStack();
        throw exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal Ger: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            endl << "  x ~ " << x.Height() << " x " << x.Width() << 
            endl << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    if( x.Width() == 1 && y.Width() == 1 )
    {
        //cout << " x.Width() == 1 && y.Width() == 1 " << endl;
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(grid);
        DistMatrix<T,MR,Star> y_MR_Star(grid);

        // Begin the algoritm
        x_MC_Star.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        y_MR_Star = y;
        BLAS::Ger( alpha, x_MC_Star.LockedLocalMatrix(),
                          y_MR_Star.LockedLocalMatrix(),
                          A.LocalMatrix()               );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        y_MR_Star.FreeConstraints();
    }
    else
    {
        // Not supported
    }
    
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Ger
( const float alpha, const DistMatrix<float,VR,Star>& x,
                     const DistMatrix<float,VC,Star>& y,
                           DistMatrix<float,MC,MR>& A );

template void Elemental::BLAS::Ger
( const double alpha, const DistMatrix<double,VR,Star>& x,
                      const DistMatrix<double,VC,Star>& y,
                            DistMatrix<double,MC,MR>& A );

