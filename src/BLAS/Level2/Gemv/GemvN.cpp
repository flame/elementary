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
Elemental::BLAS::Internal::GemvN
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& x,
  const T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemvN");
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
    if( A.Height() != yLength || A.Width() != xLength )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal GemvN: " <<
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
        DistMatrix<T,MR,Star> x_MR_Star(grid);
        DistMatrix<T,MC,Star> z_MC_Star(grid);

        // Start the algorithm
        BLAS::Scal( beta, y );
        x_MR_Star.ConformWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );
        //--------------------------------------------------------------------//
        x_MR_Star = x;
        BLAS::Gemv( Normal,
                    alpha, A.LockedLocalMatrix(), 
                           x_MR_Star.LockedLocalMatrix(),
                    (T)0,  z_MC_Star.LocalMatrix()       );
        y.ReduceScatterUpdate( (T)1, z_MC_Star );
        //--------------------------------------------------------------------//
        x_MR_Star.FreeConstraints();
        z_MC_Star.FreeConstraints();
    }
    else
    {
        // Not supported
    }

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::GemvN
( const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& x,
  const float beta,        DistMatrix<float,MC,MR>& y );

template void Elemental::BLAS::Internal::GemvN
( const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& x,
  const double beta,        DistMatrix<double,MC,MR>& y );




template<typename T>
void
Elemental::BLAS::Internal::GemvN
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,VR,Star>& x,
  const T beta,        DistMatrix<T,VC,Star>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemvN");
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
    if( A.Height() != yLength || A.Width() != xLength )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal GemvN: " <<
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
        DistMatrix<T,MR,Star> x_MR_Star(grid);
        DistMatrix<T,MC,Star> z_MC_Star(grid);


        // Start the algorithm
        BLAS::Scal( beta, y );
        x_MR_Star.ConformWith( A );
        z_MC_Star.AlignWith( A );
        z_MC_Star.ResizeTo( A.Height(), 1 );


        //--------------------------------------------------------------------//
        x_MR_Star = x;
        BLAS::Gemv( Normal,
                    alpha, A.LockedLocalMatrix(), 
                           x_MR_Star.LockedLocalMatrix(),
                    (T)0,  z_MC_Star.LocalMatrix()       );


        // Hacking:
        // We cannot invoke y.ReduceScatterUpdate( (T)1, z_MC_Star );
        // Since there is no y_VC_Star.ReduceScatterUpdate function.
        DistMatrix<T,MC,MR> y_MC_MR(grid);
        y_MC_MR.AlignWith( y );
        y_MC_MR = y;
        y_MC_MR.ReduceScatterUpdate( (T)1, z_MC_Star );
        y = y_MC_MR;

        //--------------------------------------------------------------------//
        x_MR_Star.FreeConstraints();
        z_MC_Star.FreeConstraints();
        y_MC_MR.FreeConstraints();
    }
    else
    {
        // Not supported
    }

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::GemvN
( const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,VR,Star>& x,
  const float beta,        DistMatrix<float,VC,Star>& y );

template void Elemental::BLAS::Internal::GemvN
( const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,VR,Star>& x,
  const double beta,        DistMatrix<double,VC,Star>& y );


