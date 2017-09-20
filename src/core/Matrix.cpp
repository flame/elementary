/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "ElementalEnvironment.h"
#include "ElementalMatrix.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::Matrix<T>::Print( const string msg ) const
{
#ifndef RELEASE
    PushCallStack("Matrix::Print");
#endif
    if( msg != "" )
        cout << msg << endl;

    const int height = Height();
    const int width = Width();

    for( int i=0; i<height; ++i )
    {
        for( int j=0; j<width; ++j )
        {
            cout << operator()(i,j) << " ";
        }
        cout << endl;
    }
    cout << endl;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::ResizeTo
( const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width)");
    if( height < 0 || width < 0 )
    {
        cerr << "Height and width must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
    if( _viewing )
    {
        cerr << "Does not make sense to resize matrix when viewing other data."
             << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    const int minLDim = 1;
    _height = height;
    _width  = width;
    _ldim   = max( height, minLDim );

    _memory.Require(_ldim*width);
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::ResizeTo
( const int height, const int width, const int ldim )
{
#ifndef RELEASE
    PushCallStack("Matrix::ResizeTo(height,width,ldim)");
    if( height < 0 || width < 0 )
    {
        cerr << "Height and width must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
    if( _viewing )
    {
        cerr << "Does not make sense to resize matrix when viewing other data."
             << endl;
        DumpCallStack();
        throw exception();
    }
    if( ldim < height )
    {
        cerr << "Tried to set ldim(" << ldim << ") < height (" << height << ")"
             << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = height;
    _width  = width;
    _ldim   = ldim;

    _memory.Require(ldim*width);
    _data = _memory.Buffer();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::View( Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A)");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with Matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = A.Height();
    _width  = A.Width();
    _ldim   = A.LDim();
    _data   = A.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::LockedView( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A)");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with Matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height     = A.Height();
    _width      = A.Width();
    _ldim       = A.LDim();
    _lockedData = A.LockedBuffer();
    _viewing    = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::View
( Matrix<T>& A, 
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::View(A,i,j,height,width)");
    if( i < 0 || j < 0 )
    {
        cerr << "Indices must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
    if( height < 0 || width < 0 )
    {
        cerr << "Height and width must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with Matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        cerr << "Trying to view outside of a Matrix: " 
             << "up to (" << i+height-1 << "," << j+width-1 << ") "
             << "of " << A.Height() << " x " << A.Width() << " Matrix." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height     = height;
    _width      = width;
    _ldim       = A.LDim();
    _data       = A.Buffer(i,j,height,width);
    _viewing    = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::LockedView
( const Matrix<T>& A, 
  const int i, const int j, const int height, const int width )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView(A,i,j,height,width)");
    if( i < 0 || j < 0 )
    {
        cerr << "Indices must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
    if( height < 0 || width < 0 )
    {
        cerr << "Height and width must be non-negative." << endl;
        DumpCallStack();
        throw exception();
    }
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with Matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( (i+height) > A.Height() || (j+width) > A.Width() )
    {
        cerr << "Trying to view outside of a Matrix: " 
             << "up to (" << i+height-1 << "," << j+width-1 << ") "
             << "of " << A.Height() << " x " << A.Width() << " Matrix." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height     = height;
    _width      = width;
    _ldim       = A.LDim();
    _lockedData = A.LockedBuffer(i,j,height,width);
    _viewing    = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::View1x2
( Matrix<T>& AL, Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View1x2");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with Matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AL.Height() != AR.Height() )
    {
        cerr << "1x2 must have consistent height to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AL.LDim() != AR.LDim() )
    {
        cerr << "1x2 must have consistent ldims to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AR.Buffer() != (AL.Buffer()+AL.LDim()*AL.Width()) )
    {
        cerr << "1x2 must have contiguous memory." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = AL.Height();
    _width  = AL.Width() + AR.Width();
    _ldim   = AL.LDim();
    _data   = AL.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::LockedView1x2
( const Matrix<T>& AL, const Matrix<T>& AR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView1x2");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with Matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AL.Height() != AR.Height() )
    {
        cerr << "1x2 must have consistent height to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AL.LDim() != AR.LDim() )
    {
        cerr << "1x2 must have consistent ldims to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AR.LockedBuffer() != (AL.LockedBuffer()+AL.LDim()*AL.Width()) )
    {
        cerr << "1x2 must have contiguous memory." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height     = AL.Height();
    _width      = AL.Width() + AR.Width();
    _ldim       = AL.LDim();
    _lockedData = AL.LockedBuffer();
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::View2x1( Matrix<T>& AT, 
                           Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x1");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AT.Width() != AB.Width() )
    {
        cerr << "2x1 must have consistent width to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AT.LDim() != AB.LDim() )
    {
        cerr << "2x1 must have consistent ldim to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AB.Buffer() != (AT.Buffer() + AT.Height()) )
    {
        cerr << "2x1 must have contiguous memory." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = AT.Height() + AB.Height();
    _width  = AT.Width();
    _ldim   = AT.LDim();
    _data   = AT.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::LockedView2x1( const Matrix<T>& AT, 
                                 const Matrix<T>& AB )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x1");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing with matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AT.Width() != AB.Width() )
    {
        cerr << "2x1 must have consistent width to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AT.LDim() != AB.LDim() )
    {
        cerr << "2x1 must have consistent ldim to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( AB.LockedBuffer() != (AT.LockedBuffer()+AT.Height()) )
    {
        cerr << "2x1 must have contiguous memory." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height     = AT.Height() + AB.Height();
    _width      = AT.Width();
    _ldim       = AT.LDim();
    _lockedData = AT.LockedBuffer();
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::View2x2( Matrix<T>& ATL, Matrix<T>& ATR,
                           Matrix<T>& ABL, Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::View2x2");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing a matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
    {
        cerr << "2x2 must conform to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
    {
        cerr << "2x2 must have consistent ldims to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ABL.Buffer() != (ATL.Buffer() + ATL.Height()) ||
        ABR.Buffer() != (ATR.Buffer() + ATR.Height()) ||
        ATR.Buffer() != (ATL.Buffer() + ATL.LDim()*ATL.Width()) )
    {
        cerr << "2x2 must have contiguous memory." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _ldim   = ATL.LDim();
    _data   = ATL.Buffer();
    _viewing = true;
    _lockedView = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::LockedView2x2
( const Matrix<T>& ATL, const Matrix<T>& ATR,
  const Matrix<T>& ABL, const Matrix<T>& ABR )
{
#ifndef RELEASE
    PushCallStack("Matrix::LockedView2x2");
    if( _memory.Size() > 0 )
    {
        cerr << "Viewing a matrix after allocating memory." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ATL.Width() != ABL.Width()   ||
        ATR.Width() != ABR.Width()   ||
        ATL.Height() != ATR.Height() ||
        ABL.Height() != ABR.Height()   )
    {
        cerr << "2x2 must conform to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ATL.LDim() != ATR.LDim() ||
        ATR.LDim() != ABL.LDim() ||
        ABL.LDim() != ABR.LDim()   )
    {
        cerr << "2x2 must have consistent ldims to combine." << endl;
        DumpCallStack();
        throw exception();
    }
    if( ABL.LockedBuffer() != (ATL.LockedBuffer() + ATL.Height()) ||
        ABR.LockedBuffer() != (ATR.LockedBuffer() + ATR.Height()) ||
        ATR.LockedBuffer() != (ATL.LockedBuffer() + ATL.LDim()*ATL.Width()) )
    {
        cerr << "2x2 must have contiguous memory." << endl;
        DumpCallStack();
        throw exception();
    }

#endif
    _height = ATL.Height() + ABL.Height();
    _width  = ATL.Width() + ATR.Width();
    _ldim   = ATL.LDim();
    _lockedData = ATL.LockedBuffer();
    _viewing = true;
    _lockedView = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::SetToIdentity()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToIdentity");
    if( _lockedView )
    {
        cerr << "Cannot set a locked view to identity." << endl; 
        DumpCallStack();
        throw exception();
    }
#endif
    const int height = Height();
    const int width = Width();

    SetToZero();
    for( int j=0; j<min(height,width); ++j )
        _data[j+j*_ldim] = (T)1;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::SetToZero()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToZero");
    if( _lockedView )
    {
        cerr << "Cannot set a locked view to zero." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    const int height = Height();
    const int width = Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*_ldim] = (T)0;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::Matrix<T>::SetToRandom()
{
#ifndef RELEASE
    PushCallStack("Matrix::SetToRandom");
    if( _lockedView )
    {
        cerr << "Cannot change the data of a locked view." << endl; 
        DumpCallStack();
        throw exception();
    }
#endif
    const int height = Height();
    const int width = Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*_ldim] = Random<T>();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
const Matrix<T>&
Elemental::Matrix<T>::operator=
( const Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Matrix::operator=");
    if( _lockedView )
    {
        cerr << "Cannot asign to a locked view." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    ResizeTo( A.Height(), A.Width() );

    const int height = Height();
    const int width = Width();
    const int ldim = LDim();
    const int ldimOfA = A.LDim();
    const T* data = A.LockedBuffer();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            _data[i+j*ldim] = data[i+j*ldimOfA];
#ifndef RELEASE
    PopCallStack();
#endif
    return *this;
}

template class Elemental::Matrix<int>;
template class Elemental::Matrix<float>;
template class Elemental::Matrix<double>;

