/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#ifndef ELEMENTAL_MEMORY_H
#define ELEMENTAL_MEMORY_H 1

namespace Elemental
{
    template<typename T>
    class Memory
    {
        size_t _size;
        T*     _buffer;
    public:
        Memory();
        Memory( size_t size );
        ~Memory();

        T*     Buffer() const;
        size_t Size()   const;

        void   Require( size_t size );
        void   Release();
    };
}

/*----------------------------------------------------------------------------*/

template<typename T>
inline
Elemental::Memory<T>::Memory()
: _size(0), _buffer(NULL)
{ }

template<typename T>
inline
Elemental::Memory<T>::Memory( size_t size )
: _size(size), _buffer(new T[size])
{ }

template<typename T>
inline
Elemental::Memory<T>::~Memory()
{
    delete[] _buffer;
}

template<typename T>
inline T*
Elemental::Memory<T>::Buffer() const
{
    return _buffer;
}

template<typename T>
inline size_t
Elemental::Memory<T>::Size() const
{
    return _size;
}

template<typename T>
inline void
Elemental::Memory<T>::Require
( size_t size )
{
    if( size > _size )
    {
        delete[] _buffer;
        _buffer = new T[size];
        _size = size;
    }
}

template<typename T>
inline void
Elemental::Memory<T>::Release()
{
#ifndef POOL_MEMORY
    delete[] _buffer;
    _size = 0;
    _buffer = 0;
#endif
}

#endif /* ELEMENTAL_MEMORY_H */

