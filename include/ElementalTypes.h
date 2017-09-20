/*
   Copyright (c) The University of Texas at Austin, 2009-2017.
   Copyright (c) Jack Poulson, 2009-2017.

   This file is part of Elementary and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_TYPES_H
#define ELEMENTAL_TYPES_H 1

#ifndef WITHOUT_COMPLEX
#include <complex>
#endif
#include <iostream>
#include <limits>
#include <string>

namespace Elemental
{
#ifndef WITHOUT_COMPLEX
    typedef std::complex<float>  scomplex; 
    typedef std::complex<double> dcomplex;
#endif

    enum Diagonal
    {
        NonUnit,
        Unit
    };

    inline char DiagonalToChar( const Diagonal diagonal )
    {
        char diagonalChar;
        switch( diagonal )
        {
            case NonUnit: diagonalChar = 'N'; break;
            default:      diagonalChar = 'U'; break;
        }
        return diagonalChar;
    }

    inline Diagonal CharToDiagonal( const char c )
    {
        Diagonal diagonal;
        switch( c )
        {
            case 'N': diagonal = NonUnit; break;
            case 'U': diagonal = Unit;    break;
            default:
                std::cerr << "CharToDiagonal expects char in {N,U}." 
                          << std::endl;
                throw std::exception();     
        }
        return diagonal;
    }

    enum Distribution
    {
        MC,  // Col of a matrix distribution
        MD,  // Diagonal of a matrix distribution
        MR,  // Row of a matrix distribution
        VC,  // Col-major vector distribution
        VR,  // Row-major vector distribution
        Star // Do not distribute
    };

    inline std::string DistToString( const Distribution distribution )
    {
        std::string distString;
        switch( distribution )
        {
            case MC: distString = "MC"; break;
            case MD: distString = "MD"; break;
            case MR: distString = "MR"; break;
            case VC: distString = "VC"; break;
            case VR: distString = "VR"; break;
            default: distString = "* "; break;
        }
        return distString;
    }

    inline Distribution StringToDist( const std::string s )
    {
        Distribution distribution;
        if( s == "MC" )
            distribution = MC;
        else if( s == "MD" )
            distribution = MD;
        else if( s == "MR" )
            distribution = MR;
        else if( s == "VC" )
            distribution = VC;
        else if( s == "VR" )
            distribution = VR;
        else if( s == "* " || s == " *" || s == "*" )
            distribution = Star;
        else
        {
            std::cerr 
                << "StringToDist expects string in "
                << "{\"MC\",\"MD\",\"MR\",\"VC\",\"VR\",\"* \",\" *\",\"*\"}."
                << std::endl;
            throw std::exception();     
        }
        return distribution;
    }

    enum Orientation
    {
        Normal,
        Transpose,
        ConjugateTranspose
    };

    inline char OrientationToChar( const Orientation orientation )
    {
        char orientationChar;
        switch( orientation )
        {
            case Normal:             orientationChar = 'N'; break;
            case Transpose:          orientationChar = 'T'; break;
            default:                 orientationChar = 'C'; break;
        }
        return orientationChar;
    }

    inline Orientation CharToOrientation( const char c )
    {
        Orientation orientation;
        switch( c )
        {
            case 'N': orientation = Normal;             break;
            case 'T': orientation = Transpose;          break;
            case 'C': orientation = ConjugateTranspose; break;
            default:
                std::cerr << "CharToOrientation expects char in {N,T,C}." 
                          << std::endl;
                throw std::exception();     
        }
        return orientation;
    }

    enum Shape
    {
        Lower,
        Upper
    };

    inline char ShapeToChar( const Shape shape )
    {
        char shapeChar;
        switch( shape )
        {
            case Lower: shapeChar = 'L'; break;
            default:    shapeChar = 'U'; break;
        }
        return shapeChar;
    }

    inline Shape CharToShape( const char c )
    {
        Shape shape;
        switch( c )
        {
            case 'L': shape = Lower; break;
            case 'U': shape = Upper; break;
            default:
                std::cerr << "CharToShape expects char in {L,U}." 
                          << std::endl;
                throw std::exception();     
        }
        return shape;
    }

    enum Side
    {
        Left,
        Right
    };

    inline char SideToChar( const Side side )
    {
        char sideChar;
        switch( side )
        {
            case Left:  sideChar = 'L'; break;
            default:    sideChar = 'R'; break;
        }
        return sideChar;
    }
    
    inline Side CharToSide( const char c )
    {
        Side side;
        switch( c )
        {
            case 'L': side = Left;  break;
            case 'R': side = Right; break;
            default:
                std::cerr << "CharToSide expects char in {L,R}." 
                          << std::endl;
                throw std::exception();     
        }
        return side;
    }

}

#endif /* ELEMENTAL_TYPES_H */

