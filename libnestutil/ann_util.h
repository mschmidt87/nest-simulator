/*
 *  ann_util.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ANN_UTIL_H
#define ANN_UTIL_H

#include <cstdlib>
#include <iostream>
#include <vector>

namespace nest
{

inline double
myrand()
{
  return static_cast< double >( rand() ) / RAND_MAX;
}

template< typename T >
inline void
vv( const std::vector< T >& v0, const std::vector< T >& v1, T& o )
{
  for ( size_t i = 0; i < v0.size(); ++i )
  {
    o += v0[ i ] * v1[ i ];
  }
}

template< typename T >
inline void
mv( const std::vector< std::vector< T > >& m, const std::vector< T > v, std::vector< T >& o )
{
  for ( size_t i = 0; i < m.size(); ++i )
  {
    vv( m[ i ], v, o[ i ] );
  }
}

template< typename T >
inline void
vadd( const std::vector< T >& v0, const std::vector< T >& v1, std::vector< T >& o )
{
  for ( size_t i = 0; i < v0.size(); ++i )
  {
    o[ i ] = v0[ i ] + v1[ i ];
  }
}

inline double
relu( const double x )
{
  return x > 0 ? x : 0;
}

} // namespace nest

#endif
