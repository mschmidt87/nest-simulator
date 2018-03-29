/*
 *  ann.h
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

#ifndef ANN_H
#define ANN_H

#include <numeric>
#include <vector>
#include <cstddef>

namespace nest
{

class ANN
{
// private:
public:
  size_t n_inputs_;
  size_t n_layers_;
  std::vector< size_t > n_units_per_layer_;
  std::vector< std::vector< std::vector< double > > > weight_;
  std::vector< std::vector< double > > bias_;
  // std::vector< std::vector< double > > activities_;

  size_t get_n_pre_( const size_t layer ) const;
  // void reset_activities_();
  void resize_();

public:
  ANN();
  ANN( const size_t n_inputs, const size_t n_layers, const std::vector< size_t >& n_units_per_layer );

  void forward( const std::vector< double >& input, std::vector< double >& output ) const;
  void forward( const std::vector< std::vector< double > >& input, std::vector< std::vector< double > >& output ) const;
  void randomize_parameters( const double scaling = 1. );
  void seed( const size_t seed );
  void set_weight( const std::vector< std::vector< std::vector< double > > >& weight );
  void set_bias( const std::vector< std::vector< double > >& bias );

};

} // namespace nest

#endif
