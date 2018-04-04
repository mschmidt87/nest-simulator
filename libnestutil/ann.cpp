/*
 *  ann.cpp
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

#include <cmath>
#include <cstdlib>

#include "ann.h"
#include "ann_util.h"

namespace nest
{

ANN::ANN()
  : n_inputs_( 1 )
  , n_layers_( 1 )
{
  n_units_per_layer_ = std::vector< size_t >( 1 );
  n_units_per_layer_[ 0 ] = 1;
  resize_();
}

ANN::ANN( const size_t n_inputs, const size_t n_layers, const std::vector< size_t >& n_units_per_layer )
  : n_inputs_( n_inputs )
  , n_layers_( n_layers )
{
  n_units_per_layer_ = n_units_per_layer;
  resize_();
}

size_t
ANN::get_n_pre_( const size_t layer ) const
{
  return layer == 0 ? n_inputs_ : n_units_per_layer_[ layer - 1 ];
}

void
ANN::resize_()
{
  weight_.resize( n_layers_ );
  bias_.resize( n_layers_ );
  // activities_.resize( n_layers_ );

  for ( size_t layer = 0; layer < n_layers_; ++layer )
  {
    weight_[ layer ].resize( n_units_per_layer_[ layer ] );
    for ( size_t unit = 0; unit < weight_[ layer ].size(); ++unit )
    {
      weight_[ layer ][ unit ].resize( get_n_pre_( layer ) );
    }

    bias_[ layer ].resize( n_units_per_layer_[ layer ] );
    // activities_[ layer ].resize( n_units_per_layer_[ layer ] );
  }
}

void
ANN::forward( const std::vector< double >& input, std::vector< double >& output ) const
{
  // reset_activities_();
  std::vector< std::vector< double > > activities;
  activities.resize( n_layers_ );
  for ( size_t layer = 0; layer < n_layers_; ++layer )
  {
    activities[ layer ].resize( n_units_per_layer_[ layer ] );
    const std::vector< double >& pre_activities = layer == 0 ? input : activities[ layer - 1 ];
    mv( weight_[ layer ], pre_activities, activities[ layer ] );
    vadd( activities[ layer ], bias_[ layer ], activities[ layer ] );
    if (layer != n_layers_ - 1)
    {
      for ( size_t unit = 0; unit < activities[ layer ].size(); ++unit )
      {
	activities[ layer ][ unit ] = relu( activities[ layer ][ unit ] );
      }
    }
  }
  output = activities[ n_layers_ - 1];
}


void
ANN::forward( const std::vector< std::vector< double > >& input, std::vector< std::vector< double > >& output ) const
{
  output.resize( input.size() );
  for ( size_t input_index = 0; input_index < input.size(); ++input_index )
  {
    forward( input[ input_index ], output[ input_index ] );
  }
}

void
ANN::randomize_parameters( const double scaling )
{
  for ( size_t layer = 0; layer < n_layers_; ++layer )
  {
    const double std = 1. / std::sqrt( get_n_pre_( layer ) ) * scaling;
    for ( size_t unit = 0; unit < weight_[ layer ].size(); ++unit )
    {
      for ( size_t pre = 0; pre < weight_[ layer ][ unit ].size(); ++pre )
      {
        weight_[ layer ][ unit ][ pre ] = 2 * myrand() - 1. * std;
      }
      bias_[ layer ][ unit ] = 2 * myrand() - 1. * std;
    }
  }
}

// void
// ANN::reset_activities_()
// {
//   for ( size_t layer = 0; layer < n_layers_; ++layer )
//   {
//     for ( size_t unit = 0; unit < weight_[ layer ].size(); ++unit )
//     {
//       activities_[ layer ][ unit ] = 0.;
//     }
//   }
// }

void
ANN::seed( const size_t seed )
{
  srand( seed );
}

void
ANN::set_weight( const std::vector< std::vector< std::vector< double > > >& weight )
{
  weight_ = weight;
}

void
ANN::set_bias( const std::vector< std::vector< double > >& bias )
{
  bias_ = bias;
}

} // namespace ann
