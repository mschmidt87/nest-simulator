/*
 *  ann_synapse.h
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


/* BeginDocumentation // TODO BROKEN
   Name: ann_synapse - Synapse type for static connections with
   homogeneous weight.

   Description:
     ann_synapse does not support any kind of plasticity. It simply
     stores the parameters delay, target, and receiver port for each connection
     and uses a common weight for all connections.

   Remarks:
     The common weight for all connections of this model must be set by
     SetDefaults on the model. If you create copies of this model using
     CopyModel, each derived model can have a different weight.

   Transmits: SpikeEvent, RateEvent, CurrentEvent, ConductanceEvent,
   DataLoggingRequest, DoubleDataEvent

   Parameters:
     No Parameters

   References:
     No References
   FirstVersion: April 2008
   Author: Susanne Kunkel, Moritz Helias
   SeeAlso: synapsedict, static_synapse
*/

#ifndef ANN_SYNAPSE_H
#define ANN_SYNAPSE_H

// C++ includes:
#include <string>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"

// Includes from libnestutil:
#include "ann.h"

namespace nest
{

/**
 * Class containing the common properties for all synapses with common weight.
 */
class CommonPropertiesANN : public CommonSynapseProperties
{
public:
  // data members common to all connections
  ANN ann_;

public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  CommonPropertiesANN()
    : CommonSynapseProperties()
  {
  }

  /**
   * Get all properties and put them into a dictionary.
   */
  void
  get_status( DictionaryDatum& d ) const
  {
    CommonSynapseProperties::get_status( d );

    DictionaryDatum weight = new Dictionary();
    DictionaryDatum bias = new Dictionary();
    for ( size_t layer = 0; layer < ann_.weight_.size(); ++layer )
    {
      DictionaryDatum weight_inner = new Dictionary();
      for ( size_t unit = 0; unit < ann_.weight_[ layer ].size(); ++unit )
      {
        initialize_property_doublevector( weight_inner, std::to_string( unit ) );
        append_property( weight_inner, std::to_string( unit ), ann_.weight_[ layer ][ unit ] );
      }
      ( *weight )[ Name( std::to_string( layer ) ) ] = weight_inner;
      initialize_property_doublevector( bias, std::to_string( layer ) );
      append_property( bias, std::to_string( layer ), ann_.bias_[ layer ] );
    }
    ( *d )[ Name( "ann_weight" ) ] = weight;
    ( *d )[ Name( "ann_bias" ) ] = bias;
  }

  /**
   * Set properties from the values given in dictionary.
   */
  void
  set_status( const DictionaryDatum& d, ConnectorModel& cm )
  {
    CommonSynapseProperties::set_status( d, cm );

    if ( d->known( "ann_weight" ) )
    {
      const DictionaryDatum weight = getValue< DictionaryDatum >( d, "ann_weight" );
      for ( size_t layer = 0; layer < ann_.weight_.size(); ++layer )
      {
        const DictionaryDatum weight_inner = getValue< DictionaryDatum >( weight, std::to_string( layer ) );
        for ( size_t unit = 0; unit < ann_.weight_[ layer ].size(); ++unit )
        {
          ann_.weight_[ layer ][ unit ] = getValue< std::vector< double > >( weight_inner, std::to_string( unit ) );
        }
      }
    }

    if ( d->known( "ann_bias" ) )
    {
      const DictionaryDatum bias = getValue< DictionaryDatum >( d, "ann_bias" );
      for ( size_t layer = 0; layer < ann_.weight_.size(); ++layer )
      {
        ann_.bias_[ layer ] = getValue< std::vector< double > >( bias, std::to_string( layer ) );
      }
    }
  }
};

/**
 * Class representing a static connection. A static connection has the
 * properties weight, delay and receiver port. A suitable Connector containing
 * these connections can be obtained from the template GenericConnector.
 */
template < typename targetidentifierT >
class ANNSynapse : public Connection< targetidentifierT >
{
  double weight_;

public:
  // this line determines which common properties to use
  typedef CommonPropertiesANN CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;
  using ConnectionBase::get_delay_steps;

  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( RateEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DataLoggingRequest&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( CurrentEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( ConductanceEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DoubleDataEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DSSpikeEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DSCurrentEvent&, rport )
    {
      return invalid_port_;
    }
  };


  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    double,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
  }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   * \param t_lastspike Time point of last spike emitted
   */
  void
  send( Event& e, thread t, double, const CommonPropertiesANN& cp )
  {
    std::vector< double > input = std::vector< double >( 1 );
    input[ 0 ] = weight_;
    std::vector< double > output = std::vector< double >( 1 );
    cp.ann_.forward( input, output );
    weight_ = output[ 0 ];

    e.set_weight( weight_ );
    e.set_delay( get_delay_steps() );
    e.set_receiver( *get_target( t ) );
    e.set_rport( get_rport() );
    e();
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }
};


template < typename targetidentifierT >
void
ANNSynapse< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
ANNSynapse< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
}


} // namespace

#endif /* #ifndef ANN_SYNAPSE_H */
