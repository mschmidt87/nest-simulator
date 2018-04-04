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
  // ANN ann_;
  ANN ann_;

  
public:
  /**
   * Default constructor.
   * Sets all property values to defaults.
   */
  CommonPropertiesANN()
    : CommonSynapseProperties()
  {
    // data members common to all connections
    const size_t n_inputs = 2;
    const size_t n_layers = 3;
  
    std::vector< size_t > n_units_per_layer ( n_layers );
    n_units_per_layer[ 0 ] = 10;
    n_units_per_layer[ 1 ] = 10;
    n_units_per_layer[ 2 ] = 1;
    ann_ = ANN( n_inputs, n_layers, n_units_per_layer );
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
    ( *d )[ Name( "ann_n_inputs") ] = ann_.n_inputs_;
    ( *d )[ Name( "ann_n_layers") ] = ann_.n_layers_;
    ( *d )[ Name( "ann_n_units_per_layer") ] = ann_.n_units_per_layer_;
  }

  /**
   * Set properties from the values given in dictionary.
   */
  void
  set_status( const DictionaryDatum& d, ConnectorModel& cm )
  {
    CommonSynapseProperties::set_status( d, cm );

    if ( d->known( "ann_n_inputs" ) )
    {
      ann_.n_inputs_ = getValue< long >( d, "ann_n_inputs" );
    }
    if ( d->known( "ann_n_layers" ) )
    {
      ann_.n_layers_ = getValue< long >( d, "ann_n_layers" );
    }
    if ( d->known( "ann_n_units_per_layer" ) )
    {
      ann_.n_units_per_layer_ = getValue< std::vector< size_t > >( d, "ann_n_units_per_layer" );
    }
    ann_.resize_();

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
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

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
    double t_lastspike,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
  }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   * \param t_lastspike Time point of last spike emitted
   */
  void
  send( Event& e, thread t, double t_lastspike, const CommonPropertiesANN& cp )
  {
    std::vector< double > input = std::vector< double >( cp.ann_.n_inputs_ );
    
    input[ 0 ] = weight_;

    Node* target = get_target( t );
    double dendritic_delay = get_delay();

    double t_spike = e.get_stamp().get_ms();
    // input[ 1 ] = t_spike - dendritic_delay;

    // get spike history in relevant range (t1, t2] from post-synaptic neuron
    std::deque< histentry >::iterator start;
    std::deque< histentry >::iterator finish;

    target->get_history(
      t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish );
    // facilitation due to post-synaptic spikes since last pre-synaptic spike
    double minus_dt;
    while ( start != finish )
    {
      minus_dt = t_lastspike - ( start->t_ + dendritic_delay );
      ++start;
      if ( minus_dt == 0 )
      {
	continue;
      }
      input[ 1 ] = minus_dt;
      std::vector< double > output = std::vector< double >( cp.ann_.n_units_per_layer_.back() );
      cp.ann_.forward( input, output );
      weight_ = output[ 0 ];
    }

    // depression due to new pre-synaptic spike
    input[ 1 ] = target->get_K_value( t_spike - dendritic_delay );
    std::vector< double > output = std::vector< double >( cp.ann_.n_units_per_layer_.back() );
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
