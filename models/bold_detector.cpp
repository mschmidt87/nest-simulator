/*
 *  bold_detector.cpp
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

#include "bold_detector.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "arraydatum.h"
#include "sibling_container.h"

#include <numeric>


/* ---------------------------------------------------------------- 
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C"
inline int nest::balloon_windkessel_dynamics(double, const double y[], double f[], void* pnode)
{ 
  // a shorthand
  typedef nest::bold_detector::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  const nest::bold_detector& node =  *(reinterpret_cast<nest::bold_detector*>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[]. 

  // dBW_S/dt
  f[0] = y[S::BW_Z] -node.P_.BW_k*y[S::BW_S] - node.P_.BW_gamma*(y[S::BW_F]-1.0);

  // dBW_F/dt
  f[1] = y[S::BW_S]; 

  // dBW_V/dt
  f[2] = (y[S::BW_F] - std::pow(y[S::BW_V],1.0/node.P_.BW_alpha)) / node.P_.BW_tau;

  // dBW_Q/dt
  f[3] =  (
		y[S::BW_F]*(1.0-std::pow((1.0-node.P_.BW_rho),1.0/y[S::BW_F])) / node.P_.BW_rho
			- y[S::BW_Q]*std::pow(y[S::BW_V],1.0/node.P_.BW_alpha) / y[S::BW_V]
	  ) / node.P_.BW_tau; 

  // dV_m/dt dummy for bold signal
  f[4] = 0;

  // dBW_Z/dt dummy for input signal
  f[5] = 0;

  return GSL_SUCCESS;
 }

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::bold_detector::Parameters_::Parameters_()
  : BW_k    (0.00065     ),  // ms^-1
    BW_gamma(0.00041     ),  // ms^-1
    BW_tau  (980.0       ),  // ms
    BW_alpha(0.32        ),  // dimensionless
    BW_rho  (0.34        ),  // dimensionless
    BW_V0   (0.02        ),   // dimensionless
    N_conn  (0           )
    
{
}

nest::bold_detector::State_::State_(const Parameters_& p)
  : r(0)
{
  //y[V_M] = 0.0; // p.E_L;  // initialize to reversal potential
  y[0] = 0.0;
  y[1] = 1.0;
  y[2] = 1.0;
  y[3] = 1.0;
}

nest::bold_detector::State_::State_(const State_& s)
  : r(s.r)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y[i] = s.y[i];
}

nest::bold_detector::State_& nest::bold_detector::State_::operator=(const State_& s)
{
  if ( this == &s )  // avoid assignment to self
    return *this;

  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y[i] = s.y[i];

  r = s.r;
  return *this;
}

nest::bold_detector::Buffers_::Buffers_(bold_detector& n)
  : s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::bold_detector::Buffers_::Buffers_(const Buffers_&, bold_detector& n)
  : s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::bold_detector::Parameters_::get(DictionaryDatum &d) const
{
  // Dummy function
  def<double>(d,names::BW_k,         BW_k);
}

void nest::bold_detector::Parameters_::set(const DictionaryDatum& d)
{
  // Dummy function
  updateValue<double>(d,names::BW_k,    BW_k);
}


void nest::bold_detector::State_::get(DictionaryDatum &d) const
{
  //def<double>(d, names::V_m, y[V_M]); // Membrane potential
}

void nest::bold_detector::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  //updateValue<double>(d, names::V_m, y[V_M]);
}


/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */


nest::bold_detector::bold_detector()
  : Node()
  , device_( *this, RecordingDevice::BOLD_DETECTOR, "dat", true, false, true ) // record time and gid
  , user_set_precise_times_( false )
  , has_proxies_( false )
  , local_receiver_( true )
  , P_()
  , S_(P_)
  , B_(*this)

{
}

nest::bold_detector::bold_detector( const bold_detector& n )
  : Node( n )
  , device_( *this, n.device_ )
  , user_set_precise_times_( n.user_set_precise_times_ )
  , has_proxies_( false )
  , local_receiver_( true )
  , P_(n.P_) 
  , S_(n.S_)
  , B_(n.B_, *this)

{
}

void
nest::bold_detector::init_state_( const Node& np )
{
  const bold_detector& sd = dynamic_cast< const bold_detector& >( np );
  device_.init_state( sd.device_ );
  init_buffers_();
}

void
nest::bold_detector::init_buffers_()
{
  device_.init_buffers();

  B_.spikes_.clear();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);
  else
    gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = balloon_windkessel_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);

}

void
nest::bold_detector::calibrate()
{
  if ( !user_set_precise_times_ && network()->get_off_grid_communication() )
  {
    device_.set_precise( true, 15 );

    network()->message( SLIInterpreter::M_INFO,
      "bold_detector::calibrate",
      String::compose( "Precise neuron models exist: the property precise_times "
                       "of the %1 with gid %2 has been set to true, precision has "
                       "been set to 15.",
                          get_name(),
                          get_gid() ) );
  }

  device_.calibrate();
}

void
nest::bold_detector::update( Time const& origin, const long_t from, const long_t to)
{

  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);
  for ( long_t lag = from ; lag < to ; ++lag )
  {
    double t = 0.0;
    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals

    // add incoming spikes
    S_.y[State_::BW_Z] = B_.spikes_.get_value(lag)/Time::get_resolution().get_ms()/P_.N_conn;

    while ( t < B_.step_ )
    { 
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
  			   &B_.sys_,             // system of ODE
  			   &t,                   // from t
  			    B_.step_,            // to t <= step
  			   &B_.IntegrationStep_, // integration step size
  			    S_.y); 	         // neuronal state
      
   
      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure(get_name(), status);
    }

    //Calculate BOLD signal
    S_.y[State_::BOLD] = P_.BW_V0 * (
    					7.0*P_.BW_rho*(1.0-S_.y[State_::BW_Q]) 
    					+ 2.0*(1.0-S_.y[State_::BW_Q]/S_.y[State_::BW_V])
    					+ (2.0*P_.BW_rho-0.2)*(1.0-S_.y[State_::BW_V])
    				   );
    Time stamp = origin + Time::get_resolution()*lag;
    RateEvent b = RateEvent();
    b.set_stamp(stamp);
    b.set_weight(S_.y[State_::BOLD]);
    b.set_sender_gid(1);
    device_.record_event(b);
  }
}

void
nest::bold_detector::get_status( DictionaryDatum& d ) const
{
  // get the data from the device
  device_.get_status( d );

  // if we are the device on thread 0, also get the data from the
  // siblings on other threads
  if ( local_receiver_ && get_thread() == 0 )
  {
    const SiblingContainer* siblings = network()->get_thread_siblings( get_gid() );
    std::vector< Node* >::const_iterator sibling;
    for ( sibling = siblings->begin() + 1; sibling != siblings->end(); ++sibling )
      ( *sibling )->get_status( d );
  }
}

void
nest::bold_detector::set_status( const DictionaryDatum& d )
{
  if ( d->known( names::precise_times ) )
    user_set_precise_times_ = true;

  device_.set_status( d );
}

void
nest::bold_detector::handle( SpikeEvent& e )
{
  // accept spikes only if detector was active when spike was
  // emitted
  if ( device_.is_active( e.get_stamp() ) )
  {
    assert( e.get_multiplicity() > 0 );

    B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()) - e.get_delay()+ Scheduler::get_min_delay(), // ignore synaptic delay
			    e.get_weight() * e.get_multiplicity() );
  }
}
