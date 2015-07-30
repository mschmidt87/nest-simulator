/*
 *  iaf_cond_alpha.cpp
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


#include "balloon_windkessel.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"
#include <limits>

#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<nest::balloon_windkessel> nest::balloon_windkessel::recordablesMap_;

namespace nest   // template specialization must be placed in namespace
{
  /*
   * Override the create() method with one call to RecordablesMap::insert_() 
   * for each quantity to be recorded.
   */
  template <>
  void RecordablesMap<balloon_windkessel>::create()
  {
    // use standard names whereever you can for consistency!

    insert_(names::V_m, 
    	    &balloon_windkessel::get_y_elem_<balloon_windkessel::State_::V_M>);
    insert_(names::t_ref_remaining, 
	    &balloon_windkessel::get_r_);
  }
}

/* ---------------------------------------------------------------- 
 * Iteration function
 * ---------------------------------------------------------------- */

extern "C"
inline int nest::balloon_windkessel_dynamics(double, const double y[], double f[], void* pnode)
{ 
  // a shorthand
  typedef nest::balloon_windkessel::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  const nest::balloon_windkessel& node =  *(reinterpret_cast<nest::balloon_windkessel*>(pnode));

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
    
nest::balloon_windkessel::Parameters_::Parameters_()
  : BW_k    (0.00065     ),  // ms^-1
    BW_gamma(0.00041     ),  // ms^-1
    BW_tau  (980.0       ),  // ms
    BW_alpha(0.32        ),  // dimensionless
    BW_rho  (0.34        ),  // dimensionless
    BW_V0   (0.02        )   // dimensionless
    
{
}

nest::balloon_windkessel::State_::State_(const Parameters_& p)
  : r(0)
{
  //y[V_M] = 0.0; // p.E_L;  // initialize to reversal potential
  y[0] = 0.0;
  y[1] = 1.0;
  y[2] = 1.0;
  y[3] = 1.0;
}

nest::balloon_windkessel::State_::State_(const State_& s)
  : r(s.r)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y[i] = s.y[i];
}

nest::balloon_windkessel::State_& nest::balloon_windkessel::State_::operator=(const State_& s)
{
  if ( this == &s )  // avoid assignment to self
    return *this;

  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y[i] = s.y[i];

  r = s.r;
  return *this;
}

nest::balloon_windkessel::Buffers_::Buffers_(balloon_windkessel& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::balloon_windkessel::Buffers_::Buffers_(const Buffers_&, balloon_windkessel& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::balloon_windkessel::Parameters_::get(DictionaryDatum &d) const
{
  // Dummy function
  def<double>(d,names::BW_k,         BW_k);
}

void nest::balloon_windkessel::Parameters_::set(const DictionaryDatum& d)
{
  // Dummy function
  updateValue<double>(d,names::BW_k,    BW_k);
}


void nest::balloon_windkessel::State_::get(DictionaryDatum &d) const
{
  //def<double>(d, names::V_m, y[V_M]); // Membrane potential
}

void nest::balloon_windkessel::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  //updateValue<double>(d, names::V_m, y[V_M]);
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::balloon_windkessel::balloon_windkessel()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

nest::balloon_windkessel::balloon_windkessel(const balloon_windkessel& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

nest::balloon_windkessel::~balloon_windkessel()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::balloon_windkessel::init_state_(const Node& proto)
{
  const balloon_windkessel& pr = downcast<balloon_windkessel>(proto);
  S_ = pr.S_;
}

void nest::balloon_windkessel::init_buffers_()
{
  Archiving_Node::clear_history();

  B_.spike_exc_.clear();       // includes resize
  B_.spike_inh_.clear();       // includes resize
  B_.currents_.clear();        // includes resize

  B_.logger_.reset();

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

  B_.I_stim_ = 0.0;
}

void nest::balloon_windkessel::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::balloon_windkessel::update(Time const & origin, const long_t from, const long_t to)
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
    S_.y[State_::BW_Z] = B_.spike_exc_.get_value(lag)/Time::get_resolution().get_ms();
    S_.y[State_::BW_Z] += B_.spike_inh_.get_value(lag)/Time::get_resolution().get_ms();

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
    S_.y[State_::V_M] = P_.BW_V0 * (
    					7.0*P_.BW_rho*(1.0-S_.y[State_::BW_Q]) 
    					+ 2.0*(1.0-S_.y[State_::BW_Q]/S_.y[State_::BW_V])
    					+ (2.0*P_.BW_rho-0.2)*(1.0-S_.y[State_::BW_V])
    				   );
    
    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);
  }
}

void nest::balloon_windkessel::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  if(e.get_weight() > 0.0)
    B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			    e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			 e.get_weight() * e.get_multiplicity() );  // inh. spikes are counted just like exc. spikes
}

void nest::balloon_windkessel::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 e.get_weight() * e.get_current());
}

void nest::balloon_windkessel::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

#endif //HAVE_GSL
