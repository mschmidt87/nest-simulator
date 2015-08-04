/*
 *  bold_detector.h
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

#ifndef BOLD_DETECTOR_H
#define BOLD_DETECTOR_H


#include <vector>
#include "nest.h"
#include "event.h"
#include "node.h"
#include "recording_device.h"
#include "exceptions.h"
#include "ring_buffer.h"

#ifdef HAVE_GSL

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* BeginDocumentation

Name: bold_detector - Device for detecting single spikes.

Description:
The bold_detector device is a recording device. It is used to record
spikes from a single neuron, or from multiple neurons at once. Data
is recorded in memory or to file as for all RecordingDevices.
By default, GID and time of each spike is recorded.

The spike detector can also record spike times with full precision
from neurons emitting precisely timed spikes. Set /precise_times to
achieve this.

Any node from which spikes are to be recorded, must be connected to
the spike detector using a normal connect command. Any connection weight
and delay will be ignored for that connection.

Simulations progress in cycles defined by the minimum delay. During each
cycle, the spike detector records (stores in memory or writes to screen/file)
the spikes generated during the previous cycle. As a consequence, any
spikes generated during the cycle immediately preceding the end of the simulation
time will not be recorded. Setting the /stop parameter to at the latest one
min_delay period before the end of the simulation time ensures that all spikes
desired to be recorded, are recorded.

Spike are not necessarily written to file in chronological order.

Receives: SpikeEvent

SeeAlso: bold_detector, Device, RecordingDevice
*/


namespace nest
{

class Network;

/**
   * Function computing right-hand side of ODE for GSL solver.
   * @note Must be declared here so we can befriend it in class.
   * @note Must have C-linkage for passing to GSL. Internally, it is
   *       a first-class C++ function, but cannot be a member function
   *       because of the C-linkage.
   * @note No point in declaring it inline, since it is called
   *       through a function pointer.
   * @param void* Pointer to model neuron instance.
   */
  extern "C"
  int balloon_windkessel_dynamics (double, const double*, double*, void*);

/**
 * Spike detector class.
 *
 * This class manages spike recording for normal and precise spikes. It
 * receives spikes via its handle(SpikeEvent&) method, buffers them, and
 * stores them via its RecordingDevice in the update() method.
 *
 * Spikes are buffered in a two-segment buffer. We need to distinguish between
 * two types of spikes: those delivered from the global event queue (almost all
 * spikes) and spikes delivered locally from devices that are replicated on VPs
 * (has_proxies() == false).
 * - Spikes from the global queue are delivered by deliver_events() at the
 *   beginning of each update cycle and are stored only until update() is called
 *   during the same update cycle. Global queue spikes are thus written to the
 *   read_toggle() segment of the buffer, from which update() reads.
 * - Spikes delivered locally may be delivered before or after
 *   bold_detector::update() is executed. These spikes are therefore buffered in
 *   the write_toggle() segment of the buffer and output during the next cycle.
 * - After all spikes are recorded, update() clears the read_toggle() segment
 *   of the buffer.
 *
 * @ingroup Devices
 */
class bold_detector : public Node
{

public:
  bold_detector();
  bold_detector( const bold_detector& );

  void set_has_proxies( const bool hp );
  bool
  has_proxies() const
  {
    return has_proxies_;
  }
  bool
  potential_global_receiver() const
  {
    return true;
  }
  void set_local_receiver( const bool lr );
  bool
  local_receiver() const
  {
    return local_receiver_;
  }

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  void handle( SpikeEvent& );

  port handles_test_event( SpikeEvent&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
//! Model parameters
  struct Parameters_ {
    // Balloon-Windkessel parameters
    double_t BW_k;
    double_t BW_gamma;
    double_t BW_tau;
    double_t BW_alpha;
    double_t BW_rho;
    double_t BW_V0;

    // store number of incoming connections
    long_t N_conn;
  
    Parameters_();        //!< Set default parameter values

    void get(DictionaryDatum&) const;  //!< Store current values in dictionary
    void set(const DictionaryDatum&);  //!< Set values from dicitonary
  };


public:
  struct State_ {
      
    //! Symbolic indices to the elements of the state vector y
    enum StateVecElems { BW_S, BW_F, BW_V, BW_Q , BOLD, BW_Z ,
			 STATE_VEC_SIZE };

    //! state vector, must be C-array for GSL solver
    double_t y[STATE_VEC_SIZE];
  
    //!< number of refractory steps remaining
    int_t r;

    State_(const Parameters_&);  //!< Default initialization
    State_(const State_&);
    State_& operator=(const State_&);

    void get(DictionaryDatum&) const;  //!< Store current values in dictionary

    /**
     * Set state from values in dictionary.
     * Requires Parameters_ as argument to, eg, check bounds.'
     */
    void set(const DictionaryDatum&, const Parameters_&);
  };    

private:
  void init_state_( Node const& );
  void init_buffers_();
  void calibrate();
  void finalize();

  /**
   * Update detector by recording spikes.
   *
   * All spikes in the read_toggle() half of the spike buffer are
   * recorded by passing them to the RecordingDevice, which then
   * stores them in memory or outputs them as desired.
   *
   * @see RecordingDevice
   */
  void update( Time const&, const long_t, const long_t );

  // make dynamics function quasi-member
  friend int balloon_windkessel_dynamics(double, const double*, double*, void*);


  /**
   * Buffer for incoming spikes.
   *
   * This data structure buffers all incoming spikes until they are
   * passed to the RecordingDevice for storage or output during update().
   * update() always reads from spikes_[network()->read_toggle()] and
   * deletes all events that have been read.
   *
   * Events arriving from locally sending nodes, i.e., devices without
   * proxies, are stored in spikes_[network()->write_toggle()], to ensure
   * order-independent results.
   *
   * Events arriving from globally sending nodes are delivered from the
   * global event queue by Scheduler::deliver_events() at the beginning
   * of the time slice. They are therefore written to spikes_[network()->read_toggle()]
   * so that they can be recorded by the subsequent call to update().
   * This does not violate order-independence, since all spikes are delivered
   * from the global queue before any node is updated.
   */

  struct Buffers_
  {
    Buffers_(bold_detector&); //!<Sets buffer pointers to 0
    Buffers_(const Buffers_&, bold_detector&); //!<Sets buffer pointers to 0

    RingBuffer spikes_;

    /* GSL ODE stuff */
    gsl_odeiv_step*    s_;    //!< stepping function
    gsl_odeiv_control* c_;    //!< adaptive stepsize control function
    gsl_odeiv_evolve*  e_;    //!< evolution function
    gsl_odeiv_system   sys_;  //!< struct describing system
      
    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double_t step_;           //!< step size in ms
    double   IntegrationStep_;//!< current integration time step, updated by GSL

  };

  RecordingDevice device_;
  // keep the order of these lines, seems to give best performance
  Parameters_ P_;
  State_      S_;
  Buffers_    B_;

  bool user_set_precise_times_;
  bool has_proxies_;
  bool local_receiver_;




};

inline void
bold_detector::set_has_proxies( const bool hp )
{
  has_proxies_ = hp;
}

inline void
bold_detector::set_local_receiver( const bool lr )
{
  local_receiver_ = lr;
}

inline port
bold_detector::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );

  Parameters_ ptmp = P_;
  ptmp.N_conn += 1;
  P_ = ptmp;
  return 0;
}

inline void
bold_detector::finalize()
{
  device_.finalize();
}

} // namespace

#endif /* #ifndef BOLD_DETECTOR_H */

#endif //HAVE_GSL
