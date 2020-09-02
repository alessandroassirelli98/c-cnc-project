/******************************************************************************
Finite State Machine
Project: C-CNC 2020
Description: Finite state machine for C-CNC

Generated by gv_fsm ruby gem, see https://rubygems.org/gems/gv_fsm
gv_fsm version 0.2.3
Generation date: 2020-09-01 17:39:38 +0200
Generated from: fsm2.dot
The finite state machine has:
  4 states
  2 transition functions
******************************************************************************/

#ifndef FSM_H
#define FSM_H
#include <stdlib.h>
#include <limits.h>
#include <libgen.h> // for dirname()
#include "ccnc.h"
#include "program.h"

// State data object
// By default set to void; override this typedef or load the proper
// header if you need
typedef struct {
  struct machine *m;
  program_t *p;
  struct machine_config *cfg;
  char this_dir[PATH_MAX];    // PATH MAX is the maximum path length
  char config_path[PATH_MAX]; // for the current operating system
  char viewer_path[PATH_MAX];
  char program_file[PATH_MAX];
} state_data_t;

// NOTHING SHALL BE CHANGED AFTER THIS LINE!

// List of states
typedef enum {
  STATE_INIT = 0,  
  STATE_IDLE,  
  STATE_RUN,  
  STATE_STOP,  
  NUM_STATES,
  NO_CHANGE
} state_t;

// State human-readable names
extern const char *state_names[];

// State function and state transition prototypes
typedef state_t state_func_t(state_data_t *data);
typedef void transition_func_t(state_data_t *data);

// State functions

// Function to be executed in state init
// valid return states: STATE_IDLE
state_t do_init(state_data_t *data);

// Function to be executed in state idle
// valid return states: NO_CHANGE, STATE_IDLE, STATE_RUN, STATE_STOP
state_t do_idle(state_data_t *data);

// Function to be executed in state run
// valid return states: STATE_IDLE
state_t do_run(state_data_t *data);

// Function to be executed in state stop
// valid return states: NO_CHANGE
state_t do_stop(state_data_t *data);


// List of state functions
extern state_func_t *const state_table[NUM_STATES];


// Transition functions
void setup(state_data_t *data);
void idle_to_stop(state_data_t *data);

// Table of transition functions
extern transition_func_t *const transition_table[NUM_STATES][NUM_STATES];

// state manager
state_t run_state(state_t cur_state, state_data_t *data);

#endif
