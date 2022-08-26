//   ____  _            _    
//  | __ )| | ___   ___| | __
//  |  _ \| |/ _ \ / __| |/ /
//  | |_) | | (_) | (__|   < 
//  |____/|_|\___/ \___|_|\_\
                          
//   _                _               _                    _ 
//  | |    ___   ___ | | __      __ _| |__   ___  __ _  __| |
//  | |   / _ \ / _ \| |/ /____ / _` | '_ \ / _ \/ _` |/ _` |
//  | |__| (_) | (_) |   <_____| (_| | | | |  __/ (_| | (_| |
//  |_____\___/ \___/|_|\_\     \__,_|_| |_|\___|\__,_|\__,_|
//
// Implement here block-related functions for look-ahead


#ifndef BLOCK_LA_H
#define BLOCK_LA_H

#include "defines.h"
#include "point.h"
#include "machine.h"

//   _____                      
//  |_   _|   _ _ __   ___  ___ 
//    | || | | | '_ \ / _ \/ __|
//    | || |_| | |_) |  __/\__ \
//    |_| \__, | .__/ \___||___/
//        |___/|_|              

// Opaque structure representing a G-code block
typedef struct block_la block_la_t;

#define TOL 1e-6

// Block types
typedef enum {
  RAPID = 0,
  LINE,
  ARC_CW,
  ARC_CCW,
  NO_MOTION
} block_la_type_t;

//   _____                 _   _                 
//  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
//  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
//  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
//  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

// LIFECYCLE ===================================================================

block_la_t *block_la_new(const char *line, block_la_t *prev, machine_t *cfg);
void block_la_free(block_la_t *b);
void block_la_print(block_la_t *b, FILE *out);
void block_la_print_velocity_target(block_la_t *b, FILE *out);
int block_la_print_velocity_profile(block_la_t *b, FILE *out);

// ALGORITHMS ==================================================================

// Parsing the G-code string. Returns an integer for success/failure
int block_la_parse(block_la_t *b);

// Compute the tangents to the line at the start and end points
int block_la_compute_tangents(block_la_t *b);

// For the lookahead approach first all the blocks must be parsed
int block_la_compute_velocities(block_la_t *b);

// Computes the forward pass i.e. only accelerations
int block_la_forward_pass(block_la_t *b);

// Computes the backward pass i.e. only decelerations
int block_la_backward_pass(block_la_t *b);

// Compute timings and velocities without taking into account the timesteps
int block_la_compute_raw_profile(block_la_t *b);

// Rescale the block velocity and timings by a factor k
// v* = v/k, t* = k t
int block_la_quantize_profile(block_la_t *b, data_t k);


// Evaluate the value of lambda at a certaint time
// also return speed in the parameter v
data_t block_la_lambda(const block_la_t *b, data_t time, data_t *v);

// Interpolate lambda over three axes
point_t *block_la_interpolate(block_la_t *b, data_t lambda);


// GETTERS =====================================================================

data_t block_la_length(const block_la_t *b);
data_t block_la_dtheta(const block_la_t *b);
data_t block_la_dt(const block_la_t *b);
data_t block_la_r(const block_la_t *b);
block_la_type_t block_la_type(const block_la_t *b);
char *block_la_line(const block_la_t *b);
size_t block_la_n(const block_la_t *b);
point_t *block_la_center(const block_la_t *b);
block_la_t *block_la_next(const block_la_t *b);
block_la_t *block_la_prev(const block_la_t *b);
point_t *block_la_target(const block_la_t *b);


#endif // BLOCK_H