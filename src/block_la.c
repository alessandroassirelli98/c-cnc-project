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
// Implement here block_la-related functions for look-ahead


#include "block_la.h"
#include <ctype.h>

//   ____            _                 _   _
//  |  _ \  ___  ___| | __ _ _ __ __ _| |_(_) ___  _ __  ___
//  | | | |/ _ \/ __| |/ _` | '__/ _` | __| |/ _ \| '_ \/ __|
//  | |_| |  __/ (__| | (_| | | | (_| | |_| | (_) | | | \__ \
//  |____/ \___|\___|_|\__,_|_|  \__,_|\__|_|\___/|_| |_|___/

// Trapezoidal velocity profile
typedef struct {
  data_t l;             // nominal feedrate and length
  point_t *w_s, *w_e; //vector tangent to the circle (if g02/g03) or direction of the line
  data_t vi, v, vf;   // initial and final feedrate
  data_t vi_fwd, vn, vf_fwd, s_inter;     //Only for plottin/debugging prposes
  data_t s1, s2;  //curvilinear absissa of feed transition
  data_t d_t1, d_tm, d_t2;       // trapezoid times
  data_t dt;               // total time and scaling factor
  data_t a1, a2;               // profile acceleration deceleration
  // | Short block | S2 Acc | S1 Acc | S2 Dec | S1 DEC
  // basically depending on the value of s1 and s2 got from the forward and backward pass
  // the bitmask is set. If a 1 appears in one of the columns it means that that pass
  // "has found" a direction to proceed. The backward pass has priority over the forward one, 
  // since it is computed afterwards and only after it we know if the targets are reachable
  int mask;              // bitmask for feed profile
} block_la_profile_t;

// block_la object structure
typedef struct block_la{
  char *line;            // G-code line
  block_la_type_t type;     // type of block_la_la
  size_t n;              // block_la_la number
  size_t tool;           // tool number
  data_t feedrate;       // nominal feedrate
  data_t act_feedrate;   // actual feedrate (possibly reduced along arcs)
  data_t spindle;        // spindle rate
  point_t *target;       // destination point
  point_t *delta;        // distance vector w.r.t. previous point
  point_t *center;       // arc center (if it is an arc)
  data_t length;         // total length
  data_t i, j, r;        // center coordinates and radius (if it is an arc)
  data_t theta0, dtheta; // arc initial angle and arc angle
  data_t acc;            // actual acceleration
  machine_t *machine;    // machine configuration
  block_la_profile_t *prof; // velocity profile
  struct block_la *prev;    // next block_la (linked list)
  struct block_la *next;    // previous block_la
} block_la_t;

// STATIC FUNCTIONS (for internal use only) ====================================
static int block_la_set_fields(block_la_t *b, char cmd, char *arg);
static point_t *point_zero(block_la_t *b);
static int block_la_arc(block_la_t *b);
static float calc_final_velocity(block_la_t *b);
static void calc_s1_s2(block_la_t *b, data_t *s1, data_t *s2, int sign);
static void compute_direction(point_t *v, point_t *p0, point_t *p1);

//   _____                 _   _
//  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
//  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
//  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
//  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

// LIFECYCLE ===================================================================

block_la_t *block_la_new(const char *line, block_la_t *prev, machine_t *cfg) {
  assert(line && cfg); // prev is NULL if this is the first block_la
  block_la_t *b = (block_la_t *)calloc(1, sizeof(block_la_t));
  if (!b) {
    perror("Could not allocate block_la");
    return NULL;
  }

  if (prev) { // copy the memory from the previous block_la
    memcpy(b, prev, sizeof(block_la_t));
    b->prev = prev;
    prev->next = b;
  } else { // this is the first block_la
    // nothing to do
  }

  // non-modal g-code parameters: I, J, R
  b->i = b->j = b->r = 0.0;

  // fields to be calculated
  b->length = 0.0;
  b->target = point_new();
  b->delta = point_new();
  b->center = point_new();

  // allocate memory for profile struct
  b->prof = (block_la_profile_t *)calloc(1, sizeof(block_la_profile_t));
  if (!b->prof) {
    perror("Could not allocate profile structure");
    return NULL;
  }

  b->prof->mask = 0b11111;
  b->prof->w_s = point_new();
  b->prof->w_e = point_new();

  b->machine = cfg;
  b->type = NO_MOTION;
  b->acc = machine_A(b->machine);
  b->line = strdup(line);
  if (! b->line) {
    perror("Could not allocate line");
    return NULL;
  }

  return b;
}

void block_la_free(block_la_t *b) {
  assert(b);
  if (b->line)
    free(b->line);
  if (b->prof){
    free(b->prof->w_s);
    free(b->prof->w_e);
    free(b->prof);
  }
  point_free(b->target);
  point_free(b->center);
  point_free(b->delta);
  free(b);
  b = NULL; 
}

void block_la_print(block_la_t *b, FILE *out) {
  assert(b && out);
  char *start, *end;
  // if this is the first block_la, p0 is the origin
  // otherwise is the target of the previous block_la
  point_t *p0 = point_zero(b);
  // inspect origin and target points
  point_inspect(p0, &start);
  point_inspect(b->target, &end);
  // print out block_la description
  fprintf(out, "%03lu %s->%s F%7.1f S%7.1f T%2lu (G%02d)\n", b->n, start, end, b->feedrate, b->spindle, b->tool, b->type);
  free(end);
  free(start);
}

void block_la_print_velocity_target(block_la_t *b, FILE *out){
    assert(b);
    fprintf(out, "%03lu, velocity targets :-> [%f, %f, %f]\n", b->n, (b->prof->vi)*60, (b->prof->v)*60, (b->prof->vf)*60);
}

// print the velocity target of the block
int block_la_print_velocity_profile(block_la_t *b, const char *out){
  assert(b);
  FILE *f;
  int rv = 0;
  if(!(f = fopen(out, "ab"))){
    fprintf(stderr, "ERROR: cannot open %s to write the space domain velocity profile\n", out);
    rv++;
  }
  fprintf(f, "%03lu, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %f \n", b->n, b->length, b->prof->d_t1, b->prof->d_t2 ,b->prof->vi, b->prof->v,\
                                                      b->prof->vf, b->prof->vn,  b->prof->vi_fwd, b->prof->vf_fwd,\
                                                      b->prof->s_inter, b->prof->s1, b->prof->s2, b->prof->mask,b->acc);
  fclose(f);
  return rv;
}

// ALGORITHMS ==================================================================

// Parsing the G-code string. Returns an integer for success/failure
int block_la_parse(block_la_t *b) {
  assert(b);
  char *word, *line, *tofree;
  point_t *p0;
  int rv = 0;

  tofree = line = strdup(b->line);
  if (!line) {
    perror("Could not allocate momory for tokenizing line");
    return 1;
  }
  // Tokenizing loop
  while ((word = strsep(&line, " ")) != NULL) {
    // word[0] is the command
    // word+1 is the pointer to the argument as a string
    rv += block_la_set_fields(b, toupper(word[0]), word + 1);
  }
  free(tofree);

  // inherit modal fields from the previous block_la
  p0 = point_zero(b);
  point_modal(p0, b->target);
  point_delta(p0, b->target, b->delta);
  b->length = point_dist(p0, b->target);
  
  // deal with motion block_las
  switch (b->type) {
  case LINE:
    // calculate feed profile
    b->acc = machine_A(b->machine);
    b->act_feedrate = b->feedrate;
    b->prof->v = b->prof->vn = (b->feedrate) / 60.0;

    break;
  case ARC_CW:
  case ARC_CCW:
    // calculate arc coordinates
    if (block_la_arc(b)) {
      rv++;
      break;
    }
    // set corrected feedrate and acceleration
    // centripetal acc = v^2/r, must be <= A
    // INI file gives A in mm/s^2, feedrate is given in mm/min
    // We divide by two because, in the critical condition where we have 
    // the maximum feedrate, in the following equation for calculating the 
    // acceleration, it goes to 0. In fact, if we accept the centripetal 
    // acceleration to reach the maximum acceleration, then the tangential 
    // acceleration would go to 0.
    // A more elegant solution would be to calculate a minimum time soltion 
    // for the whole arc, but it is outside the scope.
    b->act_feedrate =
        MIN(b->feedrate, sqrt(machine_A(b->machine)/2.0 * b->r) * 60);
    b->prof->v = (b->act_feedrate) / 60.0;
    // tangential acceleration: when composed with centripetal one, total
    // acceleration must be <= A
    // a^2 <= A^2 + v^4/r^2
    b->acc = sqrt(pow(machine_A(b->machine), 2) - pow(b->act_feedrate / 60, 4) / pow(b->r, 2));
    // deal with complex result
    if (isnan(b->acc)) {
      eprintf("Cannot compute arc: insufficient acceleration");
      rv++;
    }
    // calculate feed profile
    break;
  default:
    break;
  }

  b->prof->l = b->length;
  // return number of parsing errors
  return rv;
}

void block_la_compute_tangents(block_la_t *b){
  data_t theta, ct, st, x, y;
  switch (b->type)
  {
  case LINE:
    // Calculate the direction of the segment and store it
    compute_direction(b->prof->w_s, point_zero(b), b->target);
    compute_direction(b->prof->w_e, point_zero(b), b->target);
    break;

  case ARC_CW:
  case ARC_CCW:
    // Here we get the vector from the starting and endopoint to the center
    // the tangent to the circle is rotated of 90 deg wrt to z axis
    theta = (b->type == ARC_CW) ? M_PI/2 : -M_PI/2;
    ct = cos(theta);
    st = sin(theta);
    
    compute_direction(b->prof->w_s, point_zero(b), b->center);
    compute_direction(b->prof->w_e, b->target, b->center);

    x = point_x(b->prof->w_s);
    y = point_y(b->prof->w_s);
    point_set_x(b->prof->w_s, ct * x - st * y);
    point_set_y(b->prof->w_s, st * x + ct * y);

    x = point_x(b->prof->w_e);
    y = point_y(b->prof->w_e);
    point_set_x(b->prof->w_e, ct * x - st * y);
    point_set_y(b->prof->w_e, st * x + ct * y);
    break;
  default:
    break;
}
}

void block_la_compute_velocities(block_la_t *b){
    assert(b); 
    if ( b->type != RAPID && b->type != NO_MOTION){
        data_t vi, vf;
        if(!b->prev || b->prev->type == RAPID){ // if first block or after a rapid one
            vi = 0;
            vf = calc_final_velocity(b);
        }
        else if (!b->next || b->next->type == RAPID){ // if last block or before a rapid one
            vi = b->prev->prof->vf;
            vf = 0;
        }
        else{ // otherwise
            vi = b->prev->prof->vf;
            vf = calc_final_velocity(b);
        }
        b->prof->vi = vi;
        b->prof->vf = vf;
    }
}

// Computes the forward pass i.e. only accelerations
// Note that here we do not "cover" the deceleration bits,
// since they have higher priority
void block_la_forward_pass(block_la_t *b){
  assert(b);
  if(b->type == RAPID | b->type == NO_MOTION) return;

  data_t v_max;
  data_t s1, s2;

  // Maximum final velocity that can be reached in the current block
  v_max = sqrt(2*b->acc * b->length + pow((b->prof->vi),2));

  if (b->prof->vf > v_max){
    b->prof->vf = v_max;
    b->next->prof->vi = v_max;
  }

  calc_s1_s2(b, &s1, &s2, 1);

  // If it founds s1 positive, then we may start with a acceleration
  if (s1 >= 0){
    b->prof->mask &= 0b11111;
    b->prof->s1 = s1;
  }
  else b->prof->mask &= 0b11011;

  // If it founds s2 < L, then we may end with a deceleration
  if (s2 <= b->length){
    b->prof->mask &= 0b11111; 
    b->prof->s2 = s2;
  }
  else b->prof->mask &= 0b10111;

  b->prof->vi_fwd = b->prof->vi;
  b->prof->vf_fwd = b->prof->vf;

  return;
}

void block_la_backward_pass(block_la_t *b){
  assert(b);
  if(b->type == RAPID | b->type == NO_MOTION) return;

  data_t v_max;
  data_t s1, s2;
  int temp_mask = 0b11111;

  // Check if the final velocity can be reached, 
  // If not, then lower it
  v_max = sqrt(2*b->acc * b->length + pow((b->prof->vf),2));
  if(b->prof->vi > v_max){
    b->prof->vi = v_max;
    b->prev->prof->vf = v_max;
    // Then we have only deceleration
    b->prof->mask &= 0b10011;
  }

  calc_s1_s2(b, &s1, &s2, -1);

  // If it founds s1 positive, then we start with a deceleration
  if (s1 >= 0){
    temp_mask &= 0b11011;
    b->prof->s1 = s1;
  }
  else temp_mask &= 0b11110;

  // If it founds s2 < L, then we end with a deceleration
  if (s2 <= b->length){
    temp_mask &= 0b10111; 
    b->prof->s2 = s2;
  }
  else temp_mask &= 0b11101;
  
  // If the backward pass didn't find any point, then it's all about accelerations
  temp_mask &= (temp_mask == 0b11111) ? 0b11100 : 0b11111; 
  b->prof->mask &= temp_mask;

  // In case of short block
  // This means the velocity is always under the desired one
  if (b->prof->mask == 0b10110 && b->prof->s1 > b->prof->s2){ 
    b->prof->s_inter = b->length/2 + (pow(b->prof->vf, 2) - pow(b->prof->vi, 2))/(4*b->acc);
    
    // Cleaning of the mask: remove useless flag of acceleration or deceleration 
    if (fabs (b->prof->s_inter - b->length) < TOL) b->prof->mask = 0b10100; // If the intersection is in L, then we have only acceleration
    else if (fabs (b->prof->s_inter - 0) < TOL) b->prof->mask = 0b10010; // If it's in zero, then only deceleration

  }
  // This means the velocity is always above the desired one
  else if (b->prof->mask == 0b11001 && b->prof->s1 > b->prof->s2){ //DA block
    b->prof->s_inter = (pow(b->prof->vi, 2) - pow(b->prof->vf, 2))/(4*b->acc) + b->length/2;

    // Cleaning of the mask: remove useless flag of acceleration or deceleration 
    if (fabs (b->prof->s_inter - b->length) < TOL) b->prof->mask = 0b10001; // If the intersection is in L, then we have only deceleration
    else if (fabs (b->prof->s_inter - 0) < TOL) b->prof->mask = 0b11000; // If it's in zero, then only acceleration
  }

  else if (b->prof->s1 < b->length && b->prof->s2 > 0) b->prof->mask &= 0b01111; //otherwise it's not a short block

  // Last check: if the target velocity vi has been changed because it could not be reached,
  // the target vf of the previous block is decreased wrt to the forward pass
  // then the acceleration starting point we computed in the forward pass must be changed in order to
  // comply with the new point.
  calc_s1_s2(b, &s1, &s2, 1);
  b->prof->s2 = (!(b->prof->mask & 0b00010)) ? s2 : b->prof->s2; 
  b->prof->s1 = (!(b->prof->mask & 0b00001)) ? s1 : b->prof->s1;
  return;
}

// Compute timings and velocities without taking into account the timesteps
int block_la_compute_raw_profile(block_la_t *b){
  assert(b);
  // Short block: no maintenance and the first bit of the mask is true
  if(b->type == RAPID || b->type == NO_MOTION) return 0;

  int rv = 0;
  data_t sign, delta, s_check, v, a1, a2;

  data_t vi = b->prof->vi;
  data_t s_int = b->prof->s_inter;
  data_t s1 = b->prof->s1;
  data_t s2 = b->prof->s2;
  data_t l = b->length;
  data_t d_t1 = b->prof->d_t1;
  data_t d_t2 = b->prof->d_t2;
  data_t d_tm = b->prof->d_tm;

  a1 = a2 = b->acc;
  if (b->prof->mask & 0b10000){
    // If the first acceleration bit is set, it means we have to accelerate, otherwise decelerate
    sign = ((b->prof->mask & 0b10101) == 0b10100) ? 1 : -1;
    a1 *= sign;
    // TODO is it the best way? sometimes we have issue with numerical precision,
    // and this is negative, then the sqrt fails
    delta = 2 * a1 * s_int + pow(vi, 2);
    delta = (delta < 0 && delta > - TOL) ? 0 : delta;
    d_t1 = 1/a1 * (-vi + sqrt(delta));
    v= vi + a1 * d_t1;

    // If the second acceleration bit is set, it means we have to accelerate, otherwise decelerate
    sign = ((b->prof->mask & 0b11010) == 0b11000) ? 1 : -1;
    a2 *= sign;
    delta = 2 * a2 * (l - s_int) + pow(v, 2);
    delta = (delta < 0 && delta > - TOL) ? 0 : delta;
    d_t2 = 1/a2 * (-v + sqrt(delta));
  }

  // long block
  else {
    sign = ((b->prof->mask & 0b00101) == 0b00100) ? 1 : -1;
    a1 *= sign;

    delta = 2 * a1 * s1 + pow(vi, 2);
    delta = (delta < 0 && delta > - TOL) ? 0 : delta;
    d_t1 = 1/a1 * (-vi + sqrt(delta));
    
    v = vi +  a1 * d_t1;
    d_tm = (s2 - s1)/v;
    
    sign = ((b->prof->mask & 0b01010) == 0b01000) ? 1 : -1;
    a2 *= sign;

    delta = 2 * a2 * (l - s2) + pow(v, 2);
    delta = (delta < 0 && delta > - TOL) ? 0 : delta;
    d_t2 = 1/a2 * (-v + sqrt(delta));
  }

  s1 = vi * d_t1 + a1 / 2.0 * pow(d_t1, 2);
  v = vi + a1 * d_t1;
  s2 = s1 + v * d_tm;
  s_check = s2 + v * d_t2 + a2 / 2 * pow(d_t2, 2);

  if (isnan(d_t1)){
    fprintf(stderr, "ERROR : cannot compute t1 for block: %03lu\n", b->n);
    rv++;
    }

  if (isnan(d_t2)){
    fprintf(stderr, "ERROR : cannot compute  t2 for block: %03lu\n", b->n);
    rv++;

    }

  if (fabs(b->length - s_check) > TOL){
    fprintf(stderr, "Mistmach of interpolated and nominal length %03lu Error is : %f\n", b->n, b->length - s_check);
    rv++;
  }

  b->prof->a1 = a1;
  b->prof->a2 = a2;
  b->prof->d_t1 = d_t1;
  b->prof->d_t2 = d_t2;
  b->prof->d_tm = d_tm;
  b->prof->dt = d_t1 + d_t2 + d_tm;
  b->prof->v = v;


  return rv;
}

// Rescale the block velocity and timings by a factor k
// v* = v/k, t* = k t
void block_la_quantize_profile(block_la_t *b, data_t k){

  b->prof->d_t1 *= k;
  b->prof->d_t2 *= k;
  b->prof->d_tm *= k;
  b->prof->dt *= k;

  b->prof->vi /= k;
  b->prof->v /=k;
  b->prof->vf /= k;

  b->prof->a1 /= pow(k, 2);
  b->prof->a2 /= pow(k, 2);

}

// Evaluate the value of lambda at a certaint time
data_t block_la_lambda(const block_la_t *b, data_t t, data_t *v) {
  assert(b);
  data_t s1, s2, v1;
  data_t r;
  data_t d_t1 = b->prof->d_t1;
  data_t d_t2 = b->prof->d_t2;
  data_t d_tm = b->prof->d_tm;
  data_t a1 = b->prof->a1;
  data_t a2 = b->prof->a2;
  data_t vi = b->prof->vi;

  if (t < 0) {
    r = 0.0;
    *v = 0.0;
  }
  else if (t < d_t1) { // acceleration
    r = vi * t + a1 * pow(t, 2) / 2.0;
    *v = vi +a1 * t;
  }
  else if (t < (d_t1 + d_tm)) { // maintenance
    s1 = vi * d_t1 + a1 / 2.0 * pow(d_t1, 2);
    v1 = vi + a1 * d_t1;

    r = s1 + v1*(t - d_t1);
    *v = v1;
  }
  else if (t < (d_t1 + d_tm + d_t2)) { // deceleration
    data_t t_2 = d_t1 + d_tm;

    s1 = vi * d_t1 + a1 / 2.0 * pow(d_t1,2);
    v1 = vi + a1 * d_t1;
    s2 = s1 + v1 * d_tm;

    r =  s2 + v1 * (t - t_2) + a2 / 2 * pow(t - t_2, 2);
    *v = v1 + a2 * (t - t_2);
  }

  else {
    r = b->prof->l;
    *v = b->prof->vf;
  }

  r /= b->prof->l;
  *v *= 60; // convert to mm/min

  return r;
}

// CAREFUL: this function allocates a point
point_t *block_la_interpolate(block_la_t *b, data_t lambda) {
  assert(b);
  point_t *result = machine_setpoint(b->machine);
  point_t *p0 = point_zero(b);

  if (b->type == LINE) {
    point_set_x(result, point_x(p0) + point_x(b->delta) * lambda);
    point_set_y(result, point_y(p0) + point_y(b->delta) * lambda);
  }
  else if (b->type == ARC_CW || b->type == ARC_CCW) {
    point_set_x(result, point_x(b->center) + 
      b->r * cos(b->theta0 + b->dtheta * lambda));
    point_set_y(result, point_y(b->center) + 
      b->r * sin(b->theta0 + b->dtheta * lambda));
  }
  else {
    fprintf(stderr, "Unexpected block type!\n");
    return NULL;
  }
  point_set_z(result, point_z(p0) + point_z(b->delta) * lambda);

  return result;
}


// GETTERS =====================================================================

#define block_la_getter(typ, par, name) \
typ block_la_##name(const block_la_t *b) { assert(b); return b->par; }

block_la_getter(data_t, length, length);
block_la_getter(data_t, dtheta, dtheta);
block_la_getter(data_t, prof->dt, dt);
block_la_getter(block_la_type_t, type, type);
block_la_getter(char *, line, line);
block_la_getter(size_t, n, n);
block_la_getter(data_t, r, r);
block_la_getter(point_t *, center, center);
block_la_getter(block_la_t *, next, next);
block_la_getter(block_la_t *, prev, prev);
block_la_getter(point_t *, target, target);
block_la_getter(data_t, act_feedrate, act_feedrate);


 

//   ____  _        _   _         __                  
//  / ___|| |_ __ _| |_(_) ___   / _|_   _ _ __   ___ 
//  \___ \| __/ _` | __| |/ __| | |_| | | | '_ \ / __|
//   ___) | || (_| | |_| | (__  |  _| |_| | | | | (__ 
//  |____/ \__\__,_|\__|_|\___| |_|  \__,_|_| |_|\___|
// Definitions for the static functions declared above

// Calculate the arc coordinates
static int block_la_arc(block_la_t *b) {
  data_t x0, y0, z0, xc, yc, xf, yf, zf, r;
  point_t *p0 = point_zero(b);
  x0 = point_x(p0);
  y0 = point_y(p0);
  z0 = point_z(p0);
  xf = point_x(b->target);
  yf = point_y(b->target);
  zf = point_z(b->target);

  if (b->r) { // if the radius is given
    data_t dx = point_x(b->delta);
    data_t dy = point_y(b->delta);
    r = b->r;
    data_t dxy2 = pow(dx, 2) + pow(dy, 2);
    data_t sq = sqrt(-pow(dy, 2) * dxy2 * (dxy2 - 4 * r * r));
    // signs table
    // sign(r) | CW(-1) | CCW(+1)
    // --------------------------
    //      -1 |     +  |    -
    //      +1 |     -  |    +
    int s = (r > 0) - (r < 0);
    s *= (b->type == ARC_CCW ? 1 : -1);
    xc = x0 + (dx - s * sq / dxy2) / 2.0;
    yc = y0 + dy / 2.0 + s * (dx * sq) / (2 * dy * dxy2);
  }
  else { // if I,J are given
    data_t r2;
    r = hypot(b->i, b->j);
    xc = x0 + b->i;
    yc = y0 + b->j;
    r2 = hypot(xf - xc, yf - yc);
    if (fabs(r - r2) > machine_error(b->machine)) {
      fprintf(stderr, "Arc endpoints mismatch error (%v)\n", r - r2);
      return 1;
    }
    b->r = r;
  }
  point_set_x(b->center, xc);
  point_set_y(b->center, yc);
  b->theta0 = atan2(y0 - yc, x0 - xc);
  b->dtheta = atan2(yf - yc, xf - xc) - b->theta0;
  // we need the net angle so we take the 2PI complement if negative
  if (b->dtheta <0) 
    b->dtheta = 2 * M_PI + b->dtheta;
  // if CW, take the negative complement
  if (b->type == ARC_CW)
    b->dtheta = -(2 * M_PI - b->dtheta);
  //
  b->length = hypot(zf - z0, b->dtheta * b->r);
  // from now on , it's safer to drop radius angle
  b->r = fabs(b->r);
  b->length = hypot(zf - z0, b->dtheta * b->r);

  return 0;
}

// Return a reliable previous point, i.e. machine zero if this is the first 
// block_la
static point_t *point_zero(block_la_t *b) {
  assert(b);
  return b->prev ? b->prev->target : machine_zero(b->machine);
}

// Compute the terminal velocity of the current block b
// based on the nominal velocity of the current block and the next one
// and also on the angle between the two segments 
static float calc_final_velocity(block_la_t *b){
  assert(b);

    data_t dp, ctheta, m0, m1, v;
    data_t x0 = point_x(b->prof->w_e);
    data_t y0 = point_y(b->prof->w_e);
    data_t z0 = point_z(b->prof->w_e);
    data_t x1 = point_x(b->next->prof->w_s);
    data_t y1 = point_y(b->next->prof->w_s);
    data_t z1 = point_z(b->next->prof->w_s);

    // Dot product of the current tangent vector and the previous one
    dp = x0 * x1 + y0 * y1 + z0 * z1;
    // Magnitude of the two vectors
    m0 = sqrt(pow(x0, 2) + pow(y0, 2) + pow(z0, 2));
    m1 = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));

    // Calculate cos(theta) as dot_prod(a, b) / |ab|
    ctheta = (dp / (m0 * m1)); 

    v = fabs(((b->prof->v + b->next->prof->v) / 2) * ctheta);
    v =  (ctheta >= 1/2) ? v : 0; // If the angle is > 60 deg, the the target is zero vel

    return v;
}

// Compute s1 and s2 of the block, sign should be 1 for forward, -1 for backward
static void calc_s1_s2(block_la_t *b, data_t *s1, data_t *s2, int sign){
  assert(b && s1 && s2);
  *s1 = (pow(b->prof->v, 2) - pow(b->prof->vi, 2)) / (2 * sign *b->acc);
  *s2 = b->length + (pow(b->prof->v, 2) - pow(b->prof->vf, 2)) / (2 * sign * b->acc);
}

// Calculate the components of the vector p1-p0
static void compute_direction(point_t *v, point_t *p0, point_t *p1){
  data_t i, j, k;
  i = point_x(p1) - point_x(p0);
  j = point_y(p1) - point_y(p0);
  k = point_z(p1) - point_z(p0);

  point_set_xyz(v, i, j, k);
}

// Parse a single G-code word (cmd+arg)
static int block_la_set_fields(block_la_t *b, char cmd, char *arg) {
  assert(b && arg);
  switch (cmd)
  {
  case 'N':
    b->n = atol(arg);
    break;
  case 'G':
    b->type = (block_la_type_t)atoi(arg);
    break;
  case 'X':
    point_set_x(b->target, atof(arg));
    break;
  case 'Y':
    point_set_y(b->target, atof(arg));
    break;
  case 'Z':
    point_set_z(b->target, atof(arg));
    break;
  case 'I': 
    b->i = atof(arg);
    break;
  case 'J':
    b->j = atof(arg);
    break;
  case 'R':
    b->r = atof(arg);
    break;
  case 'F':
    b->feedrate = atof(arg);
    break;
  case 'S':
    b->spindle = atof(arg);
    break; 
  case 'T':
    b->tool = atol(arg);   
    break;
  default:
    fprintf(stderr, "ERROR: Usupported G-code command %c%s\n", cmd, arg);
    return 1;
    break;
  }
  // cannot have R and IJ on the same block_la
  if (b->r && (b->i || b->j)) {
    fprintf(stderr, "ERROR: Cannot mix R and IJ\n");
    return 1;
  }
  return 0;
}




//   _____ _____ ____ _____   __  __       _       
//  |_   _| ____/ ___|_   _| |  \/  | __ _(_)_ __  
//    | | |  _| \___ \ | |   | |\/| |/ _` | | '_ \
//    | | | |___ ___) || |   | |  | | (_| | | | | |
//    |_| |_____|____/ |_|   |_|  |_|\__,_|_|_| |_|
//
#ifdef BLOCK_LA_MAIN
int main() {
  block_la_t *b1 = NULL, *b2 = NULL, *b3 = NULL;
  machine_t *cfg = machine_new(NULL);

  b1 = block_la_new("N10 G00 X90 Y90 Z100 t3", NULL, cfg);
  block_la_parse(b1);
  b2 = block_la_new("N20 G01 Y100 X100 F1000 S2000", b1, cfg);
  block_la_parse(b2);
  b3 = block_la_new("N30 G01 Y200", b2, cfg);
  block_la_parse(b3);

  block_la_print(b1, stdout);
  block_la_print(b2, stdout);
  block_la_print(b3, stdout);

  block_la_free(b1);
  block_la_free(b2);
  block_la_free(b3);
  return 0;
}
#endif