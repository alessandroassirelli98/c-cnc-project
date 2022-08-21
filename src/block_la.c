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
  data_t vi, v1, v, v_inter, v2, vf;   // initial and final feedrate
  data_t vi_fwd, vf_fwd;     //Only for debug prposes
  data_t s1, s2, s_inter;  //curvilinear absissa of feed transition
  int mask;              // bitmask for feed profile
  data_t d_t1, d_tm, d_t2;       // trapezoid times
  data_t dt;               // total time
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
  data_t alpha_s, alpha_e; //angle of the tangent to the circle (if g02/g03) or angle of the line
  data_t acc;            // actual acceleration
  machine_t *machine;    // machine configuration
  block_la_profile_t *prof; // velocity profile
  struct block_la *prev;    // next block_la (linked list)
  struct block_la *next;    // previous block_la
} block_la_t;

// STATIC FUNCTIONS (for internal use only) ====================================
static int block_la_set_fields(block_la_t *b, char cmd, char *arg);
static point_t *point_zero(block_la_t *b);
static void block_la_compute(block_la_t *b);
static int block_la_arc(block_la_t *b);
static float calc_final_velocity(block_la_t *b);
static void calc_s1_s2(block_la_t *b, data_t *s1, data_t *s2, int sign);
static data_t quantize(data_t t, data_t tq, data_t *dq);

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

  //b->prof->v1 = b->prof->v2 = b->prof->vi = b->prof->vf = 0;
  b->prof->mask = 0b11111;

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
  if (b->prof)
    free(b->prof);
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

void block_la_print_velocity_profile(block_la_t *b){
    assert(b);
    printf("%03lu, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %f \n", b->n, b->length, b->prof->vi, b->prof->v,\
                                                       b->prof->vf, b->prof->vi_fwd, b->prof->vf_fwd, b->prof->v1, b->prof->v_inter, b->prof->v2,\
                                                       b->prof->s_inter, b->prof->s1, b->prof->s2, b->prof->mask,b->acc);
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
    b->prof->v = (b->feedrate) / 60.0;
    // Calculate the angle of the segment
    b->alpha_s = b->alpha_e = \
              atan2(point_y(b->target) - point_y(point_zero(b)), point_x(b->target) - point_x(point_zero(b)));
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
  // return number of parsing errors
  return rv;
}

// For the lookahead approach first all the blocks must be parsed
int block_la_calculate_velocities(block_la_t *b){
    assert(b); 
    if ( b->type != RAPID ){
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
    return 0;
}


int block_la_forward_pass(block_la_t *b){
  assert(b);
  if(b->type == RAPID) return 0;

  data_t v_max;
  data_t *s1, *s2;
  s1 = calloc(1, sizeof(data_t));
  s2 = calloc(1, sizeof(data_t));

  // Maximum final velocity that can be reached in the current block
  v_max = sqrt(2*b->acc * b->length + pow((b->prof->vi),2));

  if (b->prof->vf > v_max){
    b->prof->vf = v_max;
    b->next->prof->vi = v_max;
  }

  calc_s1_s2(b, s1, s2, 1);

  if (*s1 >= 0){
    b->prof->mask &= 0b11111;
    b->prof->s1 = *s1;
  }
  else b->prof->mask &= 0b11011;

  if (*s2 <= b->length){
    b->prof->mask &= 0b11111; 
    b->prof->s2 = *s2;
  }
  else b->prof->mask &= 0b10111;

  b->prof->vi_fwd = b->prof->vi;
  b->prof->vf_fwd = b->prof->vf;


  free(s1);
  free(s2);
  return 0;
}

int block_la_backward_pass(block_la_t *b){
  assert(b);
  if(b->type == RAPID) return 0;

  data_t v_max;
  data_t *s1, *s2;
  int temp_mask = 0b11111;

  s1 = calloc(1, sizeof(data_t));
  s2 = calloc(1, sizeof(data_t));

  // Check if the final velocity can be reached, 
  // If not, then lower it
  v_max = sqrt(2*b->acc * b->length + pow((b->prof->vf),2));
  if(b->prof->vi > v_max){
    b->prof->vi = v_max;
    b->prev->prof->vf = v_max;
    b->prof->mask &= 0b10011;
  }

  calc_s1_s2(b, s1, s2, -1);

  // If it founds s1 positive, then we start with a deceleration
  if (*s1 >= 0){
    temp_mask &= 0b11011;
    b->prof->s1 = *s1;
  }
  else temp_mask &= 0b11110;

  // If it founds s2 < L, then we end with a deceleration
  if (*s2 <= b->length){
    temp_mask &= 0b10111; 
    b->prof->s2 = *s2;
  }
  else temp_mask &= 0b11101;
  
  // If the backward pass didn't find any point, then it's all abaut acceleration
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
  calc_s1_s2(b, s1, s2, 1);
  b->prof->s2 = (!(b->prof->mask & 0b00010)) ? *s2 : b->prof->s2; 
  b->prof->s1 = (!(b->prof->mask & 0b00001)) ? *s1 : b->prof->s1;

  free(s1);
  free(s2);
  return 0;
}


// Compute timings and velocities without taking into account the timesteps
int block_la_compute_timings(block_la_t *b){
  assert(b);
  data_t sign, v_check, dq;
  // Short block: no maintenance and the first bit of the mask is true
  if(b->type == RAPID) return 0;
  if (b->prof->mask & 0b10000){
    // If the first acceleration bit is set, it means we have to accelerate, otherwise decelerate
    sign = ((b->prof->mask & 0b10101) == 0b10100) ? 1 : -1;
    //plus_minus = (b->prof->s_inter == 0) ? 1 : 1;

    b->prof->d_t1 = 1/(sign * b->acc) * (-b->prof->vi + sqrt(2 * sign * b->acc * b->prof->s_inter + pow(b->prof->vi, 2)));
    b->prof->v_inter = b->prof->vi + sign * b->acc * b->prof->d_t1;

    // If the second acceleration bit is set, it means we have to accelerate, otherwise decelerate
    sign = ((b->prof->mask & 0b11010) == 0b11000) ? 1 : -1;
    b->prof->d_t2 = 1/(sign * b->acc) * (-b->prof->v_inter + sqrt(2 * sign * b->acc * (b->length - b->prof->s_inter) + pow(b->prof->v_inter, 2)));
    v_check = b->prof->v_inter + sign * b->acc * (b->prof->d_t2);

    if (fabs(v_check - b->prof->vf) > TOL){
      fprintf(stderr, "Error on block %03lu v_final computed is different from the one assigned\
      the error is: %f \n", b->n, (b->prof->vf - v_check));
      //exit(EXIT_FAILURE);
    }



  }

  // long block
  else{
    data_t v, a, d;
    sign = ((b->prof->mask & 0b00101) == 0b00100) ? 1 : -1;
    b->prof->d_t1 = 1/(sign * b->acc) * (-b->prof->vi + sqrt(2 * sign * b->acc * b->prof->s1 + pow(b->prof->vi, 2)));
    b->prof->v1 = b->prof->vi + sign * b->acc * b->prof->d_t1;

    b->prof->d_tm = (b->prof->s2 - b->prof->s1)/b->prof->v;
    b->prof->v2 = b->prof->v;

    sign = ((b->prof->mask & 0b01010) == 0b01000) ? 1 : -1;
    b->prof->d_t2 = 1/(sign * b->acc) * (-b->prof->v + sqrt(2 * sign * b->acc * (b->length - b->prof->s2)+ pow(b->prof->v, 2)));
    v_check = b->prof->v + sign * b->acc * (b->prof->d_t2);
    if (fabs(v_check - b->prof->vf) > TOL){
      fprintf(stderr, "Error on block %03lu v_final computed is different from the one assigned\
      the error is: %f\n", b->n, (b->prof->vf - v_check));
      //exit(EXIT_FAILURE);
    }

    v = (2*b->length - b->prof->vi * b->prof->d_t1 - b->prof->vf * b->prof->d_t2)/(b->prof->d_t1 + b->prof->d_t2 + 2);

  }

  return 0;
}

// Evaluate the value of lambda at a certaint time
data_t block_la_lambda(const block_la_t *b, data_t t, data_t *v) {

}

// CAREFUL: this function allocates a point
point_t *block_la_interpolate(block_la_t *b, data_t lambda) {
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

// Calculate the integer multiple of sampling time; also prvide the rounding
// amount in dq
static data_t quantize(data_t t, data_t tq, data_t *dq) {
  data_t q;
  q = ((size_t)(t / tq) + 1) * tq;
  *dq = q - t;
  return q;
}

// Calcultare the velocity profile
static void block_la_compute(block_la_t *b) {
  
}

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

  // Calculate the tangent to the arc in the initial point
  b->alpha_s = M_PI/2 + atan2(y0 - yc, x0 - xc);

  // Calculate the tangent to the arc in the final point
  b->alpha_e = M_PI/2 + atan2(yf - yc, xf - xc);
  return 0;
}

// Return a reliable previous point, i.e. machine zero if this is the first 
// block_la
static point_t *point_zero(block_la_t *b) {
  assert(b);
  return b->prev ? b->prev->target : machine_zero(b->machine);
}

// Compute the terminal velocity of the current block b
static float calc_final_velocity(block_la_t *b){
    data_t alpha;
    alpha = b->next->alpha_s - b->alpha_e;
    return fabs((b->prof->v + b->next->prof->v) / 2 * cos(alpha));
}

// Compute s1 and s2 of the block, sign should be 1 for forward, -1 for backward
static void calc_s1_s2(block_la_t *b, data_t *s1, data_t *s2, int sign){
  assert(s1 && s2);

  *s1 = (pow(b->prof->v, 2) - pow(b->prof->vi, 2)) / (2 * sign *b->acc);
  *s2 = b->length + (pow(b->prof->v, 2) - pow(b->prof->vf, 2)) / (2 * sign * b->acc);
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