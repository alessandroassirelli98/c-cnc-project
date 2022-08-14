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
  data_t a, d;             // acceleration
  data_t f, l;             // nominal feedrate and length
  data_t fi, f1, f2, ff;   // initial and final feedrate
  data_t t1, t2, tf;       // trapezoid times
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
static float block_la_segments_angle(point_t *p0, point_t *p1, point_t *p2);
static float calc_final_velocity(block_la_t *b);
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

    b->prof->f1 = b->prof->f2 = b->prof->fi = b->prof->ff = 0;

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

void block_print_velocity_profile(block_la_t *b, FILE *out){
    assert(b);
    fprintf(out, "%03lu, velocity targets :-> [%f, %f, %f]\n", b->n, b->prof->fi, b->prof->f, b->prof->ff);
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
    b->prof->f = b->feedrate;
    break;
  case ARC_CW:
  case ARC_CCW:
    // calculate arc coordinates
    if (block_la_arc(b)) {
      rv++;
      break;
    }
    // set corrected feedrate and acceleration
    // centripetal acc = f^2/r, must be <= A
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
    b->prof->f = b->feedrate;
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
        data_t fi, ff;
        if(!b->prev){ // if first block
            fi = 0;
            ff = calc_final_velocity(b);
        }
        else if (!b->next){ // if last block
            fi = b->prev->prof->ff;
            ff = 0;
        }
        else{ // otherwise
            fi = b->prev->prof->ff;
            ff = calc_final_velocity(b);
        }
        b->prof->fi = fi;
        b->prof->ff = ff;
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
      fprintf(stderr, "Arc endpoints mismatch error (%f)\n", r - r2);
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
  return 0;
}

// Return a reliable previous point, i.e. machine zero if this is the first 
// block_la
static point_t *point_zero(block_la_t *b) {
  assert(b);
  return b->prev ? b->prev->target : machine_zero(b->machine);
}

// Return the angle alpha of the two segments passing throuh p0, p1, p2
static float block_la_segments_angle(point_t *p0, point_t *p1, point_t *p2){
    data_t theta1, theta2;

    theta1 = atan2(point_y(p1) - point_y(p0), point_x(p1) - point_x(p0));
    theta2 = atan2(point_y(p2) - point_y(p1), point_x(p2) - point_x(p1));
    return (theta2-theta1);
}

static float calc_final_velocity(block_la_t *b){
    data_t alpha;
    alpha = block_la_segments_angle(point_zero(b),b->target,b->next->target);
    return fabs((b->prof->f + b->next->prof->f) / 2 * cos(alpha));
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