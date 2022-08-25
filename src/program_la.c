//   ____                                      
//  |  _ \ _ __ ___   __ _ _ __ __ _ _ __ ___  
//  | |_) | '__/ _ \ / _` | '__/ _` | '_ ` _ \ 
//  |  __/| | | (_) | (_| | | | (_| | | | | | |
//  |_|   |_|  \___/ \__, |_|  \__,_|_| |_| |_|
//                   |___/                     
//   _                _               _                    _ 
//  | |    ___   ___ | | __      __ _| |__   ___  __ _  __| |
//  | |   / _ \ / _ \| |/ /____ / _` | '_ \ / _ \/ _` |/ _` |
//  | |__| (_) | (_) |   <_____| (_| | | | |  __/ (_| | (_| |
//  |_____\___/ \___/|_|\_\     \__,_|_| |_|\___|\__,_|\__,_|

#include "program_la.h"


//   ____            _                 _   _                 
//  |  _ \  ___  ___| | __ _ _ __ __ _| |_(_) ___  _ __  ___ 
//  | | | |/ _ \/ __| |/ _` | '__/ _` | __| |/ _ \| '_ \/ __|
//  | |_| |  __/ (__| | (_| | | | (_| | |_| | (_) | | | \__ \
//  |____/ \___|\___|_|\__,_|_|  \__,_|\__|_|\___/|_| |_|___/
                                                          
// program_la object structure
typedef struct program_la {
  char *filename;                  // file name
  FILE *file;                      // file handle
  block_la_t *first, *last, *current; // block_la pointers
  size_t n;                        // total number of block_las
} program_la_t;


//   _____                 _   _
//  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
//  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
//  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
//  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

// LIFECYCLE ===================================================================

program_la_t *program_la_new(const char *filename) {
  // create memory object
  program_la_t *p = (program_la_t *)calloc(1, sizeof(program_la_t));
  if (!p) {
    perror("Could not create program_la");
    return NULL;
  }
  if (!filename || strlen(filename) == 0) {
    fprintf(stderr, "Improper or empty file name\n");
    return NULL;
  }
  // Initialize fields
  p->filename = strdup(filename);
  p->first = NULL;
  p->last = NULL;
  p->current = NULL;
  p->n = 0;
  return p;
}

// deallocate
void program_la_free(program_la_t *p) {
  assert(p);
  block_la_t *b, *tmp;
  // free the linked list of block_las
  if (p->n > 0) {
    b = p->first;
    do {
      tmp = b;
      b = block_la_next(b);
      block_la_free(tmp);
    } while (b);
  }
  free(p->filename);
  free(p);
  p = NULL;
}

// print a program_la description
void program_la_print(const program_la_t *p, FILE *output) {
  assert(p);
  block_la_t *b = p->first;
  do {
    block_la_print(b, output);
    b = block_la_next(b);
  } while (b);
}


// PROCESSING ==================================================================

// parse the program_la
// return either EXIT_SUCCESS or EXIT_FAILURE
int program_la_parse(program_la_t *p, machine_t *cfg) {
  assert(p && cfg);
  char *line = NULL;
  ssize_t line_len = 0;
  size_t n = 0;
  block_la_t *b;

  // open the file
  p->file = fopen(p->filename, "r");
  if (!p->file) {
    fprintf(stderr, "ERROR: cannot open the file %s\n", p->filename);
    return EXIT_FAILURE;
  }

  // read the file, one line at a time, and create a new block_la for
  // each line
  p->n = 0;
  while ( (line_len = getline(&line, &n, p->file)) >= 0 ) {
    // remove trailing newline (\n) replacing it with a terminator
    if (line[line_len-1] == '\n') {
      line[line_len-1] = '\0'; 
    }
    if(!(b = block_la_new(line, p->last, cfg))) {
      fprintf(stderr, "ERROR: creating the block_la %s\n", line);
      return EXIT_FAILURE;
    }
    if (block_la_parse(b)) {
      fprintf(stderr, "ERROR: parsing the block_la %s\n", line);
      return EXIT_FAILURE;
    }
    if (p->first == NULL) p->first = b;
    p->last = b;
    p->n++;
  }
  fclose(p->file);
  free(line);
  program_la_reset(p);
  return EXIT_SUCCESS;
}

int program_la_look_ahead(program_la_t *p, machine_t *m){
  block_la_t *b = p->first;
  block_la_t *bp;

  eprintf("Computing velocities ...\n");
  while (b){
    block_la_calculate_velocities(b);
    block_la_print_velocity_target(b, stderr);
    b = block_la_next(b);
  }

  eprintf("Computing forward pass ...\n");
  b = p->first;
  while (b){
    block_la_forward_pass(b);
    b = block_la_next(b);
  }

  eprintf("Computing backward pass ...\n");
  b = p->last;
  while (b){
    block_la_backward_pass(b);
    b = block_la_prev(b);
  }

  eprintf("Computing timings ...\n");
  data_t k;
  data_t t = 0;
  data_t tt = 0;
  data_t t_star = 0;
  
  b = p->first;
  while (b){
    if (block_la_compute_raw_profile(b)){
      eprintf("ERROR: in computing timings \n");
      exit(EXIT_FAILURE);
    }
    block_la_print_velocity_profile(b);
    t += block_la_dt(b);

    // If next block is a zero velocity one
    if (!block_la_next(b) || \
        block_la_type(block_la_next(b)) == RAPID || \
        block_la_type(block_la_next(b)) == NO_MOTION ){

      t_star = (size_t)(t / machine_tq(m) + 1) * machine_tq(m);
      k = t_star / t;
      t = 0;

      // Loop back to the previous zero velocity block
      bp = b;
      while(bp && block_la_type(bp) != RAPID && block_la_type(bp) != NO_MOTION){
        block_la_quantize_profile(bp, k);
        bp = block_la_prev(bp);
      } 
    }

    tt += t_star;
    b = block_la_next(b);
  }
  eprintf("Total time: %f\n", tt);
  
  free(b);
  free(bp);
  return 0;
}

// linked-list navigation functions
block_la_t *program_la_next(program_la_t *p) {
  assert(p);
  if (p->current == NULL) p->current = p->first;
  else p->current = block_la_next(p->current);
  return p->current;
}

void program_la_reset(program_la_t *p) {
  assert(p);
  p->current = NULL;
}



// GETTERS =====================================================================

#define program_la_getter(typ, par, name) \
typ program_la_##name(const program_la_t *p) { assert(p); return p->par; }

program_la_getter(char *, filename, filename);
program_la_getter(block_la_t *, first, first);
program_la_getter(block_la_t *, current, current);
program_la_getter(block_la_t *, last, last);
program_la_getter(size_t, n, length);




