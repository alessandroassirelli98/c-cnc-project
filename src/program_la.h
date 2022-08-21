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
//
// Add here program_la-relatred functions for implementing the look-ahead behaviour


#ifndef program_la_LA_H
#define program_la_LA_H

#include "defines.h"
#include "block_la.h"
#include "machine.h"

//   _____                      
//  |_   _|   _ _ __   ___  ___ 
//    | || | | | '_ \ / _ \/ __|
//    | || |_| | |_) |  __/\__ \
//    |_| \__, | .__/ \___||___/
//        |___/|_|              

// Opaque structure
typedef struct program_la program_la_t;


//   _____                 _   _                 
//  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
//  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
//  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
//  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

// LIFECYCLE ===================================================================

// create a new program_la from the given filename
program_la_t *program_la_new(const char *filename);

// deallocate
void program_la_free(program_la_t *program_la);

// print a program_la description
void program_la_print(const program_la_t *program_la, FILE *output);

// PROCESSING ==================================================================

// parse the program_la
// return either EXIT_SUCCESS or EXIT_FAILURE
int program_la_parse(program_la_t *program_la, machine_t *cfg);

int program_la_look_ahead(program_la_t *p, machine_t *m);

// linked-list navigation functions
block_la_t *program_la_next(program_la_t *program_la);
void program_la_reset(program_la_t *program_la);


// GETTERS =====================================================================

char *program_la_filename(const program_la_t *p);
size_t program_la_length(const program_la_t *p);
block_la_t *program_la_current(const program_la_t *p);
block_la_t *program_la_first(const program_la_t *p);
block_la_t *program_la_last(const program_la_t *p);


#endif // end double inclusion guard