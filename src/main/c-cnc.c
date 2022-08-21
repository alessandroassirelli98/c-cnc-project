#include "../defines.h"
#include "../machine.h"
#include "../program_la.h"
#include "../block_la.h"
#include "../point.h"
//#include "../fsm.h"

#define eprintf(...) fprintf(stderr, __VA_ARGS__)

#if 0
int main(int argc, char const *argv[]) {
  ccnc_state_data_t state_data = {
    .ini_file = "settings.ini",
    .prog_file = argv[1],
    .machine = NULL,
    .prog = NULL
  };
  ccnc_state_t cur_state = CCNC_STATE_INIT;
  do {
    cur_state = ccnc_run_state(cur_state, &state_data);
    wait_next(machine_tq(state_data.machine) * 1E9 / machine_rt_pacing(state_data.machine));
  } while (cur_state != CCNC_STATE_STOP);
  ccnc_run_state(cur_state, &state_data);
  return 0;
}



#else
int main(int argc, char const *argv[]) {

  point_t *sp = NULL;
  block_la_t *b = NULL;
  program_la_t *p = NULL;
  data_t t, tt, tq, lambda, f;
  machine_t *machine = machine_new("settings.ini");
  if (!machine) {
    eprintf("Error creating machine instance\n");
    exit(EXIT_FAILURE);
  }
  tq = machine_tq(machine);
  p = program_la_new(argv[1]);
  if (!p) {
    eprintf("Could not create program, exiting.\n");
    exit(EXIT_FAILURE);
  }
  if (program_la_parse(p, machine) == EXIT_FAILURE) {
    eprintf("Could not parse program in %s, exiting.\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  program_la_print(p, stderr);
  program_la_look_ahead(p, machine);


  machine_free(machine);
  return 0;
}

#endif