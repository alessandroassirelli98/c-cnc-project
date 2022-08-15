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

  block_la_t *b1 = NULL, *b2 = NULL, *b3 = NULL;
  machine_t *machine = machine_new("settings.ini");

  b1 = block_la_new("N10 G00 X0 Y0 Z0", NULL, machine);
  block_la_parse(b1);
  b2 = block_la_new("N20 G1 X100 Y0 Z0 F1000 S2000", b1, machine);
  block_la_parse(b2);
  b3 = block_la_new("N30 G02 X100 Y50 I0 J25", b2, machine);
  block_la_parse(b3);



  block_la_print(b1, stdout);
  block_la_calculate_velocities(b1);
  block_print_velocity_profile(b1, stdout);
  block_la_print(b2, stdout);
  block_la_calculate_velocities(b2);
  block_print_velocity_profile(b2, stdout);
  block_la_print(b3, stdout);
  block_la_calculate_velocities(b3);
  block_print_velocity_profile(b3, stdout);



  machine_free(machine);
  return 0;
}

#endif