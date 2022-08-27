#include "../defines.h"
#include "../machine.h"
#include "../program_la.h"
#include "../block_la.h"
#include "../point.h"
#include "../fsm.h"

#define eprintf(...) fprintf(stderr, __VA_ARGS__)

#if 1
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
  char *velocity_profile_save;
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

  velocity_profile_save = argv[2];
  program_la_look_ahead(p, machine, velocity_profile_save);

  tt = 0;
  printf("n,t_tot,t_blk,lambda,s,feed,x,y,z\n");
  while ((b = program_la_next(p))) {
    if (block_la_type(b) == RAPID || block_la_type(b) > ARC_CCW) {
      continue;
    }
    eprintf("Interpolating the block %s\n", block_la_line(b));
    // interpolation loop
    // careful: we check t <= block_dt(b) + tq/2.0 for double values are
    // never exact, and we may have that adding many tq carries over a small
    // error that accumuates and may result in n*tb being greater than Dt
    // (if so, we would miss the last step
    for (t = 0; t <= block_la_dt(b) + tq; t += tq, tt += tq) {
      lambda = block_la_lambda(b, t, &f);
      sp = block_la_interpolate(b, lambda);
      if (!sp) continue;
      printf("%lu,%f,%f,%f,%f,%f,%f,%f,%f\n", block_la_n(b), tt, t,
        lambda, lambda * block_la_length(b), f,
        point_x(sp), point_y(sp), point_z(sp));
      wait_next(5e6);
    }
  }

  machine_free(machine);
  return 0;
}

#endif